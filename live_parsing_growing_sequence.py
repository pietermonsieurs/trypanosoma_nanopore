#!/usr/bin/env python3

import gzip
import os
import time
import re

from collections import deque
from collections import defaultdict
from Bio import SeqIO

from rich.live import Live
from rich.table import Table
from rich.text import Text


# ============================================================
# CONFIG
# ============================================================

DATA_DIR = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run_test/"

SAMPLES = (
    "barcode09",
    "barcode10",
    "barcode11",
    "barcode12",
    "barcode13"
)

POLL_SECONDS = 0.05
GROWTH_STEP = 1
GROWTH_DELAY_SECONDS = 0.01 
TAIL_NT = 100
DISPLAY_FRACTION = 0.10
KEEP_COMPLETED = 10
MAX_GROWING = 20

SANGER_BASE_COLORS = {
    "A": "green",
    "C": "blue",
    "G": "white",
    "T": "red",
    "U": "red",
    "N": "yellow",
}

ADAPTED_FASTQ_PATTERN = re.compile(r"adapt|adapted|trim", re.IGNORECASE)


def colorize_sequence(sequence):
    colored = Text()

    for base in sequence:
        colored.append(base, style=SANGER_BASE_COLORS.get(base.upper(), "magenta"))

    return colored


def select_latest_fastq(sample_path):
    fastq_files = [
        file_name
        for file_name in os.listdir(sample_path)
        if file_name.endswith(".fastq.gz")
    ]

    if not fastq_files:
        return None

    adapted_fastq_files = [
        file_name for file_name in fastq_files
        if ADAPTED_FASTQ_PATTERN.search(file_name)
    ]

    candidate_files = adapted_fastq_files or fastq_files

    try:
        return max(
            candidate_files,
            key=lambda file_name: os.path.getmtime(os.path.join(sample_path, file_name))
        )
    except OSError:
        return None


def build_table(last_reads):
    table = Table(title="Nanopore Live Sequence Growth (Top: Last 5 Completed, Then Active)")
    table.add_column("Growing Sequence (10% view, last 100 nt)", overflow="fold")
    table.add_column("Length", justify="right")
    table.add_column("Barcode")

    for row in last_reads:
        table.add_row(
            colorize_sequence(row.get("growing_sequence", "")),
            row.get("length", "..."),
            row.get("barcode", "..."),
        )

    return table


def reorder_and_prune_visible_reads(visible_reads):
    completed = [row for row in visible_reads if row["is_done"]]
    active = [row for row in visible_reads if not row["is_done"]]

    completed.sort(key=lambda row: row.get("completed_at", 0), reverse=True)
    completed = completed[:KEEP_COMPLETED]

    visible_reads[:] = completed + active


def discover_new_reads(file_mtimes, processed_reads, pending_reads):
    for barcode in SAMPLES:
        sample_path = os.path.join(DATA_DIR, barcode)

        if not os.path.isdir(sample_path):
            continue

        fastq_file = select_latest_fastq(sample_path)

        if fastq_file is None:
            continue

        fastq_path = os.path.join(sample_path, fastq_file)

        try:
            current_mtime = os.path.getmtime(fastq_path)
        except OSError:
            continue

        previous_mtime = file_mtimes[barcode].get(fastq_file)
        if previous_mtime is not None and current_mtime <= previous_mtime:
            continue

        with gzip.open(fastq_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_id = record.id

                if read_id in processed_reads[barcode]:
                    continue

                processed_reads[barcode].add(read_id)

                sequence = str(record.seq)
                pending_reads.append(
                    {
                        "sequence": sequence,
                        "progress": 0,
                        "growing_sequence": "",
                        "display_sequence": "",
                        "length": "...",
                        "barcode": "...",
                        "true_barcode": barcode,
                        "is_done": False,
                        "completed_at": None,
                    }
                )

        file_mtimes[barcode].clear()
        file_mtimes[barcode][fastq_file] = current_mtime


def refill_visible_reads(visible_reads, pending_reads):
    active_count = sum(1 for row in visible_reads if not row["is_done"])

    while pending_reads and active_count < MAX_GROWING:
        visible_reads.append(pending_reads.popleft())
        active_count += 1


def parse_live_data():
    file_mtimes = defaultdict(dict)
    processed_reads = defaultdict(set)
    pending_reads = deque()
    visible_reads = []
    last_poll = 0.0

    with Live(build_table(visible_reads), refresh_per_second=20, screen=True) as live:
        while True:
            now = time.time()
            if (now - last_poll) >= POLL_SECONDS:
                discover_new_reads(file_mtimes, processed_reads, pending_reads)
                last_poll = now

            refill_visible_reads(visible_reads, pending_reads)

            for row in visible_reads:
                if row["is_done"]:
                    continue

                seq = row["sequence"]
                seq_len = len(seq)
                display_target_len = max(1, int(seq_len * DISPLAY_FRACTION))
                row["progress"] = min(display_target_len, row["progress"] + GROWTH_STEP)

                visible_end = row["progress"]
                visible_start = max(0, visible_end - TAIL_NT)
                row["display_sequence"] = seq[visible_start:visible_end]

                row["growing_sequence"] = row["display_sequence"]

                if row["progress"] >= display_target_len:
                    row["length"] = str(seq_len)
                    row["barcode"] = row["true_barcode"]
                    row["is_done"] = True
                    row["completed_at"] = time.time()

            # Keep only newest completed rows in the dashboard while active rows continue growing.
            reorder_and_prune_visible_reads(visible_reads)
            refill_visible_reads(visible_reads, pending_reads)

            live.update(build_table(visible_reads), refresh=True)
            time.sleep(GROWTH_DELAY_SECONDS)


if __name__ == "__main__":
    parse_live_data()

