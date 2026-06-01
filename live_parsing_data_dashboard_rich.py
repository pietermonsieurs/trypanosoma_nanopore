#!/usr/bin/env python3

import gzip
import os
import time
import numpy as np

from collections import defaultdict
from Bio import SeqIO

from rich.live import Live
from rich.table import Table


# ============================================================
# CONFIG
# ============================================================

DATA_DIR = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run/fastq_pass/"

SAMPLES = (
    "barcode09",
    "barcode10",
    "barcode11",
    "barcode12",
    "barcode13"
)


# ============================================================
# STORAGE
# ============================================================

processed_files = defaultdict(list)
processed_reads = defaultdict(set)

total_base_pairs = defaultdict(int)
read_lengths = defaultdict(list)

read_count = defaultdict(int)
quality_score_sum = defaultdict(int)
quality_base_count = defaultdict(int)
reads_below_1kb = defaultdict(int)
reads_above_10kb = defaultdict(int)
reads_above_50kb = defaultdict(int)


# ============================================================
# N50
# ============================================================

def calculate_n50(lengths):

    if not lengths:
        return 0

    lengths = sorted(lengths, reverse=True)

    total = sum(lengths)

    cumulative = 0

    for length in lengths:

        cumulative += length

        if cumulative >= total / 2:
            return length

    return 0


# ============================================================
# HELPERS
# ============================================================

def calculate_mean_qscore(sample):

    if quality_base_count[sample] == 0:
        return 0.0

    return quality_score_sum[sample] / quality_base_count[sample]


def calculate_error_rate(mean_qscore):

    return 10 ** (-mean_qscore / 10) if mean_qscore > 0 else 1.0


def format_percentage(numerator, denominator):

    if denominator == 0:
        return "0.0%"

    return f"{(numerator / denominator) * 100:.1f}%"


def format_error_rate(mean_qscore):

    error_rate_percent = calculate_error_rate(mean_qscore) * 100

    if error_rate_percent >= 0.1:
        return f"{error_rate_percent:.2f}%"

    return f"{error_rate_percent:.3g}%"



# ============================================================
# DASHBOARD
# ============================================================

def build_dashboard():

    table = Table(title="Nanopore Live Monitor")

    table.add_column("Sample")
    table.add_column("Total BP")
    table.add_column("% Dataset")
    table.add_column("Reads")
    table.add_column("Median")
    table.add_column("N50")
    table.add_column("Mean Q")
    table.add_column("Mean Error")
    table.add_column("<1 kb")
    table.add_column(">10 kb")
    table.add_column(">50 kb")

    total_dataset_bp = sum(total_base_pairs.values())

    for sample in SAMPLES:

        lengths = read_lengths[sample]

        median = int(np.median(lengths)) if lengths else 0

        n50 = calculate_n50(lengths)
        mean_qscore = calculate_mean_qscore(sample)

        table.add_row(
            sample,
            f"{total_base_pairs[sample]:,}",
            format_percentage(total_base_pairs[sample], total_dataset_bp),
            f"{read_count[sample]:,}",
            str(median),
            str(n50),
            f"{mean_qscore:.2f}",
            format_error_rate(mean_qscore),
            format_percentage(reads_below_1kb[sample], read_count[sample]),
            format_percentage(reads_above_10kb[sample], read_count[sample]),
            format_percentage(reads_above_50kb[sample], read_count[sample])
        )

    return table


# ============================================================
# PARSER
# ============================================================

def parse_live_data():

    with Live(
        build_dashboard(),
        refresh_per_second=1,
        screen=True
    ) as live:

        while True:

            for sample_dir in os.listdir(DATA_DIR):

                if sample_dir not in SAMPLES:
                    continue

                sample_path = os.path.join(
                    DATA_DIR,
                    sample_dir
                )

                if not os.path.isdir(sample_path):
                    continue

                files = [
                    f for f in os.listdir(sample_path)
                    if f.endswith(".fastq.gz")
                ]

                if not files:
                    continue

                new_files = [
                    f for f in files
                    if f not in processed_files[sample_dir]
                ]

                for fastq_file in sorted(new_files):

                    processed_files[sample_dir].append(
                        fastq_file
                    )

                    print(
                        f"Processing new file: "
                        f"{sample_dir} -> {fastq_file}"
                    )

                    fastq_path = os.path.join(
                        sample_path,
                        fastq_file
                    )

                    with gzip.open(fastq_path, "rt") as handle:

                        for record in SeqIO.parse(handle, "fastq"):

                            read_id = record.id

                            if read_id in processed_reads[sample_dir]:
                                continue

                            processed_reads[sample_dir].add(
                                read_id
                            )

                            seq = str(record.seq)
                            qualities = record.letter_annotations.get(
                                "phred_quality",
                                []
                            )

                            read_len = len(seq)

                            total_base_pairs[sample_dir] += read_len

                            if qualities:
                                quality_score_sum[sample_dir] += sum(qualities)
                                quality_base_count[sample_dir] += len(qualities)

                            read_lengths[sample_dir].append(
                                read_len
                            )

                            read_count[sample_dir] += 1

                            if read_len < 1000:
                                reads_below_1kb[sample_dir] += 1

                            if read_len > 10000:
                                reads_above_10kb[sample_dir] += 1

                            if read_len > 50000:
                                reads_above_50kb[sample_dir] += 1

                            live.update(
                                build_dashboard()
                            )

            time.sleep(5)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    parse_live_data()