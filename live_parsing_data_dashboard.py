#!/usr/bin/env python3

import gzip
import os
from collections import defaultdict

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import SeqIO


# DATA_DIR = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run/fastq_pass/"
DATA_DIR = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run_test/"
SAMPLES = ("barcode09", "barcode10", "barcode11", "barcode12")


def calculate_n50(lengths):
    if not lengths:
        return 0

    lengths_sorted = sorted(lengths, reverse=True)
    half = sum(lengths_sorted) / 2

    cumulative = 0
    for length in lengths_sorted:
        cumulative += length
        if cumulative >= half:
            return length

    return 0


@st.cache_resource
def get_stats_store():
    return {
        "processed_reads": defaultdict(set),
        "total_bp": defaultdict(int),
        "read_lengths": defaultdict(list),
        "latest_read_length": defaultdict(int),
        "read_count": defaultdict(int),
        "latest_sequence": defaultdict(str),
        "file_mtimes": defaultdict(dict),
    }


def parse_updated_files(data_dir, samples, stats):

    for sample in samples:
        sample_path = os.path.join(data_dir, sample)
        if not os.path.isdir(sample_path):
            continue

        fastq_files = sorted(
            file_name
            for file_name in os.listdir(sample_path)
            if file_name.endswith(".fastq.gz")
        )

        for fastq_name in fastq_files:
            fastq_path = os.path.join(sample_path, fastq_name)

            try:
                current_mtime = os.path.getmtime(fastq_path)
            except OSError:
                continue

            previous_mtime = stats["file_mtimes"][sample].get(fastq_name)
            if previous_mtime is not None and current_mtime <= previous_mtime:
                continue

            try:
                with gzip.open(fastq_path, "rt") as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        read_id = record.id
                        if read_id in stats["processed_reads"][sample]:
                            continue

                        read_len = len(record.seq)
                        stats["processed_reads"][sample].add(read_id)
                        stats["latest_read_length"][sample] = read_len
                        stats["total_bp"][sample] += read_len
                        stats["read_lengths"][sample].append(read_len)
                        stats["read_count"][sample] += 1
                        stats["latest_sequence"][sample] = str(record.seq)[:1000]
            except Exception:
                continue

            stats["file_mtimes"][sample][fastq_name] = current_mtime


def sample_histogram(lengths, sample_name):
    if not lengths:
        st.caption("Waiting for reads...")
        return

    df = pd.DataFrame({"Read Length": lengths})
    fig = px.histogram(
        df,
        x="Read Length",
        nbins=40,
        title=None,
    )
    fig.update_layout(height=160, margin=dict(l=0, r=0, t=0, b=0))
    fig.update_xaxes(title_text=None)
    fig.update_yaxes(title_text=None)
    st.plotly_chart(fig, use_container_width=True)


st.set_page_config(layout="wide", page_title="Nanopore 4-Sample Live Dashboard")
st.markdown(
    "<h4 style='margin-bottom: 0.3rem;'>Nanopore Live Dashboard (4 Samples)</h4>",
    unsafe_allow_html=True,
)

refresh_seconds = 10
st.markdown(f"<meta http-equiv='refresh' content='{refresh_seconds}'>", unsafe_allow_html=True)
# st.caption("Auto-refresh: 10s. Processed reads/files are kept in cached memory between reruns.")

stats = get_stats_store()
parse_updated_files(DATA_DIR, SAMPLES, stats)

header_cols = st.columns([1, 1, 1, 1, 1])
with header_cols[0]:
    st.markdown("**Parameter**")
for index, sample in enumerate(SAMPLES, start=1):
    with header_cols[index]:
        st.markdown(f"**{sample}**")

rows = [
    ("Total Reads", lambda sample: f"{stats['read_count'][sample]:,}"),
    ("Latest Length", lambda sample: f"{stats['latest_read_length'][sample]:,}"),
    ("N50", lambda sample: f"{calculate_n50(stats['read_lengths'][sample]):,}"),
    ("Total BP", lambda sample: f"{stats['total_bp'][sample]:,}"),
    ("First 1000 bp", lambda sample: stats["latest_sequence"][sample][:1000] if stats["latest_sequence"][sample] else "N/A"),
]

for row_label, value_fn in rows:
    row_cols = st.columns([1, 1, 1, 1, 1])
    with row_cols[0]:
        st.markdown(f"**{row_label}**")
    for index, sample in enumerate(SAMPLES, start=1):
        with row_cols[index]:
            # st.metric(label="", value=value_fn(sample), label_visibility="collapsed")
            st.write(value_fn(sample))

hist_cols = st.columns([1, 1, 1, 1, 1])
with hist_cols[0]:
    st.markdown("**Read Length Distribution**")
for index, sample in enumerate(SAMPLES, start=1):
    with hist_cols[index]:
        sample_histogram(stats["read_lengths"][sample], sample)