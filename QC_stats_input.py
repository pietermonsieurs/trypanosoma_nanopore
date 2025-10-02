#!/usr/bin/env python

import gzip
import sys
import statistics
import glob
import re



def fastq_reader(handle):
    while True:
        header = handle.readline().decode().strip()
        if not header:
            break
        seq = handle.readline().decode().strip()
        handle.readline()  # skip '+'
        handle.readline()  # skip qualities
        yield seq

def summarize_fastq(fastq_file):
    lengths = []
    n_reads = 0

    # open .gz with gzip, normal with open()
    open_func = gzip.open if str(fastq_file).endswith(".gz") else open
    with open_func(fastq_file, "rb") as f:
        for seq in fastq_reader(f):
            n_reads += 1
            lengths.append(len(seq))

    mean_len = statistics.mean(lengths) if lengths else 0
    sd_len = statistics.stdev(lengths) if len(lengths) > 1 else 0
    total_size = sum(lengths)

    # print(f"File: {fastq_file}")
    # print(f"  Total reads: {n_reads}")
    # print(f"  Mean read length: {mean_len:.2f}")
    # print(f"  StdDev read length: {sd_len:.2f}")
    # print(f"  Total bases: {total_size}\n")

    return(n_reads, mean_len, sd_len, total_size)


for i in range(1, 7):
    barcode = f"{i:02d}"
    fastq_dir = f'/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass/barcode{barcode}/'
    fastq_files = glob.glob(f"{fastq_dir}/*.fastq.gz")   # adjust path
    for fq_file in fastq_files:
        n_reads, mean_len, sd_len, total_size = summarize_fastq(fq_file)

        ## extract the hour
        match = re.search(r'_(\d+)\.fastq\.gz$', fq_file)
        if match:
            hour = int(match.group(1))   # int() removes leading zeros
            # print(hour)   # 79

        print(f"barcode{barcode},{hour},{n_reads},{mean_len},{sd_len},{total_size}")

