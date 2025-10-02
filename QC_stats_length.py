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

def get_lengths(fastq_file, barcode):
    open_func = gzip.open if str(fastq_file).endswith(".gz") else open

    with open_func(fastq_file, "rb") as f:
        for seq in fastq_reader(f):
            print(f"{barcode},{len(seq)}")

data_dir = "/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/"
for i in range(1, 7):
    barcode = f"{i:02d}"
    fastq_file = f"{data_dir}/barcode{barcode}.fastq.gz"
    get_lengths(fastq_file, f"barcode{barcode}")

    