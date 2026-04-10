#!/usr/bin/env python

import os

assembly_dir = "/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly/"
print_to_file = 1  # set to 1 to print in CSV format

## print header if printing to file
if print_to_file == 1:
    print("Barcode,Num_Contigs,Total_Length,Mean_Length,Median_Length,N50")
    
for barcode in range(1,7):
    barcode_str = f"{barcode:02d}"
    assembly_file = f"{assembly_dir}/barcode{barcode_str}/assembly.fasta"
    
    if not os.path.exists(assembly_file):
        not print_to_file and print(f"File not found: {assembly_file}, skipping...")
        continue

    with open(assembly_file, "r") as f:
        total_length = 0
        n_contigs = 0
        lengths = []
        
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    seq_len = len(seq)
                    lengths.append(seq_len)
                    total_length += seq_len
                    n_contigs += 1
                    seq = ""
            else:
                seq += line
        if seq:
            seq_len = len(seq)
            lengths.append(seq_len)
            total_length += seq_len
            n_contigs += 1
        
        if lengths:
            mean_length = total_length / n_contigs
            lengths.sort()
            median_length = lengths[n_contigs // 2] if n_contigs % 2 == 1 else (lengths[n_contigs // 2 - 1] + lengths[n_contigs // 2]) / 2
            n50 = lengths[0]
            cumsum = 0
            for length in reversed(lengths):
                cumsum += length
                if cumsum >= total_length / 2:
                    n50 = length
                    break
        else:
            mean_length = 0
            median_length = 0
            n50 = 0
        
        if print_to_file == 0:
            print(f"Barcode: {barcode_str}")
            print(f"Number of contigs: {n_contigs}")
            print(f"Total assembly length: {total_length}")
            print(f"Mean contig length: {mean_length:.2f}")
            print(f"Median contig length: {median_length}")
            print(f"N50: {n50}")
            print("-" * 40)
        else:
            print(f"{barcode_str},{n_contigs},{total_length},{mean_length:.2f},{median_length},{n50}")