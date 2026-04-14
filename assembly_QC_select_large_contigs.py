#!/usr/bin/env python3

## script to select large contigs from assembly QC output. Based on the 
## assembly fasta file only select those contigs that are larger than a 
# specified length. This is useful to create the backbone of the reference
## genome

from Bio import SeqIO

## parameter settings
fasta_input_dir = '/user/antwerpen/205/vsc20587/aitg_scratch/jvdabbeele/Tbg_de_novo/'
fasta_ouput_dir = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/mauve/'

## set the minimum contig length
min_contig_length = 50000

## define the input fasta file
input_fasta = 'consensus_TbgI_merged.fasta'

## read the fasta file and select the large contigs
for record in SeqIO.parse(f"{fasta_input_dir}{input_fasta}.fasta", 'fasta'):
    if len(record.seq) >= min_contig_length:
        SeqIO.write(record, f"{fasta_output_dir}{input_fasta}.minL'{min_contig_length}.fasta", 'fasta')

