#!/usr/bin/env python3

## do selection of only the sequences originating from the nuclear genome, and omit
## the smaller contigs (potnentially minichromosomes) that are not part of the nuclear genome
## Nulcear chromosomes start with "Tb927" and are longer than 100,000 bp

# module load Python
# module load BioPython

from Bio import SeqIO

## parameter settings
fasta_file = "/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome.fasta"
output_file = "/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome_nuclear.fasta"

## read the fasta file and filter sequences
with open(fasta_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id.startswith("Tb927") and len(record.seq) > 100000:
            SeqIO.write(record, outfile, "fasta")

