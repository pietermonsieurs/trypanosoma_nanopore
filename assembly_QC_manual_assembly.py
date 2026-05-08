#!/usr/bin/env python3

from Bio import SeqIO

## select contigs and scaffolds from three different assemblies, and 
## then merge them together to get a better assembly
data_dir = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/'
tbgI_file = f'{data_dir}/consensus_TbgI_merged.fasta'
sspace_file = f'{data_dir}/sspace/TbgI_merged_sspace/scaffolds.fasta'
gapless_file = f'{data_dir}/gapless/gapless.fa'

## read in the three assemblies into a seqIO object
tbgI_assembly = SeqIO.to_dict(SeqIO.parse(tbgI_file, "fasta"))
sspace_assembly = SeqIO.to_dict(SeqIO.parse(sspace_file, "fasta"))
gapless_assembly = SeqIO.to_dict(SeqIO.parse(gapless_file, "fasta"))

## select chromosome by chromome manually by selecting one of the contigs out of
## the Bio::Seq object
chr1 = gapless_assembly['scaffold67'].seq
chr2 = f"{tbgI_assembly['contig_34'].seq}{tbgI_assembly['contig_61'].seq}{tbgI_assembly['scaffold_53'].reverse_complement().seq}"
chr3 = sspace_assembly['scaffold10|size1896602'].seq
chr4 = tbgI_assembly['contig_76'].seq
chr5 = sspace_assembly['scaffold8|size1681482'].seq
chr6 = tbgI_assembly['scaffold_63'].seq
chr7 = sspace_assembly['scaffold4|size3255367'].seq
chr8 = sspace_assembly['scaffold7|size2874005'].seq
chr9 = sspace_assembly['scaffold2|size3656834'].seq
chr10 = gapless_assembly['scaffold53'].seq
chr11 = sspace_assembly['scaffold3|size6165291'].seq

## wrtie the merged assembly to a fasta file, with the chromosome names as headers 
output_file = f"{data_dir}/merged_assembly.fasta"
output_fasta = f">chr1\n{chr1}\n>chr2\n{chr2}\n>chr3\n{chr3}\n>chr4\n{chr4}\n>chr5\n{chr5}\n>chr6\n{chr6}\n>chr7\n{chr7}\n>chr8\n{chr8}\n>chr9\n{chr9}\n>chr10\n{chr10}\n>chr11\n{chr11}"
with open(output_file, 'w') as f:
    f.write(output_fasta)
f.close()


