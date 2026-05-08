#!/usr/bin/env python3

import os
import sys
import Bio.SeqIO

## Step 1: link the ids of the fasta file in the initial input file 
## with the generic ids that sspace gives when copying the fasta file
## to an intermediate folder
my_debug = False
source_file = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/consensus_TbgI_merged.fasta'
sspace_fasta_file = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/sspace/TbgI_merged_sspace/intermediate_files/contigs.fa'
sspace_output_file = '/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/sspace/TbgI_merged_sspace/scaffold_evidence.txt'

## extract IDs - there are in the same order so can be combined afterwards
## in a dictionary
source_ids = []
for record in Bio.SeqIO.parse(source_file, 'fasta'):
    my_debug and print(record.id)  
    source_ids.append(record.id)

sspace_ids = []
for record in Bio.SeqIO.parse(sspace_fasta_file, 'fasta'):
    sspace_id = record.id.split("_")[0]
    my_debug and print(sspace_id)  
    sspace_ids.append(sspace_id)

## make dictionary with the generic sspace ids as keys and the original ids as values
id_mapping = {}
for sspace_id, source_id in zip(sspace_ids, source_ids):
    id_mapping[sspace_id] = source_id

## Step 2: parse the sspace output files and replace the generic ids with the original ids. Output 
## of sspace looks like the text below, where every scaffold starts with a new ">" and then the 
## combination of contigs is given that make up the scaffold. Per contig is given the contig ID
## (only number, so contig_16 -> 16) and is preceeded with a f or a r to indicated  whether it is 
## a forward or reverse contig. After the contig ID, the size of the contig is given, and then the
## links and gaps that are present in the scaffold. The
# the size of the contig, and the links and gaps that are present in the scaffold. The
# >scaffold1|size4100888|tigs1|Lineair
# f16|size4100888

# >scaffold2|size3656834|tigs3|Lineair
# f36|size2898242|links1|gaps4109
# f3|size15752|links2|gaps2052
# r25|size736679

## initiate the scaffold information variables outside the loop, so they can be used to store 
## the information of the current scaffold while processing the contig lines
scaffold_id = None
scaffold_size = None
scaffold_number_of_contigs = None

for line in open(sspace_output_file):
    line = line.strip()

    if line == '':
        # skip empty lines
        continue
    elif line.startswith('>'):
        # this is a new scaffold, so we can just print it as is
        parts = line.split('|')
        scaffold_id = parts[0].replace(">", "")  # e.g. >scaffold1
        scaffold_size = parts[1].replace("size", "")  # e.g. size4100888
        scaffold_number_of_contigs = parts[2].replace("tigs", "")  
    else:
        # this is a contig line, so we need to replace the generic id with the original id
        parts = line.split('|')
        contig_info = parts[0]  # e.g. f16 or r25
        contig_id = f"contig{contig_info[1:]}"  # remove the f or r to get the contig number, and add contig_ prefix to match the original id format
        my_debug and print(f"Processing contig ID: {contig_id}")  # debug print

        contig_direction = contig_info[0]  # f or r
        contig_size = parts[1].replace("size", "")  # e.g. size4100888
        links_and_gaps = parts[2:]  # e.g. links1|gaps4109

        # get the original id from the mapping dictionary
        original_id = id_mapping[contig_id]
        my_debug and print(original_id)  # debug print

        ## convert r and f to reverse and forward, to make it more clear
        if contig_direction == 'f':
            contig_direction = 'forward'
        elif contig_direction == 'r':
            contig_direction = 'reverse'

        # reconstruct the line with the original id
        new_line = f"{scaffold_id}\t{scaffold_size}\t{scaffold_number_of_contigs}\t{contig_id}\t{original_id}\t{contig_direction}\t{contig_size}\t"
        if links_and_gaps:
            new_line += '\t'.join(links_and_gaps)

        print(new_line)
