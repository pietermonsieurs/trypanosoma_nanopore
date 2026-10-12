#!/bin/bash -l

#SBATCH --job-name=map_porec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=05:00:00
#SBATCH --account=ap_itg_tryps

## modules
module load minimap2
module load BEDTools
module load SAMtools

## settings
threads=28

## input parameters
# assembly_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/purge_dups/purged.fa
# porec_reads=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_porec/poreC_combined.fastq.gz
# bam_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/purge_dups/purged.pored.bam
# bed_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/purge_dups/purged.pored.bed

## input parameter after updating the assembly file with contigs and scaffolds
## obtained from redoing Hifiasm and Flye on nanopore sequencing data that were not
## mapping to the reference genome
assembly_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_for_PoreC/MM05_final_confident.fasta
porec_reads=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_porec/poreC_combined.fastq.gz
bam_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_for_PoreC/MM05_final_confident.pored.bam
bed_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_for_PoreC/MM05_final_confident.pored.bed



## index the assembly file
samtools faidx ${assembly_file}

## run minimap2 to map the poreC reads to the assembly and convert to a sorted BAM file
minimap2 -ax map-ont --secondary=no -t ${threads} ${assembly_file} ${porec_reads} \
    | samtools view -bh -F 0x104 -@ ${threads} \
    | samtools sort -n -@ ${threads} -o ${bam_file}

bamToBed -i ${bam_file} > ${bed_file}
