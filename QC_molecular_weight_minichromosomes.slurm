#!/usr/bin/bash

#SBATCH --job-name=minichroms
#SBATCH --account ap_itg_tryps
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=28

module load minimap2
module load SAMtools


## paraameter settings
ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_sofi_10x/data/refgenome/glasgow/tco_il3000_v0.47.fasta
threads=2

## strain MSOROM7
# fastq_file=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/barcode05.fastq.gz  ## MSOROM7
# sample_name=MSOROM7

## strain IL3000
fastq_file=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/barcode06.fastq.gz  ## IL3000
sample_name=IL3000

output_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/QC/molecular_weight/

## run minimap2 and indexing + idxstats
minimap2 \
    -ax map-ont \
    $ref_genome \
    $fastq_file \
    -t ${threads} | \
    samtools sort -@$threads -o ${output_dir}/${sample_name}.bam

samtools index ${output_dir}/${sample_name}.bam
samtools idxstats ${output_dir}/${sample_name}.bam > ${output_dir}/${sample_name}_idxstats.txt