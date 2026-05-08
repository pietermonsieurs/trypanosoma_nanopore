#!/usr/bin/bash

#SBATCH --job-name=inspector_manual_assembly
#SBATCH --account ap_itg_tryps
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=28

## this can be run using the wrapper created in the file assembly_QC_manual_assembly_inspector_wrapper.sh, 
## which will loop through all the fastq files in the directory and run this script for each of them, 
## creating a separate output directory for each sample

module load calcua/2024a
module load minimap2
module load SAMtools
module load Pysam
module load statsmodels


export PATH=/user/antwerpen/205/vsc20587/data/software/Inspector/:${PATH}

threads=28
genome_assembly_fasta=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/manual_assembly.fasta

## hard coded example
# fastq_file_nanopore=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/barcode03.fastq.gz
# sample_name=ELIANE

## create variable for output directory and create it if not exists 
output_directory=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/inspector/${sample_name}/
mkdir -p ${output_directory}

## run Inspector
inspector.py \
    -c ${genome_assembly_fasta} \
    -r ${fastq_file_nanopore} \
    -d nanopore \
    -t ${threads} \
    -o ${output_directory}
