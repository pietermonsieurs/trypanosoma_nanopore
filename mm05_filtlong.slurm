#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_mpu

# module load Filtlong
export PATH=/user/antwerpen/205/vsc20587/data/software/Filtlong/bin/:${PATH}

data_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/

## setting for Tbg
# input_file=${data_dir}/Tbg_merged.fastq.gz
# output_file=${data_dir}/Tbg_merged_filtlong.fastq.gz

## setting for MM05 strain
input_file=${data_dir}/mm05_combined.fastq.gz
output_file=${data_dir}/mm05_combined_filtlong.fastq.gz

filtlong \
    --min_length 1kb \
    --keep_percent 90 \
    --target_bases 1.75gb \
    ${input_file} \
    | gzip > ${output_file}