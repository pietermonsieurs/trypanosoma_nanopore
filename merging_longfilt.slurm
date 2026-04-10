#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_mpu

# module load Filtlong
export PATH=/user/antwerpen/205/vsc20587/data/software/Filtlong/bin/:${PATH}

data_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_merging

## setting for Tbg
# input_file=${data_dir}/Tbg_merged.fastq.gz
# output_file=${data_dir}/Tbg_merged_filtlong.fastq.gz

## setting for Tbg type I- so removal of ABBA (TbgII strain)
input_file=${data_dir}/TbgI_merged.fastq.gz
output_file=${data_dir}/TbgI_merged_filtlong.fastq.gz


## setting for Tcongo
# input_file=${data_dir}/Tcongo_merged.fastq.gz
# output_file=${data_dir}/Tcongo_merged_filtlong.fastq.gz

## setting for Tcongo MSOROM7
# input_file=${data_dir}/Tcongo_MSOROM7_merged.fastq.gz
# output_file=${data_dir}/Tcongo_MSOROM7_merged_filtlong.fastq.gz


filtlong \
    --min_length 1kb \
    --keep_percent 90 \
    --target_bases 1.75gb \
    ${input_file} \
    | gzip > ${output_file}
