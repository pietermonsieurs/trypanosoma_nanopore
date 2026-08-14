#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=64 --nodes=1 --partition=zen3
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_mpu

## load modules
module load Flye
module load SAMtools

## parameters
out_dir_base=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/
fastq_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/
threads=64

## create output directory
mkdir -p ${out_dir_base}

## input and output files
fastq_file=${fastq_dir}/mm05_combined_filtlong.fastq.gz
out_dir=${out_dir_base}/mm05_flye_assembly_filtlong

flye \
    --asm-coverage 50 \
    --nano-hq ${fastq_file} \
    --out-dir ${out_dir} \
    --threads ${threads} \
    --read-error 0.05 \
    --scaffold \
    --debug \
    --iterations 4 \
    --genome-size 35m





