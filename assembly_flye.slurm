#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=64 --nodes=1 --partition=zen3
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_mpu

## load modules
module load Flye
module load SAMtools

## parameters
out_dir_base=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_update2/
fastq_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_update/
threads=64

fastq_file=${fastq_dir}/barcode${barcode}.fastq.gz
out_dir=${out_dir_base}/barcode${barcode}

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

## loop for all barcode
# for barcode in {01..06}; do sbatch --export=barcode=${barcode} /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_flye.slurm; done
# sbatch --export=barcode=02 /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_flye.slurm
# sbatch --export=barcode=05 /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_flye.slurm
# sbatch --export=barcode=04 /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_flye.slurm



