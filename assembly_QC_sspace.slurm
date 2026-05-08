#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_tryps

conda activate sspace_longread

export PATH=/user/antwerpen/205/vsc20587/data/software/sspace_longread:${PATH}

assembly_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/consensus_TbgI_merged.fasta
fastq_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/sspace/Tbg_merged_filtlong.fastq
output_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/sspace/TbgI_merged_sspace

SSPACE-LongRead.pl  \
    -c ${assembly_file} \
    -p ${fastq_file} \
    -b ${output_dir} \
    -t 28 
    