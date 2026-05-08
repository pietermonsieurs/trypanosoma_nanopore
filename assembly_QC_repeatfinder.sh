#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu

## export the PATH with tandem repeat finder
export PATH=/user/antwerpen/205/vsc20587/data/software/tandemrepeatfinder:${PATH}

## parameter settings 
fasta_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/consensus_TbgI_merged.fasta
output_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/trf/

cd ${output_dir}
trf409.linux64 ${fasta_file} 2 5 5 80 15 40 100 -d -h


## perform a selection on the size of the repeat (total size of the repeated region)
~/scratch/trypanosoma_nanopore/bin/assembly_QC_repeatfinder_parse.py consensus_TbgI_merged.fasta.2.5.5.80.15.40.100.dat consensus_TbgI_merged.fasta.2.5.5.80.15.40.100.parsed.dat --minsize 10000



~/scratch/trypanosoma_nanopore/bin/assembly_QC_repeatfinder_parse.py consensus_TbgI_merged.fasta.2.5.5.80.15.40.100.dat consensus_TbgI_merged.fasta.2.5.5.80.15.40.100.parsed.dat --minsize 10000
