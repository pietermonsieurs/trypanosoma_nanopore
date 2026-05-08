#!/usr/bin/bash

#SBATCH --job-name=assembly_QC_manual_assembly_quast
#SBATCH --account ap_itg_tryps
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=28

## load modules and adapt PATH
# module load Python/3.14.2-GCCcore-15.2.0
# PATH=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/:${PATH}
# srun -n 1 -c 16 -A ap_itg_tryps --interactive --pty bash 
# module load Miniconda3

source ~/.bashrc
conda activate sspace

## check whether variables are read in correctly
echo $sample_name
echo $nanopore_reads
echo $illumina_R1_reads
echo $illumina_R2_reads

## parameter settings
assembly_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/
assembly_file=${assembly_dir}/manual_assembly.fasta
reference_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome_nuclear.fasta
reference_gff=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927.gff

## some hard coded examples, in case you want to run this for one sample. Otherwise use the wrapper
## script assembly_QC_manual_assembly_quast_wrapper.sh
# nanopore_reads=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_merging/TbgI_merged.fastq.gz
# illumina_R1_reads=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Trypanosoma/rawreads/BATCH1/ABBA/DP8400009737BL_L01_SP2004030345_1.fq.gz
# illumina_R2_reads=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Trypanosoma/rawreads/BATCH1/ABBA/DP8400009737BL_L01_SP2004030345_2.fq.gz
# output_dir=${assembly_dir}/quast_manual_assembly/

## define output directory based on sample name
output_dir=${assembly_dir}/quast/manual_assembly_${sample_name}/

## create output directory if it doesn't exist
mkdir -p ${output_dir}

## run QUAST
quast.py ${assembly_file} \
    -o ${output_dir} \
    -t 28 \
    -r ${reference_genome} \
    -g ${reference_gff} \
    --nanopore ${nanopore_reads} \
    --pe1 ${illumina_R1_reads} \
    --pe2 ${illumina_R2_reads} \
    --eukaryote \
    --min-contig 1000 




