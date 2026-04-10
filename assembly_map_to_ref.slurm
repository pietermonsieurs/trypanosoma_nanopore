#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=20:00:00
#SBATCH -A ap_itg_tryps

module load minimap2

## specify the settings
# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome_nuclear.fasta
# ref_name=TREU927

ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-68_TbruceigambienseDAL972_Genome.fasta
ref_name=DAL972

# assembly_genome=/user/antwerpen/205/vsc20587/aitg_scratch/jvdabbeele/Tbg_de_novo/consensus_bc03.fasta
results_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_QC/

## extract the assembly name
assembly_name=$(basename $assembly_genome .fasta)

## run the mapping
minimap2 -x asm5 \
    --secondary=no \
    $ref_genome \
    $assembly_genome \
    > ${assembly_name}_vs_${ref_name}.paf

## do a filtering of the paf file to keep only regions bigger than 10,000 bp and 
## with a mapping quality bigger than 30
awk '$10 >= 10000 && $12 >= 30' ${assembly_name}_vs_${ref_name}.paf > ${assembly_name}_vs_${ref_name}_filtered.paf


## loop over the different assemblies
# fastq_dir=/user/antwerpen/205/vsc20587/aitg_scratch/jvdabbeele/Tbg_de_novo/
# for assembly in ${fastq_dir}/*.fasta; do sbatch --export=assembly_genome=${assembly} /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_map_to_ref.slurm; done


