#!/bin/bash -l

#SBATCH --job-name=map_porec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=05:00:00
#SBATCH --account=ap_itg_tryps

module purge
module load calcua/all
module load minimap2
module load SAMtools
module load Pysam

## export the path to the pretextmap binary
export PATH=/user/antwerpen/205/vsc20587/data/software/PretextMap/builddir/:${PATH}

## define variables
# scr_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki/
# scr_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/control_DAL972/
scr_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/control_TREU927/


## select the correct PREFIX depending on the assembly you want to run it for
# PREFIX=MM05_medaka_noragtag
# PREFIX=MM05_medaka_ragtag
# PREFIX=MM05_nomedaka_noragtag
# PREFIX=MM05_nomedaka_ragtag
# PREFIX=MM05_medaka_noragtag.joined

# PREFIX=TriTrypDB-68_TbruceigambienseDAL972_Genome
PREFIX=TriTrypDB-67_TbruceiTREU927_Genome_nuclear

ASM=${scr_dir}/${PREFIX}.fasta
PORE_C=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_porec/poreC_combined.fastq.gz
OUT=${scr_dir}/
THREADS=28

## create log directory
mkdir -p ${OUT}/logs

## do indexing of the assembly if not already done
[ -f "${ASM}.fai" ] || samtools faidx "$ASM"

## run minimap2 to map the poreC reads to the assembly, filter for primary 
## alignments with mapping quality >= 20, and sort by read name
minimap2 -ax map-ont --secondary=no -t "$THREADS" "$ASM" "$PORE_C" \
  | samtools view -F 0x104 -q 20 -b - \
  | samtools sort -n -@ 4 -o "${OUT}/${PREFIX}.namesorted.bam" -

## create contact map with pretext
samtools view -h "${OUT}/${PREFIX}.namesorted.bam" \
  | PretextMap -o "${OUT}/${PREFIX}.pretext" --sortby length


## create contact map, which can later be used as in put in R heatmap
# python3 /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/porec_prextext_contacts_pairs.py \
#   --bam "${OUT}/porec_vs_final.namesorted.bam" \
#   --fai "${ASM}.fai" \
#   --bin-size 20000 \
#   --out-prefix "${OUT}/MM05_v2_porec"

python3 /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/porec_prextext_contacts_pairs.py \
  --bam "${OUT}/${PREFIX}.namesorted.bam" \
  --fai "${ASM}.fai" \
  --bin-size 20000 \
  --min-mapq 20 \
  --max-monomers 250 \
  --out-prefix "${OUT}/${PREFIX}_porec"
