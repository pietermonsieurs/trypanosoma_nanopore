#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu

## load modules
module load minimap2
module load SAMtools

## test files
# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TcDm28cT2T.contigs.fa
# fastq_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_nobc/barcode01.fastq.gz
threads=28

minimap_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/minimap/
file_prefix=${fastq_file%.fastq.gz}
file_prefix=${file_prefix##*/}_${ref_genome_short}
bam_file_prefix=${minimap_dir}/${file_prefix}
echo $bam_file_prefix

## run minimap2 and feed into samtools to obtain .bam file 
minimap2 \
    -ax map-ont \
    $ref_genome \
    $fastq_file | \
    samtools sort -@$threads -o ${bam_file_prefix}.bam

samtools index ${bam_file_prefix}.bam
samtools flagstat ${bam_file_prefix}.bam > ${bam_file_prefix}.flagstat
samtools stats ${bam_file_prefix}.bam > ${bam_file_prefix}.stats
samtools depth -a ${bam_file_prefix}.bam > ${bam_file_prefix}.depth

## do a quality filtering of the .bam file at mapping quality 20
samtools view -b -q 20 ${bam_file_prefix}.bam > ${bam_file_prefix}.mq20.bam
samtools index ${bam_file_prefix}.mq20.bam
samtools flagstat ${bam_file_prefix}.mq20.bam > ${bam_file_prefix}.mq20.flagstat
samtools idxstats ${bam_file_prefix}.mq20.bam > ${bam_file_prefix}.mq20.idxstats
samtools depth -a ${bam_file_prefix}.mq20.bam > ${bam_file_prefix}.mq20.depth




## run for all T. cruzi, T.brucei and T. congolense genomes
# cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_nobc/
# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TcDm28cT2T.contigs.fa
# ref_genome_short=TcDm28cT2T
# for barcode in {01..02}; do sbatch --export=ref_genome=${ref_genome},ref_genome_short=${ref_genome_short},fastq_file=${PWD}/barcode${barcode}.fastq.gz /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/mapping_minimap.slurm; done

# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-67_TbruceiTREU927_Genome.fasta
# ref_genome_short=TREU927
# for barcode in {03..04}; do sbatch --export=ref_genome=${ref_genome},ref_genome_short=${ref_genome_short},fastq_file=${PWD}/barcode${barcode}.fastq.gz /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/mapping_minimap.slurm; done

# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-68_TbruceigambienseDAL972_Genome.fasta
# ref_genome_short=DAL972
# for barcode in {03..04}; do sbatch --export=ref_genome=${ref_genome},ref_genome_short=${ref_genome_short},fastq_file=${PWD}/barcode${barcode}.fastq.gz /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/mapping_minimap.slurm; done

# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-68_TcongolenseIL3000_Genome.fasta
# ref_genome_short=IL3000
# for barcode in {05..06}; do sbatch --export=ref_genome=${ref_genome},ref_genome_short=${ref_genome_short},fastq_file=${PWD}/barcode${barcode}.fastq.gz /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/mapping_minimap.slurm; done

# ref_genome=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/refgenomes/TriTrypDB-68_TcongolenseIL3000_2019_Genome.fasta
# ref_genome_short=IL3000_2019
# for barcode in {05..06}; do sbatch --export=ref_genome=${ref_genome},ref_genome_short=${ref_genome_short},fastq_file=${PWD}/barcode${barcode}.fastq.gz /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/mapping_minimap.slurm; done