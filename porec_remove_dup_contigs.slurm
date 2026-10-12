#!/bin/bash -l

#SBATCH --job-name=map_porec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=05:00:00
#SBATCH --account=ap_itg_tryps

module purge
module load calcua/2024a
module load minimap2
module load BEDTools

export PATH=/user/antwerpen/205/vsc20587/data/software/:${PATH}
threads=28

BASE=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki/
ASM=$BASE/MM05_medaka_noragtag.joined.fasta
OUT=$BASE/MM05_medaka_noragtag.joined.nodups
mkdir -p $OUT

# split by length (adjust 1 Mb cutoff to where your megabase chromosomes start)
seqkit seq -m 500000 $ASM > $OUT/large.fa
seqkit seq -M  499999 $ASM > $OUT/small.fa

minimap2 -x asm20 -t $threads --secondary=no $OUT/large.fa $OUT/small.fa > $OUT/small_vs_large.paf


## get some stats
awk -v OFS='\t' '$10/$11 >= 0.90 {print $1,$3,$4}' $OUT/small_vs_large.paf \
 | sort -k1,1 -k2,2n | bedtools merge \
 | awk -v OFS='\t' '{c[$1]+=$3-$2} END{for(k in c) print k,c[k]}' \
 | sort > $OUT/covered.tsv

seqkit fx2tab -nl $OUT/small.fa | sort \
 | join -t $'\t' -a1 -e0 -o 1.1,1.2,2.2 - $OUT/covered.tsv \
 | awk -v OFS='\t' '{print $1,$2,$3,$3/$2}' > $OUT/contained.tsv

 ## based on the output of contained.tsv, remove the contigs that overlap with the 
 ## chromosomes

OUT=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki/MM05_medaka_noragtag.joined.nodups

cat > $OUT/remove.txt <<EOF
scaffold_33_np1212
scaffold_34_np1212
scaffold_40_np1212
scaffold_52_np1212
EOF

# assembly without the four duplicates
seqkit grep -v -f $OUT/remove.txt $ASM > $OUT/MM05_medaka_noragtag.joined.nodups.fasta

# the removed ones, kept for provenance
seqkit grep    -f $OUT/remove.txt $ASM > $OUT/MM05_medaka_noragtag.joined.duplicates.fasta