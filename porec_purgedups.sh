#!/usr/bin/bash

#SBATCH --job-name=purge_dups
#SBATCH --account ap_itg_tryps
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=28

## installing purge_dups
# cd /user/antwerpen/205/vsc20587/data/software/
# git clone https://github.com/dfguan/purge_dups.git
# cd purge_dups/src/
# make

## loading modules
module load minimap2

## other software settings
purge_dups_bin=/user/antwerpen/205/vsc20587/data/software/purge_dups/
export PATH=${purge_dups_bin}/src/:$PATH

## parameters
# ASM=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/assembly.fasta
# READS=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/mm05_combined.fastq.gz   # the ONT reads used to build it
# OUT_DIR=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/purge_dups
# mkdir -p ${OUT_DIR}

ASM=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki/MM05_medaka_noragtag.fasta
READS=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/mm05_combined.fastq.gz   # the ONT reads used to build it
OUT_DIR=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki/purge_dups_MM05_medaka_noragtag
mkdir -p ${OUT_DIR}


# 1. Align your ONT reads back to the assembly, to measure read depth
minimap2 -x map-ont -t 16 "$ASM" "$READS" | gzip -c > ${OUT_DIR}/reads_vs_asm.paf.gz

# 2. Turn that into a coverage histogram and pick low/mid/high cutoffs
pbcstat ${OUT_DIR}/reads_vs_asm.paf.gz -O ${OUT_DIR}         # -> produces PB.stat, PB.base.cov
calcuts ${OUT_DIR}/PB.stat > ${OUT_DIR}/cutoffs

# (optional but recommended, esp. for ONT) sanity-check the histogram/cutoffs visually
python ${purge_dups_bin}/scripts/hist_plot.py -c ${OUT_DIR}/cutoffs ${OUT_DIR}/PB.stat ${OUT_DIR}/coverage_hist.png

# 3. Split the assembly at gaps, then self-align it to find duplicated contigs
split_fa "$ASM" > "$ASM.split"
minimap2 -x asm5 -DP -t 16 "$ASM.split" "$ASM.split" | gzip -c > "$ASM.split.self.paf.gz"

# 4. Call the actual duplications
purge_dups -2 -T cutoffs -c ${OUT_DIR}/PB.base.cov "$ASM.split.self.paf.gz" > ${OUT_DIR}/dups.bed

# 5. Produce the final purged assembly + the removed haplotigs
cd ${OUT_DIR}
get_seqs -e ${OUT_DIR}/dups.bed "$ASM"