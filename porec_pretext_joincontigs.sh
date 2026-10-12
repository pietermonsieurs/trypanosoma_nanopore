## script to join two contigs in a fasta assembly file, adding a gap of N's between them
## those contigs were identified as mapping next to ehach other based on the contact map
## based on PoreC data, but the assembly was not able to join them
## script to concat two contigs, add a gap of N's, and rebuild the assembly fasta file after
## removing the two original contigs from the assembly

module load SAMtools

cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_final_niki

ASM="MM05_medaka_noragtag.fasta"
OUT="MM05_medaka_noragtag.joined.fasta"
GAP_N=100
NEW_NAME="scaffold_11_5_np1212_joined"

samtools faidx "$ASM"

TMP=$(mktemp -d)

# 1. Build the joined record: scaffold_11 -> N-gap -> scaffold_5, both as-is (no revcomp needed)
{
  echo ">${NEW_NAME}"
  samtools faidx "$ASM" scaffold_11_np1212 | tail -n +2 | tr -d '\n'
  printf 'N%.0s' $(seq 1 "$GAP_N")
  samtools faidx "$ASM" scaffold_5_np1212  | tail -n +2 | tr -d '\n'
  echo
} | fold -w 60 > "$TMP/joined.fasta"

# 2. Rebuild the assembly: every sequence except the two originals, then append the join
samtools faidx "$ASM" $(cut -f1 "${ASM}.fai" | grep -vE '^(scaffold_11_np1212|scaffold_5_np1212)$') > "$OUT"
cat "$TMP/joined.fasta" >> "$OUT"

samtools faidx "$OUT"
rm -r "$TMP"

echo "sequences before: $(grep -c '^>' "$ASM")   after: $(grep -c '^>' "$OUT")"