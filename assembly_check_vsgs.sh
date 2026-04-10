#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_mpu


module load BLAST+

## Tb gambiense
vsg_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_VSG_gambiense/results/trinity/TD2_longest/combined_longest_orfs.cdhit.highest_expressed.cds
for barcode in {03..04}; do
    assembly_flye=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly/barcode${barcode}/assembly.fasta
    output_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly/blast_vsg/barcode${barcode}_vsg_blast.tsv

    makeblastdb -in ${assembly_flye} -dbtype nucl -out ${assembly_flye}

    blastn -query ${vsg_file} \
           -db ${assembly_flye} \
           -out $output_file \
           -evalue 1e-50 -outfmt 6
done


## for T. congolense
vsg_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/vsg_fasta/vsgs_tcongo.fasta
for barcode in {05..06}; do
    assembly_flye=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly/barcode${barcode}/assembly.fasta
    output_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly/blast_vsg/barcode${barcode}_vsg_blast.tsv

    makeblastdb -in ${assembly_flye} -dbtype nucl -out ${assembly_flye}

    blastn -query ${vsg_file} \
           -db ${assembly_flye} \
           -out $output_file \
           -evalue 1e-50 -outfmt 6
done




