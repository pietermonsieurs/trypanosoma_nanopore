## loop over a csv file with as columsn the sample name, the nanopore fastq file, 
## and both Illumina fastq files (forward and reverse), and run the script for each sample
## assembly_QC_manual_assembly_quast.slurm

while IFS=, read -r sample_name nanopore_fastq illumina_forward illumina_reverse
do
    echo "Processing sample: $sample_name"
    echo "Nanopore fastq: $nanopore_fastq"
    echo "Illumina forward: $illumina_forward"
    echo "Illumina reverse: $illumina_reverse"    
    sbatch --export=sample_name="$sample_name",nanopore_reads="$nanopore_fastq",illumina_R1_reads="$illumina_forward",illumina_R2_reads="$illumina_reverse" /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/bin/assembly_QC_manual_assembly_quast.slurm
    echo -e "submitted: $sample_name\n\n"
done < /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_data_illumina_nanopore_Tbg.csv

