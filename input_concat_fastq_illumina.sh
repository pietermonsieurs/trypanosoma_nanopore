## concatenate fastq files from Illumina sequencing runs and write them 
## to output directory on the scratch

output_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_illumina/

cd /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_illumina_20251105/X204SC25101290-Z01-F002/01.RawData/

for i in $(ls -d */); do
    sample=$(echo $i | cut -d'/' -f1)
    echo "Processing sample: $sample"
    cat $i/*_1.fq.gz > ${output_dir}/${sample}_R1.fastq.gz
    cat $i/*_2.fq.gz > ${output_dir}/${sample}_R2.fastq.gz
done