## run on the GPU computer, where verkko is installed
src_dir=/mnt/data/Pore-C/

## use the MM05 reads obtained with long filt. Copying from the 
## supercomputer required
cd ${src_dir}/data/
scp calcua:/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/mm05_combined_filtlong.fastq.gz ./

## do correction of the fastq file with dorado correct, so that you can use them as input
## for verkko with the --hifi option. First unzip the fastq file as it cannot work with 
## regularly compressed fastq files (only bgzipped fastq files are supported). 
gunzip -c mm05_combined_filtlong.fastq.gz > mm05_combined_filtlong.fastq

## also export the LD_LIBRARY_PATH to include the dorado library path, otherwise it will not work
export LD_LIBRARY_PATH=/opt/ont/dorado/lib:$LD_LIBRARY_PATH

## run dorada correct
dorado correct \
    --device cuda:0 \
    --threads 22 \
    mm05_combined_filtlong.fastq \
    > mm05_combined_filtlong.corrected.fasta

## concateate poreC data
porec_dir=/mnt/data/Pore_C/Pore-C/no_sample_id/20260803_1331_MN49550_FBH13206_c4bd104f/fastq_pass/
cat ${porec_dir}/*.fastq.gz > ${src_dir}/data/poreC_combined.fastq.gz

## run verkko
results_dir=/mnt/data/Pore_C/assembly_mm05/results/
data_dir=/mnt/data/Pore_C/assembly_mm05/data/

conda activate verkko
verkko -d asm \
    --hifi ${data_dir}mm05_combined_filtlong.corrected.fasta \
    --nano ${data_dir}mm05_combined_filtlong.fastq \
    --porec ${data_dir}poreC_combined.fastq.gz