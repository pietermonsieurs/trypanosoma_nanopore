out_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/

cd /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass/

## concat fastq files for the nepal samples (last 2 sample 11 and 12 are not included, 
## those are samples from Grenada and not Nepal)
for barcode in {01..06}
do
    echo ${barcode}
    cat barcode${barcode}/*.fastq.gz > ${out_dir}/barcode${barcode}.fastq.gz
done
