## check the MM05 fastq file from the first run
fastq_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_mm05/
mkdir -p $fastq_dir

## copy from batch 1 (september 2025)
cp /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/barcode04.fastq.gz $fastq_dir/mm05_batch_september2025_barcode04.fastq.gz

## copy from batch 2 (May 2026 - opendoor day)
cat /user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_TbgMM05_20260530/no_sample_id/20260530_1700_MN49550_FBG64658_4febd707/fastq_pass/barcode20/*fastq.gz > $fastq_dir/mm05_batch_may2026_barcode20.fastq.gz

## copy from batch 3 (May 2026 - opendoor day library on expired flow cell)


## concatenate all the fastq files into one
cat $fastq_dir/mm05_batch_september2025_barcode04.fastq.gz \
    $fastq_dir/mm05_batch_may2026_barcode20.fastq.gz \
    > $fastq_dir/mm05_combined.fastq.gz

    