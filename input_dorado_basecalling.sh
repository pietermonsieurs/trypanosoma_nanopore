cd /mnt/data/Community_test/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/pod5_pass/

# out_dir=/mnt/data/Community_test/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode/
out_dir=/mnt/data/Community_test/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/

for barcode in {01..03}; do
    barcode_name=barcode${barcode}
    dorado basecaller \
        -o ${out_dir}/${barcode_name} \
        --emit-fastq \
        --kit-name SQK-NBD114-24 \
        --trim all \
        sup ${barcode_name}
done


## run for the second experiment. Slightly different approach as their is no
## splitting up per barcode yet, all pod5 files in one directory so still need
## to do the demultiplexing
pod5_dir=/mnt/data/pool27112025/Vol_au_vent/pool27112025/20251127_1547_MN49550_FBD95189_6984b7a7/pod5/
output_dir=/mnt/data/pool27112025/Vol_au_vent/pool27112025/20251127_1547_MN49550_FBD95189_6984b7a7/fastq_dorado_updated/

dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
    $pod5_dir \
    --recursive \
    --device cuda:0 \
    --emit-fastq \
    --trim all \
    --kit-name SQK-NBD114-24 \
    -o ${output_dir}

