## specify the directories containing the fastq files for both runs. For the first directory, the fastq
## files are directly located in that directory, for the second run there is one subdirectory per sample, 
## and the fastq files are located in those subdirectories. 
fastq_dir_run1=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20250925/no_sample_id/20250925_1219_MN49550_FBD31148_f66f99ac/fastq_pass_nobarcode_update/
fastq_dir_run2=/user/antwerpen/205/vsc20587/aitg_data/jvdabbeele/Nanopore_20251130/20251127_1547_MN49550_FBD95189_6984b7a7/fastq_dorado_updated/

## specify the output directory
out_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_merging/

## concatenate all fastq files of Tbg
cat \
    ${fastq_dir_run1}/barcode03.fastq.gz \
    ${fastq_dir_run1}/barcode04.fastq.gz \
    ${fastq_dir_run2}/barcode09/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq \
    ${fastq_dir_run2}/barcdoe10/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq \
    ${fastq_dir_run2}/barcode11/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq \
    > ${out_dir}/Tbg_merged.fastq.gz

## sanity check wehther it contains the sum of all fastq files 


## concatenate all fastq files of Tcongo
cat \
    $fastq_dir_run1/barcode05.fastq.gz \
    $fastq_dir_run1/barcode06.fastq.gz \
    ${fastq_dir_run2}/barcode13/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq \
    > ${out_dir}/Tcongo_merged.fastq.gz


## concatenate all fastq files of Tcongo MSOROM7
cat \
    $fastq_dir_run1/barcode05.fastq.gz \
    ${fastq_dir_run2}/barcode13/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq \
    > ${out_dir}/Tcongo_MSOROM7_merged.fastq.gz


## do a quick quality control of the sequencing data, both from the individual fastq files and from the merged fastq files. 
cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_merging/
for file in *.fastq.gz; do
    echo "Processing $file"
    lenght_file=${file/.fastq.gz/.lengths.txt}
    zcat ${file} | awk 'NR%4==2{print length($0)}' > ${lenght_file}
done

zcat $fastq_dir_run1/barcode03.fastq.gz | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode03_lengths.txt
zcat $fastq_dir_run1/barcode04.fastq.gz | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode04_lengths.txt
zcat $fastq_dir_run1/barcode05.fastq.gz | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode05_lengths.txt
zcat $fastq_dir_run1/barcode06.fastq.gz | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode06_lengths.txt
more ${fastq_dir_run2}/barcode09/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode09_lengths.txt
more ${fastq_dir_run2}/barcode10/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode10_lengths.txt
more ${fastq_dir_run2}/barcode11/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode11_lengths.txt
more ${fastq_dir_run2}/barcode13/FBD95189_fastq_pass_6984b7a7_7a260689_0.fastq | awk 'NR%4==2{print length($0)}' > ${out_dir}/barcode13_lengths.txt