module load dorado

dorado trim \
  --sequencing-kit SQK-NBD114-24 \
  --emit-fastq \
  /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/trimmed/barcode01.fastq.gz \
  >   /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/trimmed/barcode01.trimmed.fastq



dorado demux \
  -o demuxed_reads/ \
  --kit-name SQK-NBD114-24 \
  --emit-fastq \
  /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/trimmed/barcode01.trimmed.fastq
