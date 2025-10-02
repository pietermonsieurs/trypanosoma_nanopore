## installation of the nanoQC software using pip
# module load Python/3.13.1-GCCcore-14.2.0
# export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.13/site-packages/:${PYTHONPATH}"
# pip3 install -I  --prefix=/user/antwerpen/205/vsc20587/data/software/python_lib/ nanoQC


## specify the deeptools binary and export the PythonPath so it can be run
## with a local binfile
module load Python/3.13.1-GCCcore-14.2.0
export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.13/site-packages/:${PYTHONPATH}"
export PATH=/user/antwerpen/205/vsc20587/data/software/python_lib/bin/:${PATH}


## run QC with nanoQC
cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq/
for barcode in {01..06}; do
    fastq_file=barcode${barcode}.fastq.gz
    echo ${fastq_file}
    nanoQC ${fastq_file} -o /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/QC/barcode${barcode}_nanoQC.html
done


