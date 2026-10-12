#!/bin/bash -l

#SBATCH --partition=broadwell          # choose partition to run the nextflow head process on
#SBATCH --job-name=nextflow            # create a short name for your job
#SBATCH --nodes=1                      # node count
#SBATCH --cpus-per-task=1              # only 1 cpu cores is needed to run the nextflow head process
#SBATCH --mem-per-cpu=4G               # memory per cpu (4G is default for most partitions)
#SBATCH --time=00:05:00                # total run time limit (HH:MM:SS)
#SBATCH --account=ap_itg_tryps    # set project account


# Load the available Nextflow module.
module load Nextflow/24.10.2

## create one file with all fastq reads, concatentating different fastq files from different directories
## should only be done once. However, wf-pore-c can also take a directory with multiple fastq files as 
## input, so this step is not strictly necessary.
porec_dir=/data/antwerpen/grp/ap_itg_tryps/Nanopore_PoreC_20260811/no_sample_id/20260803_1331_MN49550_FBH13206_c4bd104f/fastq_pass/
# src_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/
# porec_file=${src_dir}/data/fastq_porec/poreC_combined.fastq.gz
# cat ${porec_dir}/*.fastq.gz > ${porec_file}

## other input and output parameters
assembly_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/results/assembly_mm05/mm05_flye_assembly/purge_dups/purged.fa
config_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/data/fastq_porec/epi2me.config

# Set Apptainer/Singularity environment variables to define caching and tmp
# directories. These are used during the conversion of Docker images to
# Apptainer/Singularity ones.
# These lines can be omitted if the variables are already set in your `~/.bashrc` file.
export APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/cache"
export APPTAINER_TMPDIR="${VSC_SCRATCH}/apptainer/tmp"
export NXF_APPTAINER_CACHEDIR="${VSC_SCRATCH}/apptainer/nextflow_cache"

# Launch Nextflow head process.
# Provide the vsc_calcua profile to use this config and let Nextflow schedule tasks
# using the Slurm job scheduler. For local execution on a single node, see below.
# Note that multiple profiles can be stacked, and here we use the built-in test profile
# of the nf-core/rnaseq pipeline for demonstration purposes.

# wget https://raw.githubusercontent.com/nf-core/configs/master/conf/vsc_calcua.config

# nextflow run epi2me-labs/wf-pore-c \
#   -profile singularity,vsc_calcua,single_node \
#   -fastq 

nextflow run epi2me-labs/wf-pore-c \
    --fastq ${porec_dir} \
    --ref ${assembly_file} \
    --cutter NlaIII \
    --pairs --mcool --hi_c --paired_end --bed \
    --threads 16 \
    -c ${config_file} \
    -profile singularity \
    --out_dir wf-pore-c-output

# -c vsc_calcua.config \
