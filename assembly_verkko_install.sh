## 1. installation on the supercomputer 

## create space for verko environment
mkdir /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/software
cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/software

## create a yml file with the following content and save it as verkko_env.yml
vi verkko.yml

# name: verkko
# channels:
#   - conda-forge
#   - bioconda
#   - defaults
# dependencies:
#   - verkko


## use HPC friendly way to create a conda environment from the yml file
module load hpc-container-wrapper
conda-containerize new --prefix /user/antwerpen/205/vsc20587/scratch/trypanosoma_nanopore/software/verkko ./verkko.yml



## 2. installation on the GPU computer

## installation of conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh

conda create -n verkko -c conda-forge -c bioconda -c defaults verkko
conda activate verkko
