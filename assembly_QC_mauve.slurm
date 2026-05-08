#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_tryps


## installation of Mauve software
# https://darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/mauve_linux_snapshot_2015-02-13.tar.gz
# 


## parameter settings
contig_dir=/user/antwerpen/205/vsc20587/aitg_scratch/jvdabbeele/Tbg_de_novo/
ref_dir=


