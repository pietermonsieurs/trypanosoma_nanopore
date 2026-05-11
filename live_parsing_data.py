#!/usr/bin/env python3

## script to parse live data from a directory, where there is a growing file, but every hour
## a new file is created. The script will keep parsing the latest (growing) file, and keep 
## track of those nanopore sequencing reads it already processed, and only report the stats
## for the new reads. The script will also keep track of the latest file, and if a new file 
## is created, it will switch to that file. Additional complication: the script is monitoring
## one data directory (fastq_pass) which contains different directories, each for a different
## sample. The script should be able to monitor all samples, and report the stats for each 
## sample separately.

import os
from Bio import SeqIO
