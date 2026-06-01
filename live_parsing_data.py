#!/usr/bin/env python3

## script to parse live data from a directory, where there is a growing file, but every hour
## a new file is created. The script will keep parsing the latest (growing) file, and keep 
## track of those nanopore sequencing reads it already processed, and only report the stats
## for the new reads. The script will also keep track of the latest file, and if a new file 
## is created, it will switch to that file. Additional complication: the script is monitoring
## one data directory (fastq_pass) which contains different directories, each for a different
## sample. The script should be able to monitor all samples, and report the stats for each 
## sample separately.

import gzip
import os
from Bio import SeqIO
import time

def parse_live_data(data_dir):
    ## Dictionary to keep track of the list of files already seen for each sample, to detect new files
    processed_files = {}
    ## Dictionary to keep track of the reads already processed for each sample
    processed_reads = {}
    ## Dictionary to keep track of the total number of base pairs processed for each sample
    total_base_pairs = {}
    ## Dictionary to keep track of the lengths of the reads for each sample
    read_lengths = {}


    while True:
        # List all directories in the data directory
        for sample_dir in os.listdir(data_dir):
            if not sample_dir in ('barcode09', 'barcode10', 'barcode11', 'barcode12', 'barcode13'):
                continue

            print(f"Checking sample directory: {sample_dir}")
            time.sleep(1)  # Sleep briefly to avoid overwhelming the file system
            sample_path = os.path.join(data_dir, sample_dir)
            if os.path.isdir(sample_path):
                ##  List all files in the sample directory. Those are gzipped fastq files, 
                ## but we can only parse the latest one, which is growing. The others are 
                ## finished and should not be parsed anymore.
                files = [f for f in os.listdir(sample_path) if f.endswith('.fastq.gz')]
                if not files:
                    continue
                    
                for fastq_file in files:
                    print(f"Found file: {fastq_file} in sample directory: {sample_dir}")

                    if fastq_file not in processed_files.get(sample_dir, []):
                        print(f"New file detected for sample {sample_dir}: {fastq_file}")
                        processed_files[sample_dir] = processed_files.get(sample_dir, []) + [fastq_file]
                        processed_reads[sample_dir] = set()  # Reset processed reads for new file
                        time.sleep(2)

                    ## create the full path to the fastq_file
                    fastq_file = os.path.join(sample_path, fastq_file)

                    ## Parse the new file and report stats for new reads. Also keep the total number 
                    ## of reads processed for each sample, and report that as well. Additionally, make
                    ## a list containing all the read lengths to make a plot of the read length distribution 
                    ## at the end of the run. Also report the total number of bases processed for each sample, 
                    ## and the average read length after running all the records
                    with gzip.open(fastq_file, 'rt') as handle:
                        for record in SeqIO.parse(handle, 'fastq'):
                            read_id = record.id
                            if read_id not in processed_reads[sample_dir]:
                                processed_reads[sample_dir].add(read_id)
                                ## Report stats for the new read (e.g., length, quality)
                                print(f"Sample: {sample_dir}, Read ID: {read_id}, Length: {len(record.seq)}, sequence: {record.seq}")
                                # time.sleep(0.1)  # Sleep briefly to avoid overwhelming the output

                                ## Update total base pairs
                                total_base_pairs[sample_dir] = total_base_pairs.get(sample_dir, 0) + len(record.seq)    

                                ## add the read length to the list of read lengths for this sample
                                read_lengths[sample_dir] = read_lengths.get(sample_dir, []) + [len(record.seq)]

                    ## print some stats for the sample, such as total base pairs processed and average read length
                    if sample_dir in total_base_pairs and sample_dir in read_lengths:
                        average_read_length = total_base_pairs[sample_dir] / len(processed_reads[sample_dir]) if processed_reads[sample_dir] else 0
                        print(f"Sample: {sample_dir}, Total Base Pairs: {total_base_pairs[sample_dir]}, Average Read Length: {average_read_length}")
                        print(read_lengths[sample_dir])  # Print the list of read lengths for this sample
                        time.sleep(2)

        # Sleep for a short period before checking again (e.g., 10 seconds)
        print("Waiting for new data...")
        time.sleep(10)


if __name__ == "__main__":
    data_directory = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run/fastq_pass/"  # Update this path to your data directory
    parse_live_data(data_directory)
