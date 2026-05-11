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

def parse_live_data(data_dir):
    # Dictionary to keep track of the latest file for each sample
    latest_files = {}
    # Dictionary to keep track of the reads already processed for each sample
    processed_reads = {}

    while True:
        # List all directories in the data directory
        for sample_dir in os.listdir(data_dir):
            sample_path = os.path.join(data_dir, sample_dir)
            if os.path.isdir(sample_path):
                # List all files in the sample directory
                files = [f for f in os.listdir(sample_path) if f.endswith('.fastq')]
                if not files:
                    continue
                # Get the latest file based on modification time
                latest_file = max(files, key=lambda f: os.path.getmtime(os.path.join(sample_path, f)))
                latest_file_path = os.path.join(sample_path, latest_file)

                # Check if we have a new file for this sample
                if sample_dir not in latest_files or latest_files[sample_dir] != latest_file:
                    print(f"New file detected for sample {sample_dir}: {latest_file}")
                    latest_files[sample_dir] = latest_file
                    processed_reads[sample_dir] = set()  # Reset processed reads for new file

                # Parse the latest file and report stats for new reads
                with open(latest_file_path, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fastq'):
                        read_id = record.id
                        if read_id not in processed_reads[sample_dir]:
                            processed_reads[sample_dir].add(read_id)
                            # Report stats for the new read (e.g., length, quality)
                            print(f"Sample: {sample_dir}, Read ID: {read_id}, Length: {len(record.seq)}, Quality: {record.letter_annotations['phred_quality']}")

        # Sleep for a short period before checking again (e.g., 10 seconds)
        print("Waiting for new data...")
        time.sleep(10)


if __name__ == "__main__":
    data_directory = "/Users/pmonsieurs/programming/trypanosoma_nanopore/data/nanopore_run/fastq_pass/"  # Update this path to your data directory
    parse_live_data(data_directory)
