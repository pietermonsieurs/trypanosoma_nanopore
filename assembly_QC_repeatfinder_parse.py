#!/usr/bin/env python3

import argparse
import os

my_debug = 0

def parse_repeatfinder_output(file_path):
    repeats = []
    sequence_id = ""
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split()

            if line.startswith('Sequence: '):
                sequence_id = line.strip().replace("Sequence: ", "")
                continue

            if len(parts) < 9:
                my_debug and print(f"Warning: Skipping malformed line: {line.strip()}")
                continue  # Skip malformed lines

            repeat_info = {
                'sequence_id': sequence_id,
                'start': parts[0],
                'end': parts[1],
                'repeat_length': int(parts[2]),
                'copynumber': float(parts[3]),
                'consensus_size': int(parts[4]),
                'total_length': int(parts[7]),
                'repeat': parts[13],
                'repeat_full': parts[14]
            }
            print(f"Parsed repeat: {repeat_info}")  # Debug print
            print(f"line was {line.strip()}")  # Debug print
            repeats.append(repeat_info)
    return repeats

def main():
    parser = argparse.ArgumentParser(description='Parse RepeatFinder output and extract repeat information.')
    parser.add_argument('input_file', help='Path to the RepeatFinder output file')
    parser.add_argument('output_file', help='Path to the output file to save parsed repeat information')
    parser.add_argument('--minsize', type=int, default=0, help='Minimum size of the repeat regions (default: 0)')
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.input_file):
        print(f"Error: The file {args.input_file} does not exist.")
        return
    
    repeats = parse_repeatfinder_output(args.input_file)
    
    with open(args.output_file, 'w') as f:
        f.write("sequence_id\tstart\tend\trepeat_length\tcopy_number\tconsensus_size\ttotal_length\trepeat\trepeat_full\n")
        for repeat in repeats:
            print(f"Processing repeat: {repeat}")  # Debug print
            print(f"Repeat total length: {repeat['total_length']}")  # Debug print
            if int(repeat['total_length']) >= args.minsize:
                # f.write(f"{repeat['sequence_id']}\t{repeat['start']}\t{repeat['end']}\t{repeat['period']}\t{repeat['copynumber']}\t{repeat['consensus_size']}\t{repeat['total_length']}\t{repeat['repeat']}\t{repeat['repeat_full']}\n")
                f.write(f"{repeat['sequence_id']}\t{repeat['start']}\t{repeat['end']}\t{repeat['repeat_length']}\t{repeat['copynumber']}\t{repeat['consensus_size']}\t{repeat['total_length']}\t{repeat['repeat']}\n")
    
    print(f"Parsed repeat information has been saved to {args.output_file}")


if __name__ == "__main__":
    main()