#!/usr/bin/env python3

import os
import re
import argparse
from Bio import SeqIO

def get_sequence(input_files):
    '''
    Fetch all sequences from a list of FASTA files
    '''
    seq_map = {}
    for input_file in input_files:
        for seq_record in SeqIO.parse(input_file.strip(), "fasta"):
            if len(seq_record.seq) == 0:
                continue
            sequence = re.sub('[\*\-]', '', str(seq_record.seq))  # remove '*' and '-' characters
            seq_map[str(seq_record.description)] = sequence.upper()
    print(f"\nFound {len(seq_map)} sequences: \n{''.join(seq_map.keys())}")
    return seq_map


def write_sequence(seq_map, output):
    with open(output, 'w') as f:
        for name, sequence in seq_map.items():
            f.write(f">{name}\n{sequence}\n")
    

def main():
    parser = argparse.ArgumentParser(description="Command line tool to find sequences from FASTA files/folders and export to a single file")
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help = 'Files or folders that contain sequences in fasta format. Each sequence must have unique name.')
    parser.add_argument('--output', '-o', required=True,
                        help = "Output file")
    args = parser.parse_args()
    
    input = []
    for path in args.input:
        if os.path.exists(path):
            if os.path.isfile(path):
                input.append(path)
            elif os.path.isdir(path):
                seq_files = [os.path.join(path,file) for file in os.listdir(path) if file.endswith(('.fa', '.fasta', '.faa', '.fna', '.ffn', '.frn', '.fas'))]
                input += seq_files
            else:
                print(f"{path} is neither a folder nor a file")
        else:
            print(f"Cannot find {path}") 
    print(f"Working on {len(input)} files:\n {', '.join(input)}")
    write_sequence (get_sequence(input), args.output)


if __name__ == "__main__":
    main()