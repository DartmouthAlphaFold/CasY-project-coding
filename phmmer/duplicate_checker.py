#!/usr/bin/env python3
import hashlib
import re
from Bio import SeqIO
import os
import argparse


def get_sequence(sequence_files):
    seq_map = {}
    for sequence_file in sequence_files:
        for seq_record in SeqIO.parse(sequence_file.strip(), "fasta"):
            if len(seq_record.seq) == 0:
                continue
            sequence = re.sub('[\*\-]', '', str(seq_record.seq).upper())  # remove '*' and '-' characters
            seq_map[str(seq_record.id)] = sequence
    print(f'Found {len(seq_map)} sequences')
    return seq_map


def compare_strings_pairwise(seq_map):
    hash_dict = {}  # Dictionary to store hash values
    
    # Generate hash values for each string
    for name, sequence in seq_map.items():
        sha256_hash = hashlib.sha256()
        sha256_hash.update(sequence.encode('utf-8'))
        hash_value = sha256_hash.hexdigest()
        
        if hash_value in hash_dict:
            hash_dict[hash_value].append(name)
        else:
            hash_dict[hash_value] = [name]

    # Compare strings with the same hash value
    for hash_value, sequence_list in hash_dict.items():
        if len(sequence_list) > 1:
            print(f"Sequences {', '.join(sequence_list)} are identical.")
                

def main():
    parser = argparse.ArgumentParser(description="Command line tool to read fasta files and search for identical sequences using hashing")
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help = 'Files or folders that contain sequences in fasta format. Each sequence must have unique name.')
    args = parser.parse_args()
    
    input = []
    for path in args.input:
        if os.path.exists(path):
            if os.path.isfile(path):
                input.append(path)
            elif os.path.isdir(path):
                seq_files = [os.path.join(path,file) for file in os.listdir(path) if file.endswith(('.fa', '.fasta', '.faa', '.fna', 'ffn', 'frn'))]
                input += seq_files
            else:
                print(f"{path} is neither a folder nor a file")
        else:
            print(f"Cannot find {path}")
    print(f"Working on {', '.join(input)}")    
    compare_strings_pairwise(get_sequence(tuple(input)))


if __name__ == "__main__":
    main()