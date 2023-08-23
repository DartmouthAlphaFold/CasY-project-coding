import argparse
import os

def find_matching_brackets(sequence):
    stack = []
    positions = []

    for i, char in enumerate(sequence):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                positions.append((stack.pop(), i))
    return positions


def make_pairing_file(dot_file):
    with open(dot_file, 'r') as f:
        chars = f.read()

    positions = find_matching_brackets(chars)
    pairings = ""
    for start, stop in positions:
        pairings += f"{start+1}\t{stop+1}\t1\n"
        
    print(f"Pairings:\n{pairings}")
    pairing_file = f'{os.path.split(dot_file)[0]}/pairing.txt'
    with open(pairing_file, 'w') as f:
        f.write(pairings)
    return pairing_file


def build_model(sequence_file, pairing_file, output_file):
    from Bio import SeqIO
    seq_record = SeqIO.parse(sequence_file, 'fasta')[0]
    length = len(seq_record.seq)
    name = seq_record.id
    print(f'Sequence: {name}')
    
    from chimerax.core.commands import run
    run(session, 'close session; log clear')
    run(session, f'rna path {pairing_file} length {length}')
    run(session, f'rna model {sequence_file} #1 start 1')
    run(session, f'save {output_file}')
    
    
def main():
    parser = argparse.ArgumentParser(description="Simple CLI tool to build RNA structure using ChimeraX "
                                     "using secondary structure predicted by Snapgene")
    parser.add_argument('--fasta', '-f', required=True, help='RNA sequence in FASTA format')
    parser.add_argument('--pairing', '-p', required=True, help='Pairings in dot-bracket (txt) format from Snapgene')
    parser.add_argument('--output', '-o', default='./rna.pdb', help = "Output file in PDB format")
    args = parser.parse_args()
    
    pairing_file = make_pairing_file(args.pairing)
    build_model(args.fasta, pairing_file, args.output)
    
    
if __name__ == '__main__':
    main()