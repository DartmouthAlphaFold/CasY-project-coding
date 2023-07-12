#!/usr/bin/env python3
from urllib.request import *
from urllib.parse import urlencode, urlparse
import os
import argparse
import json
from bs4 import BeautifulSoup
from Bio import SeqIO
import time
import re
import sys
import hashlib


def get_sequence(args, input_files):
    '''
    Fetch all sequences from a list of FASTA files
    '''
    seq_map = {}
    for input_file in input_files:
        for seq_record in SeqIO.parse(input_file.strip(), "fasta"):
            if len(seq_record.seq) == 0:
                continue
            sequence = re.sub('[\*\-]', '', str(seq_record.seq))  # remove '*' and '-' characters
            seq_map[str(seq_record.id)] = sequence.upper()
    if args.debug: print(seq_map.items())
    return seq_map


# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers
    
    
def phmmer(args, seq_map):
    '''
    Send API to search protein sequences against databases using pHMMER
    Receive and export subjects found into metadata files
    '''
    
    # get params
    sequence_db=args.database
    range=args.range
    output_format=args.output_format
    output_dir=args.output
    
    opener = build_opener(SmartRedirectHandler())
    install_opener(opener)

    outputs= []
    for name, sequence in seq_map.items():
        parameters = {
            'seq': sequence,
            'seqdb': sequence_db
        }
        enc_params = urlencode(parameters).encode("utf-8")

        # post the search request to the server
        request = Request('https://www.ebi.ac.uk/Tools/hmmer/search/phmmer', enc_params)

        # get the url where the results can be fetched from
        # https://www.ebi.ac.uk/Tools/hmmer/results/{uuid}/score
        try:
            results = urlopen(request)
        except Exception as e:
            print(e)
            sys.exit()
            
        results_url = results.get('location')
        session_id = urlparse(results_url).path.split("/")[-2]

        # modify the range and format in results here
        if output_format == 'html':
            res_params = {'range': f'1,{str(range)}'}
        else:
            res_params = {
                'output': output_format,
                'range': f'1,{str(range)}'
            }

        # add the parameters to your request for the results
        enc_res_params = urlencode(res_params)
        modified_res_url = results_url + '?' + enc_res_params

        # send a GET request to the server
        results_request = Request(modified_res_url)
        try:
            data = urlopen(results_request)
        except Exception as e:
            print (e)
            sys.exit()

        # write the results
        if (data.status == 200):
            output_file = f'{output_dir}/{name}.{session_id}.{output_format}'
            with open(output_file, 'w') as f:
                f.write(data.read().decode('utf-8'))
            outputs.append({'query': name,
                            'length': len(sequence),
                            'filename': output_file})
        else:
            print(f'{data.status}\n{data.msg}\n{data.reason}\n{data.url}')
            
    return outputs


def read_json(args, outputs):
    '''
    Read the metadata in JSON format and export satisfied hits to FASTA format
    '''
    
    for output in outputs:
        query_name = output['query']
        query_length = output['length']
        filename = output['filename']
        
        with open(filename, 'r') as f:
            results = json.load(f)['results']
            num_hits = results['stats']['nhits']
            
            if (num_hits > 0):
                sig_num_hits = 0
                for hit in results['hits']:
                    # skip identical hit (identical hit should only have 1 domain with 100% identity)
                    if (len(hit['domains']) == 1 and hit['domains'][0]['aliId'] == 1.0):
                        continue
                    
                    # skip insignificant hits
                    if (float(hit['evalue']) >= args.evalue):
                        continue
                    
                    # skip hits that have low domain coverage of the query
                    domains_coverage = hit['domains'][-1]['alisqto'] - hit['domains'][0]['alisqfrom']
                    if domains_coverage <= 0.5*query_length:
                        continue
                    
                    sig_num_hits += 1
                    # set URL to retrieve sequence based on database used
                    if args.database in ['uniprotrefprot', 'uniprotkb', 'swissprot', 'alphafold']:
                        hit_seq_url = hit['extlink'] + ".fasta" 
                    # https://www.uniprot.org/uniprot/{id}, then redirected
                    elif args.database == 'pdb':
                        # hit_seq_url = "https://www.rcsb.org/fasta/entry" + hit['name'].split("_")[0]
                        hit_seq_url = hit['extlink'].split("_")[0] + "/fasta"
                        # https://www.ebi.ac.uk/pdbe/entry/pdb/{id}/fasta
                    else:
                        continue   # ensembl, rp, etc
                    
                    # send a GET request to retrieve the full subject sequence
                    hit_seq_request = Request(hit_seq_url)
                    try:
                        hit_seq_response = urlopen(hit_seq_request)
                    except Exception as e:
                        print(f"{filename}: {e.code} {e.msg}")
                    
                    # fetch the sequence, avoid rewritting recorded hits
                    if (hit_seq_response.status == 200):
                        seq_output = f"{args.output}/{hit['name']}.fa"
                        if not os.path.exists(seq_output):
                            with open(seq_output, 'w') as f:
                                content = hit_seq_response.read().decode('utf-8')
                                f.write(content)
                    else:
                        if args.debug: print(f'{hit_seq_response.status}\n{hit_seq_response.msg}\n{hit_seq_response.reason}\n{hit_seq_response.url}\n')
                    
                    if args.debug:
                        print(f"{hit['name']}\n{hit['species']}\n{hit['desc']}\n{hit['evalue']}\n"
                        f"{hit['score']}\n{hit['extlink']}\n") 
                print(f"{query_name}: {sig_num_hits}/{num_hits} hits were chosen")
                
            else:
                print(f"{query_name}: no hits were found")


def check_replicates(args):
    '''
    Check if any duplicate exists among query and subject sequences
    '''
    
    # combine queries and subjects into a map
    subject_files = [os.path.join(args.output, file) for file in os.listdir(args.output) if file.endswith(('.fa', '.fasta', '.faa', '.fna', 'ffn', 'frn'))]
    query_files = args.input.split(",")
    seq_map = get_sequence(args, subject_files + query_files)
    print(f"Matching {len(seq_map)} sequences")
    
    hash_dict = {}  # Dictionary to store hash values
    # generate hash values for each string
    for name, sequence in seq_map.items():
        sha256_hash = hashlib.sha256()
        sha256_hash.update(sequence.encode('utf-8'))
        hash_value = sha256_hash.hexdigest()
        
        if hash_value in hash_dict:
            hash_dict[hash_value].append(name)
        else:
            hash_dict[hash_value] = [name]

    # compare strings with the same hash value
    for hash_value, sequence_list in hash_dict.items():
        if len(sequence_list) > 1:
            print(f"Sequences {', '.join(sequence_list)} are identical.")
            

def str_to_bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        import argparse
        raise argparse.ArgumentTypeError('Boolean value expected.')            
            
            
def main():
    parser = argparse.ArgumentParser(description="Sequence search with profile HMMER on multiple databases "
                                     "using the API available at https://www.ebi.ac.uk/Tools/hmmer/search/phmmer.")
    parser.add_argument('--input', '-i', required=True,
                        help="Files in fasta format containing the query protein sequence(s). Each file can contain multiple sequences.\n"
                        "Files should be separated by commas. All FASTA sequences must have unique names.")
    parser.add_argument('--output', '-o', required=True,
                        help="Output folder")
    parser.add_argument('--evalue', '-e', default = 0.01, type = float,
                        help = "Threshold of e-value. (default: 0.01)")
    parser.add_argument('--range', '-r', default=50, type = int,
                        help="Maximum number of top subjects for each query (default: 50)")
    parser.add_argument('--database', '-db', default='uniprotkb',
                        choices=['uniprotrefprot', 'uniprotkb', 'swissprot', 'pdb', 'alphafold', 
                                 'ensembl', 'ensemblgenomes', 'rp75', 'rp55', 'rp35', 'rp15',
                                 'merops', 'qfo', 'chembl'],
                        help="Reference Proteomes, UniProtKB, SwissProt, PDB, AlphaFold,"
                        "Ensembl (all/human/mouse/zebrafish), Ensembl Genomes (all/bacteria/fungi/metazoa/plants/protists),"
                        "Representative Proteomes (rp75/rp55/rp35/rp15), MEROPS, Quest for Orthologs, ChEMBL (default: uniprotkb)")
    parser.add_argument('--output_format', '-f', default='json',
                        choices=['html', 'json', 'xml', 'yaml', 'text'],
                        help="Format of the metadata file. (default: json)")
    parser.add_argument('--debug', type=str_to_bool, default = False,
                        help="Print detailed info for debugging. (default: False)")
    args = parser.parse_args()
    
    # make output directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    start_time = time.time()
    inputs = args.input.split(",")  # get input files
    sequences = get_sequence(args, inputs)
    outputs = phmmer(args, seq_map=sequences)
    if args.output_format == 'json':
        read_json(args, outputs=outputs)
    check_replicates(args)
    end_time = time.time()
    print(f"Finished in {round(end_time - start_time, 3)} seconds")
    
    
if __name__=="__main__":
    main()
    
    # args = argparse.ArgumentParser().parse_args()
    # args.evalue = 0.01
    # args.debug = False
    # args.database = 'pdb'
    # args.output = "C:/XResearch/Mgnify_CasY1"
    # read_json(args, ['C:/XResearch/Mgnify_CasY1/CasY1.55DD4408-1DB9-11EE-A6A1-7958D2021FDD.json'])
    
    