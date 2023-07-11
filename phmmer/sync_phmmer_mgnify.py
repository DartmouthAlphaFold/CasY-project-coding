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
import logging
import hashlib
import threading


def get_sequence(args, input_files):
    '''
    Fetch all sequences from a list of FASTA files
    '''
    seq_map = {}
    file_map = {}
    for input_file in input_files:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            if len(seq_record.seq) == 0:
                continue
            sequence = re.sub('[\*\-]', '', str(seq_record.seq))  # remove '*' and '-' characters
            seq_map[str(seq_record.id)] = sequence.upper()
            file_map[str(seq_record.id)] = input_file
            
    if args.debug: print(seq_map.items())
    print(f"\nReading {len(seq_map)} sequences: \n {', '.join(seq_map.keys())}")
    return seq_map, file_map


# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers
    
    
def phmmer(args, name, sequence):
    '''
    Send API to search protein sequence against databases using pHMMER
    Receive and export subjects found into a metadata file
    '''
    
    # get params
    sequence_db=args.database
    range=args.range
    output_dir=args.output
    
    opener = build_opener(SmartRedirectHandler())
    install_opener(opener)

    parameters = {
        'seq': sequence,
        'seqdb': sequence_db
    }
    enc_params = urlencode(parameters).encode("utf-8")

    # post the search request to the server
    request = Request('https://www.ebi.ac.uk/metagenomics/sequence-search/search/phmmer', enc_params)
    try:
        results = urlopen(request)
    except Exception as e:
        print(f"Failed to search for {name}: {e.msg}")
        return
        
    results_url = results.get('location')
    session_id = urlparse(results_url).path.split("/")[-2]

    # modify the range in results
    res_params = {
        'output': 'json',
        'range': f'1,{str(range)}'
    }

    # add the parameters to your request for the results
    enc_res_params = urlencode(res_params)
    modified_res_url = results_url + '?' + enc_res_params
    # https://www.ebi.ac.uk/metagenomics/sequence-search/results/{uuid}/score
    
    # send a GET request to the server
    results_request = Request(modified_res_url)
    try:
        data = urlopen(results_request)
    except Exception as e:
        print(f"Failed to get {name} results: {e.msg}")
        return

    # write the results
    if (data.status == 200):
        output_file = f'{output_dir}/metadata/{name}.{str(len(sequence))}.{session_id}.json'
        with open(output_file, 'w') as f:
            f.write(data.read().decode('utf-8'))
    else:
        print(f'{data.status}\n{data.msg}\n{data.reason}\n{data.url}')
        
    
def thread_phmmer(args, seq_map):
    '''
    Create threads and search for sequences simultaneously
    '''
    threads = []  # list to hold all the threads
    for name, sequence in seq_map.items():
        # create a thread for each query
        thread = threading.Thread(target=phmmer, args=(args, name, sequence,))
        thread.start()
        threads.append(thread)
    
    print(f'\nPerforming {len(threads)} searches:')
    # wait for all threads to finish
    for thread in threads:
        thread.join()
    

def read_json(args):
    '''
    Read the metadata in JSON format and return a map of unique sequence names
    '''
    json_files = [file for file in os.listdir(os.path.join(args.output,'metadata')) if file.endswith('.json')]
    dict_by_name = {}  # dictionary that holds all unique sequence names
    dict_with_id = {}
    print(f'\nReading {len(json_files)} search results:')
    
    for file in json_files:
        data = file.split('.')
        query_name = data[0]
        query_length = int(data[1])
        
        with open(os.path.join(args.output, 'metadata', file), 'r') as f:
            try:
                results = json.load(f)['results']
            except Exception as e:
                print(f"Cannot read results for {query_name}. Skipping to next file.")
                print(e)
                print("\n")
                continue
            num_hits = results['stats']['nhits']
            
            if (num_hits > 0):
                sig_num_hits = 0
                # downloaded_num_hits = 0
                
                for hit in results['hits']:
                    # skip insignificant hits
                    if (float(hit['evalue']) >= args.evalue):
                        break
                    
                    # skip identical hit (identical hit only has 1 domain with 100% identity)
                    if (len(hit['domains']) == 1 and hit['domains'][0]['aliId'] == 1.0):
                        continue
                    
                    # skip hits that have low domain coverage of the query
                    # domains_coverage = hit['domains'][-1]['alisqto'] - hit['domains'][0]['alisqfrom']
                    # if domains_coverage <= 0.5*query_length:
                    #     continue
                    
                    sig_num_hits += 1
                    if hit['name'] in dict_by_name:
                        dict_by_name[hit['name']] += 1
                    else:
                        dict_by_name[hit['name']] = 1  
                        dict_with_id[hit['name']] = hit['sindex']           
                    
                    if args.debug:
                        print(f"{hit['name']}\n{hit['species']}\n{hit['desc']}\n{hit['evalue']}\n"
                        f"{hit['score']}\n{hit['extlink']}\n") 

                print(f"{query_name}: {sig_num_hits}/{num_hits} hits")
            else:
                print(f"{query_name}: no hits were found")
                
    print(f'\nFound {len(dict_by_name)} hits in total:')
    for key, value in dict_by_name.items():
        print(f'{key}: {value} times') 
    return dict_with_id


def download_seq (args, seq_id, seq_ac):
    '''
    Send API to download the sequence to a FASTA file
    '''
    # build the URL
    hit_seq_url = f"https://www.ebi.ac.uk/metagenomics/sequence-search/seq?seq_ac={seq_ac}&seq_id={seq_id}"
    
    # send a GET request to obtain the full sequence
    hit_seq_request = Request(hit_seq_url)
    try:
        hit_seq_response = urlopen(hit_seq_request)
    except Exception as e:
        print(f"Failed to download {seq_ac}: {e.msg}")
        return
    
    # write sequence
    if (hit_seq_response.status == 200):
        seq_output = f"{args.output}/fasta/{seq_ac}.fa"
        
        with open(seq_output, 'w') as f:
            # soup = BeautifulSoup(hit_seq_response.read().decode('utf-8'), 'html.parser')
            # content = soup.find('pre').string
            content = hit_seq_response.read().decode('utf-8')
            content = re.sub(r'<pre.*?>|</pre>', '', content)
            f.write(content)
        print(f"\t{seq_ac} downloaded")
    else:
        if args.debug: print(f'{hit_seq_response.status}\n{hit_seq_response.msg}\n{hit_seq_response.reason}\n{hit_seq_response.url}\n')
        
    
def thread_download_seq (args, dict_with_id):
    '''
    Create threads and download sequences simultaneously
    '''
    threads = []  # list to hold all the threads
    for seq_ac, seq_id in dict_with_id.items():
        # create a thread for each query
        thread = threading.Thread(target=download_seq, args=(args, seq_id, seq_ac))
        thread.start()
        threads.append(thread)
    
    print(f"\nDownloading {len(threads)} sequences:")
    # wait for all threads to finish
    for thread in threads:
        thread.join()


def filter_duplicates(args):
    '''
    Filter duplicated sequences among subjects and between queries and subjects
    Export into a FASTA file 
    '''
    # combine queries and subjects into a map
    query_files = args.input.split(",")
    subject_folder = os.path.join(args.output, 'fasta')
    subject_files = [os.path.join(subject_folder, file) for file in os.listdir(subject_folder) if file.endswith('.fa')]
    seq_map, file_map = get_sequence(args, query_files + subject_files)
    print(f"\nComparing {len(seq_map)} sequences:")
    
    hash_dict = {}  # dictionary to store hash values
    # generate hash values for each sequence
    for name, sequence in seq_map.items():
        sha256_hash = hashlib.sha256()
        sha256_hash.update(sequence.encode('utf-8'))
        hash_value = sha256_hash.hexdigest()
        
        if hash_value in hash_dict:
            hash_dict[hash_value].append(name)
        else:
            hash_dict[hash_value] = [name]

    # create the result file
    result_file = os.path.join(args.output, 'results.fa')
    with open(result_file, 'w'): pass
    
    count = 0
    # compare sequences with the same hash value
    for hash_value, sequence_list in hash_dict.items():
        if len(sequence_list) > 1:
            print(f"{', '.join(sequence_list)} are identical")

        # write the 1st sequence from each hash value, skip if it is one of the queries
        first_id = sequence_list[0]
        if (file_map[first_id] in query_files):
            continue
        
        print(f'\nWriting {first_id}\n')
        valid_seq_file = file_map[first_id]
        with open(valid_seq_file, 'r') as input:
            content = input.read()
        with open(result_file, 'a') as output:
            output.write(content)
        count += 1
    print(f"Finished writing {count} sequences")
            

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
    parser = argparse.ArgumentParser(description="A simple CLI tool for pHMMER-based sequence searching on MGnify database, "
                                     "using the API available at https://www.ebi.ac.uk/metagenomics/sequence-search/search/phmmer"
                                     ".\nCurrent database version is 2018_06.")
    parser.add_argument('--input', '-i', required=True,
                        help="Files in FASTA format containing the query protein sequence(s). Each file can contain multiple sequences.\n"
                        "Files should be separated by commas. All FASTA sequences must have unique names.")
    parser.add_argument('--output', '-o', required=True,
                        help="Output folder")
    parser.add_argument('--evalue', '-e', default = 0.01, type = float,
                        help = "Threshold of e-value. (default: 0.01)")
    parser.add_argument('--range', '-r', default=50, type = int,
                        help="Maximum number of top subjects for each query (default: 50)")
    parser.add_argument('--database', '-db', default='full',
                        choices=['all', 'partial', 'full', 
                                 'aquatic', 'marine', 'freshwater', 'soil',
                                 'human', 'human-digestive', 'human-non-digestive','animal', 
                                 'engineered', 'other'],
                        help="Sequence type: all, partial, or full length sequences."
                        "Environment: aquatic, marine, freshwater, or soil."
                        "Host-associated biome: human, human disgestive system, human non-digestive, or animal."
                        "Other biomes/environments: engineered or other. (default: full)")
    parser.add_argument('--debug', type=str_to_bool, default = False,
                        help="Print detailed info for debugging. (default: False)")
    args = parser.parse_args()
    
    # make output directories
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    if not os.path.exists(os.path.join(args.output, 'metadata')):
        os.mkdir(os.path.join(args.output, 'metadata'))
    if not os.path.exists(os.path.join(args.output, 'fasta')):
        os.mkdir(os.path.join(args.output, 'fasta'))
    
    inputs = args.input.split(",")  # get input files
    
    start_time = time.time()
    seq_map, inputfile_map = get_sequence(args, inputs)
    thread_phmmer(args, seq_map)
    download_seq_ids = read_json(args)
    thread_download_seq(args, download_seq_ids)
    filter_duplicates(args)
    end_time = time.time()
    print(f"Finished in {round(end_time - start_time, 3)} seconds")
    
    
if __name__=="__main__":
    main()
    
    # args = argparse.ArgumentParser().parse_args()
    # args.evalue = 0.01
    # args.debug = False
    # args.database = 'full'
    # args.output = "./phmmer/mgnify_sync_Cas12cd_test_v2"
    # args.input = "./Amino_Acids/FASTA/Cas12cd.fa"
    # seq_id_list = read_json(args)
    # thread_download_seq(args, seq_id_list)
    # filter_duplicates(args)
    
    