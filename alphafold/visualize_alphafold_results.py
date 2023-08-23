#!/usr/bin/env python3

'''
Python script for visualizing AlphaFold 2.3.2 results on Dartmouth HPC cluster
By Duc Nguyen '24, Colby College
Updated on 2023-08-23
'''

import os
import glob
import pickle
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
 
class ARG:
    def __init__(self, repo):
        self.input_dir = repo
        self.output_dir = repo
        self.name = os.path.basename(os.path.normpath(repo))
        
def get_pae_plddt(model_dicts, model_type):
    out = {}
    if model_type == "multimer":
        for i,d in enumerate(model_dicts):
            out[d['model_name']] = {'plddt': d['plddt'],
                                   'pae':d['predicted_aligned_error'],
                                   'ptm':d['ptm'],
                                   'iptm':d['iptm']}
    elif model_type == "monomer_ptm":    
        for i,d in enumerate(model_dicts):
            out[f'model_{i+1}'] = {'plddt': d['plddt'], 
                                'pae':d['predicted_aligned_error'],
                                'ptm':d['ptm']}
    else:
        for i,d in enumerate(model_dicts):
            out[f'model_{i+1}'] = {'plddt': d['plddt']}
    return out
 
def generate_output_images(feature_dict, model_dicts, ranking_dict, 
                           out_dir, name, pae_plddt_per_model, model_type):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
    
    ###################### PLOT MSA WITH COVERAGE ####################
    plt.figure(figsize=(14, 7), dpi=300)
    # plt.subplot(1, 2, 1)
    plt.title(f"Sequence coverage ({name})")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage.png")
    
    ###################### PLOT LDDT PER POSITION ####################
    plt.figure(figsize=(14, 7), dpi=300)
    # plt.subplot(1, 2, 2)
    plt.title(f"Predicted LDDT per position ({name})")
    
    s = 0
    for model_name, value in pae_plddt_per_model.items():
        if model_type == 'multimer':
            plt.plot(value["plddt"],
                    label=f"{model_name}, iptm+ptm: {round(list(ranking_dict['iptm+ptm'].values())[s], 6)}, plddt: {round(np.mean(value['plddt']), 6)}")
        elif model_type == 'monomer_ptm':
            plt.plot(value["plddt"], 
                    label=f"{model_name}, plddt: {round(list(ranking_dict['plddts'].values())[s], 6)}, ptm: {round(float(value['ptm']), 6)}")
        else:
            plt.plot(value["plddt"], 
                    label=f"{model_name}, plddt: {round(list(ranking_dict['plddts'].values())[s], 6)}")
        s += 1
        
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}LDDT.png")
    
    ###################### EXPORT LDDT PER POSITION PER MODEL ####################
    plddt_df = pd.DataFrame(columns = ['position'] + list(pae_plddt_per_model.keys()))
    for model_name, value in pae_plddt_per_model.items():
        plddt_df[model_name] = value["plddt"]
    plddt_df['position'] = np.arange(1, plddt_df.shape[0] + 1)
    plddt_df.to_csv(f"{out_dir}/plddts.csv", index = False, mode = 'w')    
    
    ################# PLOT THE PREDICTED ALIGNED ERROR################
    if (model_type == "monomer"):
        return
    num_models = len(model_dicts)
    plt.figure(figsize=(14, 14), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.clf()
        # plt.subplot(1, num_models, n + 1)
        plt.title(f"{name} {model_name}")
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
        plt.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_PAE.png")
        
    # plt.savefig(f"{out_dir}/{name+('_' if name else '')}PAE.pdf")
    ##################################################################

def main():
    
    if (len(sys.argv) == 1):
        sys.exit("Output/Input directory are required")
    
    main_repo = sys.argv[1]
    repo = [os.path.join(main_repo, entry.name) for entry in os.scandir(main_repo) if entry.is_dir() and entry.name!='msas'] # This is a list of all output repositories
    if len(repo) == 0:
        sys.exit(f'Found 0 folder containing AlphaFold results inside {main_repo}.\nPlease check the path again.')
    
    
    for r in repo:
        args = ARG(r)
        with open(os.path.join(r, "ranking_debug.json"), 'r') as f:
            ranking_dict = json.load(f)
            
        feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))

        file_list = [os.path.basename(f) for f in os.listdir(args.input_dir)]
        is_multimer = len(set(['result_model_1_multimer_pred_0.pkl', 'result_model_1_multimer_v2_pred_0.pkl', 'result_model_1_multimer_v3_pred_0.pkl']) & set(file_list)) > 0
        is_monomer_ptm = 'result_model_1_ptm_pred_0.pkl' in file_list
        model_type = "multimer" if is_multimer else ("monomer_ptm" if is_monomer_ptm else "monomer")
        
        if model_type != 'multimer':
            model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}{"_ptm" if is_monomer_ptm else ""}_pred_0.pkl','rb'))
                            for f in range(1,6)] 
        else:
            pkl_files = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir) if file.startswith('result_model_') and file.endswith('.pkl')]
            model_dicts = []
            for pkl_file in sorted(pkl_files):
                model_dict = pickle.load(open(pkl_file, 'rb'))
                model_dict['model_name'] = os.path.splitext(os.path.basename(pkl_file))[0].replace('result_', '')
                model_dicts.append(model_dict)
        

        pae_plddt_per_model = get_pae_plddt(model_dicts, model_type)
        generate_output_images(feature_dict, model_dicts, ranking_dict, 
                            args.output_dir if args.output_dir else args.input_dir, 
                            args.name, pae_plddt_per_model, model_type)
    
if __name__ == "__main__":
    main()