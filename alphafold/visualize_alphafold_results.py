import os
import glob
import pickle
import json
import numpy as np
import matplotlib.pyplot as plt
import sys
 
class ARG:
    def __init__(self, repo):
        self.input_dir = repo
        self.output_dir = repo
        self.name = os.path.basename(os.path.normpath(repo))
        
def get_pae_plddt(model_dicts, model_type):
    out = {}
    if (model_type == "monomer"):
        for i,d in enumerate(model_dicts):
            out[f'model_{i+1}'] = {'plddt': d['plddt']}
    else:    
        for i,d in enumerate(model_dicts):
            out[f'model_{i+1}'] = {'plddt': d['plddt'], 
                                'pae':d['predicted_aligned_error']}
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
        plt.plot(value["plddt"], 
                 label=f"{model_name}, plddts: {round(list(ranking_dict['plddts'].values())[s], 6)}")
        s += 1
        
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}LDDT.png")
    
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
    repo = [os.path.join(main_repo, entry.name) for entry in os.scandir(main_repo) if entry.is_dir()] # This is a list of all output repositories
    
    for r in repo:
        args = ARG(r)
        with open(os.path.join(r, "ranking_debug.json"), 'r') as f:
            ranking_dict = json.load(f)
            
        feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))
    
        is_multimer = ('result_model_1_multimer_v2_pred_0.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
        is_monomer_ptm = ('result_model_1_ptm_pred_0.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
        
        if is_multimer==False:
            model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}{"_multimer_v2" if is_multimer else ""}{"_ptm" if is_monomer_ptm else ""}_pred_0.pkl','rb'))
                        for f in range(1,6)]
        else:
            model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}{"_multimer_v2" if is_multimer else ""}{"_ptm" if is_monomer_ptm else ""}_pred_{g}.pkl','rb')) 
                        for f in range(1,6) for g in range(5)]

        model_type = "multimer" if is_multimer else ("monomer_ptm" if is_monomer_ptm else "monomer")
        pae_plddt_per_model = get_pae_plddt(model_dicts, model_type)
        generate_output_images(feature_dict, model_dicts, ranking_dict, 
                            args.output_dir if args.output_dir else args.input_dir, 
                            args.name, pae_plddt_per_model, model_type)
    
if __name__ == "__main__":
    main()