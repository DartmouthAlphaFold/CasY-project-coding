import os
import glob
import pickle
import json
import numpy as np
import sys
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})

# text size
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIG_SIZE = 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

# color
LOWEST_CONF_COLOR = '#FF7D45'
LOW_CONF_COLOR = '#FFDB13'
HIGH_CONF_COLOR = '#65CBF3'
HIGHEST_CONF_COLOR = '#0053D6'


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
    # plt.figure(figsize=(14, 7), dpi=300)
    # # plt.subplot(1, 2, 1)
    # plt.title(f"Sequence coverage ({name})")
    # plt.imshow(final,
    #            interpolation='nearest', aspect='auto',
    #            cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    # plt.plot((msa != 21).sum(0), color='black')
    # plt.xlim(-0.5, msa.shape[1] - 0.5)
    # plt.ylim(-0.5, msa.shape[0] - 0.5)
    # plt.colorbar(label="Sequence identity to query", )
    # plt.xlabel("Positions")
    # plt.ylabel("Sequences")
    # plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage.png")
    
    ###################### PLOT pLDDT PER POSITION (colored scatterplot) ####################
    plt.figure(figsize=(15, 5), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        data = value["plddt"]
        
        plt.clf()
        fig, ax = plt.subplots()
        ax.axhline(y=50, linestyle='--', lw = 0.5, color = "black")
        ax.axhline(y=70, linestyle='--', lw = 0.5, color = "black")
        ax.axhline(y=90, linestyle='--', lw = 0.5, color = "black")

        # Color each section of the main line plot
        for i in range(len(data)):
            if data[i] < 50:
                color = '#FF7D45'
            elif data[i] < 70:
                color = '#FFDB13'
            elif data[i] < 90:
                color = '#65CBF3'
            else:
                color = '#0053D6'
            # ax.plot([i, i+1], [data[i], data[i+1]], color=color)
            ax.scatter(i, data[i], color = color, s = 1)
        ax.set_title(f"Predicted LDDT per position \n {name} {model_name} (mean pLDDT = {round(list(ranking_dict['plddts'].values())[n], 6)})")
        ax.set_ylim(0, 100)
        ax.set_yticks((0,50,70,90,100))
        ax.set_ylabel("Predicted LDDT")
        ax.set_xlabel("Positions")
        fig.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_pLDDT_scatter.pdf")
    
    
    ###################### PLOT pLDDT PER POSITION (line plot, dashed line) ####################
    plt.figure(figsize=(15, 5), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.clf()
        fig, ax = plt.subplots()
        # dashed line
        ax.axhline(y=50, linestyle='--', lw = 0.25, color = "black")
        ax.axhline(y=70, linestyle='--', lw = 0.25, color = "black")
        ax.axhline(y=90, linestyle='--', lw = 0.25, color = "black")
        # main plot
        ax.plot(value["plddt"], lw = 0.75)
        ax.set_title(f"Predicted LDDT per position \n {name} {model_name} (mean pLDDT = {round(list(ranking_dict['plddts'].values())[n], 6)})")
        ax.set_ylim(0, 100)
        ax.set_yticks((0,50,70,90,100))
        ax.set_ylabel("Predicted LDDT")
        ax.set_xlabel("Positions")
        fig.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_pLDDT_line.pdf")
    
    
    ###################### PLOT pLDDT PER POSITION (colored bar plot) ####################
    plt.figure(figsize=(15, 5), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        data = value["plddt"]
        plt.clf()
        fig, ax = plt.subplots()
        
        bar_colors = []
        for i in range(len(data)):
            color = ""
            if data[i] < 50:
                color = '#FF7D45'
            elif data[i] < 70:
                color = '#FFDB13'
            elif data[i] < 90:
                color = '#65CBF3'
            else:
                color = '#0053D6'
            bar_colors.append(color)
            
        ax.bar(np.arange(len(data)) + 1, data, color = bar_colors)
        ax.set_title(f"Predicted LDDT per position \n {name} {model_name} (mean pLDDT = {round(list(ranking_dict['plddts'].values())[n], 6)})")
        ax.set_ylim(0, 100)
        ax.set_ylabel("Predicted LDDT")
        ax.set_xlabel("Positions")
        fig.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_pLDDT_bar.pdf")
        
    
    ###################### PLOT pLDDT PER POSITION (line plot with colored dots) ####################
    plt.figure(figsize=(15, 5), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        data = np.array(value["plddt"])[:65]
        
        plt.clf()
        fig, ax = plt.subplots()
        # Line plot
        ax.plot(np.arange(len(data)) + 1, data, linestyle='dotted', lw = 0.75, color = "black")
        # Color each dot
        for i in range(len(data)):
            if data[i] > 90:
                color = HIGHEST_CONF_COLOR
            elif data[i] > 70:
                color = HIGH_CONF_COLOR
            elif data[i] > 50:
                color = LOW_CONF_COLOR
            else:
                color = LOWEST_CONF_COLOR
            ax.scatter(i+1, data[i], color = color, s = int(1500/len(data) + 1))
            
        ax.set_title(f"Predicted LDDT per position \n {name} {model_name} (mean pLDDT = {round(list(ranking_dict['plddts'].values())[n], 6)})")
        ax.set_ylim(0, 100)
        ax.set_yticks((0,50,70,90,100))
        ax.set_ylabel("Predicted LDDT")
        ax.set_xlim(0,)
        ax.set_xlabel("Positions")
        # Legend
        highest_conf_patch = mpatches.Patch(color=HIGHEST_CONF_COLOR, label='Very high (pLDDT > 90)')
        high_conf_patch = mpatches.Patch(color=HIGH_CONF_COLOR, label='Confident (90 > pLDDT > 70)')
        low_conf_patch = mpatches.Patch(color=LOW_CONF_COLOR, label='Low (70 > pLDDT > 50)')
        lowest_conf_patch = mpatches.Patch(color=LOWEST_CONF_COLOR, label='Very low (pLDDT < 50)')
        ax.legend(handles=[highest_conf_patch, high_conf_patch, low_conf_patch, lowest_conf_patch],
                  loc='lower right')
        fig.tight_layout()
        fig.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_pLDDT_scatter_line.pdf")
    
    
    ################# PLOT THE PREDICTED ALIGNED ERROR################
    # if (model_type == "monomer"):
    #     return
    # # num_models = len(model_dicts)
    # plt.figure(figsize=(14, 14), dpi=300)
    # for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
    #     plt.clf()
    #     # plt.subplot(1, num_models, n + 1)
    #     plt.title(f"{name} {model_name} (mean pLDDT = {round(list(ranking_dict['plddts'].values())[n], 6)})")
    #     plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
    #     plt.colorbar()
    #     plt.savefig(f"{out_dir}/{name+('_' if name else '')+model_name}_PAE.png")
        
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