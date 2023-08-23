import sys
import os
import random
import numpy as np
import pandas as pd
import time
from functools import wraps
import pvclust
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
import seaborn as sns

# cd <folder of thie script>
# runscript matchmaker.py <input folder> <output folder>
# e.g. runscript matchmaker.py "C:\XResearch\Archive_PDB\selected" "C:\XResearch\Matchmaking\new_selected" 

def timer(progress_func):
    '''
    Measure how much time each function executed
    '''
    @wraps(progress_func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = progress_func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f'Step "{progress_func.__name__}" took {execution_time:.3f} seconds')
        return result
    return wrapper


def get_proteins(input_dir):
    '''
    Get the list of protein filess from input folder
    '''
    files = [f for f in os.listdir(input_dir) if f.endswith(('.pdb', '.cif'))]
    return files


def read_html_file(input_file: str, verbose = False):
    '''
    Read one HTML file and return a result map
    '''
    
    # get all text
    file = open(f"{input_file}")
    soup = BeautifulSoup(file.read(), "html.parser")
    contents = soup.contents[0]
    contents_text = contents.get_text()
    sentences = contents_text.split("\n\n")
    
    # find the sentence containing results
    raw_result = None
    for sentence in sentences:
        if "Matchmaker" in sentence:
            raw_result = sentence
            break
    if raw_result == None:
        return None
    
    # map results
    alignment_score_index = raw_result.find("score = ")
    alignment_score = float(raw_result[alignment_score_index+len("score = "):].split("\n")[0])
    if verbose: print(f"Sequence alignment score: {alignment_score}")
    
    pairs_index = raw_result.find("RMSD between ")
    pruned_atom_pairs_count = int(raw_result[pairs_index+len("RMSD between "):].split(" ")[0])
    if verbose: print(f"Number of pruned atom pairs: {pruned_atom_pairs_count}")
    
    pruned_rmsd_index = raw_result.find("RMSD between ")
    pruned_rmsd_info = raw_result[pruned_rmsd_index+len("RMSD between "):].split(";")[0].strip()
    pruned_rmsd = float(pruned_rmsd_info.split(" ")[-2])
    if verbose: print(f"RMSD between pruned atom pairs: {pruned_rmsd}")
    
    atom_pairs_index = raw_result.find("across all ")
    atom_pairs_info = raw_result[atom_pairs_index + len("across all "):].split(";")[0].strip().split(")")[0]
    all_atom_pairs_count = int(atom_pairs_info.split(" ")[0])
    if verbose: print(f"Number of all atom pairs: {all_atom_pairs_count}")
    
    all_rmsd = float(atom_pairs_info.split(" ")[-1])
    if verbose: print(f"RMSD between all atom pairs: {all_rmsd}")
    
    return {"alignment_score": alignment_score,
            "pruned_atom_pairs": pruned_atom_pairs_count,
            "pruned_rmsd": pruned_rmsd,
            "all_atom_pairs": all_atom_pairs_count,
            "all_rmsd": all_rmsd}


@timer
def matchmaker(proteins, input_dir, output_dir, ssFraction):
    '''
    Match proteins, print output to HTML format, then read HTML file to print out matrices
    Return 
    '''
    protein_names = list(map(lambda x: x.split(".")[0], proteins))
    
    alignment_score_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    pruned_atom_pairs_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    pruned_rmsd_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    all_atom_pairs_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    all_rmsd_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    pruned_over_all_atom_mat = pd.DataFrame(index = protein_names, columns = protein_names)
    
    from chimerax.core.commands import run
    run(session, 'close session; log clear')
    for i in range(len(proteins)):
        for j in range(i, len(proteins)):
            # match and print to html
            run(session, f'open {input_dir}\{proteins[i]} {input_dir}\{proteins[j]}')
            run(session, f'match #1 to #2 ssFraction {ssFraction}')
            html_file = f'{output_dir}\{os.path.splitext(proteins[i])[0]}.vs.{os.path.splitext(proteins[j])[0]}.html'
            run(session, f'log save {html_file}')
            run(session, 'close session; log clear')
            
            # read html files
            results_map = read_html_file(html_file)
            alignment_score = results_map['alignment_score']
            pruned_atom_pairs = results_map['pruned_atom_pairs']
            pruned_rmsd = results_map['pruned_rmsd']
            all_atom_pairs = results_map['all_atom_pairs']
            all_rmsd = results_map['all_rmsd']
            
            # write into matrices
            protein1_name = proteins[j].split(".")[0]
            protein2_name = proteins[i].split(".")[0]
            alignment_score_mat.loc[protein1_name, protein2_name] = alignment_score
            pruned_atom_pairs_mat.loc[protein1_name, protein2_name] = pruned_atom_pairs
            pruned_rmsd_mat.loc[protein1_name, protein2_name] = pruned_rmsd
            all_atom_pairs_mat.loc[protein1_name, protein2_name] = all_atom_pairs
            all_rmsd_mat.loc[protein1_name, protein2_name] = all_rmsd
            pruned_over_all_atom_mat.loc[protein1_name, protein2_name] = round(pruned_atom_pairs/all_atom_pairs, 3)
            
            # remove the html file
            if os.path.exists(html_file):
                os.remove(html_file) 
    
    # normalize the alignment score to itself
    norm_alignment_score = alignment_score_mat / np.diagonal(alignment_score_mat)[:, np.newaxis]
    norm_alignment_score_mat = pd.DataFrame(norm_alignment_score, index = protein_names, columns = protein_names)
    
    # normalize the alignment score by length
    norm_alignment_score_2 = alignment_score_mat / all_atom_pairs_mat
    norm_alignment_score_mat_2 = pd.DataFrame(norm_alignment_score_2, index = protein_names, columns = protein_names)
    
    # write to excel    
    matrix_file = f'{output_dir}/alignment.xlsx' 
    writer = pd.ExcelWriter(matrix_file, engine="xlsxwriter")
    alignment_score_mat.to_excel(writer, sheet_name="Alignment score")
    norm_alignment_score_mat.to_excel(writer, sheet_name="Norm alignment score diag")
    norm_alignment_score_mat_2.to_excel(writer, sheet_name="Norm alignment score length")
    all_rmsd_mat.to_excel(writer, sheet_name="RMSD (All pairs)")
    all_atom_pairs_mat.to_excel(writer, sheet_name="All atom pairs")
    pruned_rmsd_mat.to_excel(writer, sheet_name="RMSD between pruned atom pairs")
    pruned_atom_pairs_mat.to_excel(writer, "Pruned atom pairs")   
    pruned_over_all_atom_mat.to_excel(writer, "Pruned over All ratio")
    writer.close()    
     
    return matrix_file
    

def clustermap(matrix_file, output_dir, clusters=None):
    '''
    Plot heatmaps and hierarchical clustering tree for each distance matrices
    using the `seaborn` library
    '''
    excel_data = pd.read_excel(matrix_file, sheet_name=None, header=0, index_col=0)
    
    if len(clusters) == 0:
        row_colors = None
    else:
        palette = 'wrgbymck' if len(clusters.unique()) <= 8 else random.shuffle(sns.color_palette("Paired"))
        colors = dict(zip(clusters.unique(), palette))
        row_colors = clusters.map(colors)
    
    for sheet_name, sheet_data in excel_data.items():
        full_matrix = pd.DataFrame(np.tril(sheet_data) + np.tril(sheet_data, -1).T)
        full_matrix.index = sheet_data.index
        full_matrix.columns = sheet_data.columns    
            
        # UPGMA
        sns.clustermap(full_matrix, method="average", col_cluster=False, cbar_pos=(0.11, 0.84, 0.03, 0.15),
                #  annot=True, fmt = "g", annot_kws={"size": 150/full_matrix.shape[0]}, 
                   figsize=(13, 10), cmap='YlOrRd', row_colors=row_colors)
        
        # write file
        plot_name = f'{output_dir}/{sheet_name}.png'
        plt.savefig(plot_name)


@timer
def pvclustR(script_path, matrix_file, output_dir, clusters = None):
    '''
    Perform hierarchical clustering via multiscale bootstrap resampling using 
    an R script that calls the `pvclust` library
    '''
    import subprocess
    try:
        subprocess.run(["C:/Program Files/R/R-4.2.2/bin/Rscript", script_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing the R script: {e}")


@timer
def pvclustpy(matrix_file, output_dir, clusters = None):
    excel_data = pd.read_excel(matrix_file, sheet_name=None, header=0, index_col=0)
    
    random.seed(0)
    for sheet_name, sheet_data in excel_data.items():
        full_matrix = pd.DataFrame(np.tril(sheet_data) + np.tril(sheet_data, -1).T)
        full_matrix.index = sheet_data.index
        full_matrix.columns = sheet_data.columns  
        pv = pvclust.PvClust(full_matrix, method="ward", metric="euclidean", nboot=10000, parallel=False)
        pv.plot(f"{output_dir}/Dendrogram ({sheet_name}).pdf", labels=full_matrix.index,
                param_display='AU', sig_level = 95, orientation = "left")
        print(f"\n...{sheet_name}...")
        pv.print_result(digits=5)
        pv.seplot(f"{output_dir}/SEplot ({sheet_name}).pdf", annotate=True)


def main():
    if (len(sys.argv) < 3):
        sys.exit("Input and output folders are required")
    
    input_dir = str(sys.argv[1])
    output_dir = str(sys.argv[2])
    
    # distance matrix
    proteins = get_proteins(input_dir)
    ssFraction = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    matrix_file = matchmaker(proteins, input_dir, output_dir, ssFraction)
    
    # build dendrogram, colored by provided cluster map
    clusters_file = f'{input_dir}/clusters.csv'
    if os.path.exists(clusters_file):
        clusters = pd.read_csv(clusters_file, header=0, index_col=0).iloc[:,-1]
    else:
        clusters = pd.Series()
    # clustermap(matrix_file, output_dir, clusters)
    pvclustpy(matrix_file, output_dir, clusters=None)


def test():
    input_dir = "C:/XResearch/Archive_PDB/selected"
    output_dir = "C:/XResearch/Matchmaking/selected/pvclust_parallel"
    matrix_file = "C:/XResearch/Matchmaking/selected/alignment.xlsx"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # script_path = "C:/XResearch/Coding/chimerax/pvclust.R"
    # pvclust(script_path, matrix_file, output_dir)
    pvclustpy(matrix_file, output_dir)
    
    
print(f'Running in {__name__}')
if __name__.startswith('ChimeraX_sandbox'):
    main()