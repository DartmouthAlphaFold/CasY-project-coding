import sys
import os
import pandas as pd
from bs4 import BeautifulSoup

# cd C:/XResearch
# runscript matchmaker.py <input folder> <output folder>
# e.g. runscript matchmaker.py "C:\XResearch\Proteins\CasY-Predictions-Rank0" "C:\XResearch\Matchmaking\CasYvsCas12c" 


def get_proteins(input_dir, chosen_proteins = None):
    '''
    Get the list of proteins from protein files
    '''
    files = [f for f in os.listdir(input_dir)]
    if chosen_proteins is not None:
        true_list = []
        for chosen in chosen_proteins:
            if chosen in files:
                true_list.append(chosen)
        return true_list
        
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
    if verbose: print(f"Nunber of pruned atom pairs: {pruned_atom_pairs_count}")
    
    pruned_rmsd_index = raw_result.find("RMSD between ")
    pruned_rmsd_info = raw_result[pruned_rmsd_index+len("RMSD between "):].split(";")[0].strip()
    pruned_rmsd = float(pruned_rmsd_info.split(" ")[-2])
    if verbose: print(f"RMSD between pruned atom pairs: {pruned_rmsd}")
    
    atom_pairs_index = raw_result.find("across all ")
    atom_pairs_info = raw_result[atom_pairs_index + len("across all "):].split(";")[0].strip().split(")")[0]
    all_atom_pairs_count = int(atom_pairs_info.split(" ")[0])
    if verbose: print(f"Nunber of all atom pairs: {all_atom_pairs_count}")
    
    all_rmsd = float(atom_pairs_info.split(" ")[-1])
    if verbose: print(f"RMSD between all atom pairs: {all_rmsd}")
    
    return {"alignment_score": alignment_score,
            "pruned_atom_pairs": pruned_atom_pairs_count,
            "pruned_rmsd": pruned_rmsd,
            "all_atom_pairs": all_atom_pairs_count,
            "all_rmsd": all_rmsd}


def match(proteins, input_dir, output_dir):
    '''
    Match proteins, print output to HTML format, then read HTML file to create matrices
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
            run(session, 'match #1 to #2')
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
            
    # write to excel        
    writer = pd.ExcelWriter(f'{output_dir}/alignment.xlsx', engine="xlsxwriter")
    alignment_score_mat.to_excel(writer, sheet_name="Alignment score")
    all_rmsd_mat.to_excel(writer, sheet_name="RMSD (All pairs)")
    all_atom_pairs_mat.to_excel(writer, sheet_name="All atom pairs")
    pruned_rmsd_mat.to_excel(writer, sheet_name="RMSD between pruned atom pairs")
    pruned_atom_pairs_mat.to_excel(writer, "Pruned atom pairs")   
    pruned_over_all_atom_mat.to_excel(writer, "Pruned over All ratio")
    writer.close()     
    

def main():
    if (len(sys.argv) < 3):
        sys.exit("Input and output folders are required")
    
    input_dir = str(sys.argv[1])
    output_dir = str(sys.argv[2])
    
    if (len(sys.argv) > 3):
        chosen_proteins = list(map(str, sys.argv[3:]))
        proteins = get_proteins(input_dir, chosen_proteins)
    else:
        proteins = get_proteins(input_dir)
    
    match(proteins, input_dir, output_dir)

main()