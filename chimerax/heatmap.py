import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import os
import sys

# py heatmap.py <input_alignment_file> <output_directory>
# e.g. py heatmap.py C:\XResearch\Matchmaking\selected\alignment.xlsx C:\XResearch\Matchmaking\selected

def heatmap(dataframe, title = None, rownames = None, colnames = None, outdir = None, row_colors = None):
    full_matrix = pd.DataFrame(np.tril(dataframe) + np.tril(dataframe, -1).T)
    full_matrix.index = rownames
    full_matrix.columns = colnames    
    # UPGMA
    sns.clustermap(full_matrix, method="average", col_cluster=False, cbar_pos=(0.11, 0.84, 0.03, 0.15),
                #  annot=True, fmt = "g", annot_kws={"size": 150/full_matrix.shape[0]}, 
                   figsize=(13, 10), cmap='YlOrRd', row_colors=row_colors)
    
    # write file
    plot_name = f'{outdir}/{title}.png'
    if os.path.exists(plot_name):
        os.remove(plot_name)
    plt.savefig(plot_name)
    
    
def handle_excel_file(inputfile, outdir):
    excel_data = pd.read_excel(inputfile, sheet_name=None, header=0, index_col=0)
    
    for sheet_name, sheet_data in excel_data.items():
        dist_matrix = sheet_data.iloc[:,:-1]
        clusters = sheet_data.iloc[:,-1]
        print(clusters)
        
        clusters = pd.read_csv("C:\XResearch\Archive_PDB\selected\clusters.csv", header=None, index_col=0).iloc[:,-1]
        print(clusters)
        palette = 'wrgbymck' if len(clusters.unique()) <= 8 else random.shuffle(sns.color_palette("Paired"))
        
        colors = dict(zip(clusters.unique(), palette))
        print(colors)
        row_colors = clusters.map(colors)
        print(row_colors)
        
        heatmap(dist_matrix, title=sheet_name, rownames=sheet_data.index, colnames=sheet_data.columns[:-1],
                row_colors = row_colors, 
                outdir=outdir)


def main():
    if (len(sys.argv) != 3):
        sys.exit("Input file and output folder are required")
    inputfile = sys.argv[1]
    outdir = os.path.join(sys.argv[2], "heatmap")
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    
    handle_excel_file(inputfile, outdir)
        
if __name__ == "__main__":
    main()

