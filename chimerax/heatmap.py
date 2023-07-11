import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from scipy.spatial import distance

# py heatmap.py <input_alignment_file> <output_directory>
# e.g. py heatmap.py C:\XResearch\Matchmaking\CasYvsCas12c\alignment.xlsx C:\XResearch\Matchmaking\CasYvsCas12c

def heatmap(dataframe, title = None, rownames = None, colnames = None, outdir = None):
    full_matrix = pd.DataFrame(np.tril(dataframe) + np.tril(dataframe, -1).T)
    full_matrix.index = rownames
    full_matrix.columns = colnames    
    # UPGMA
    sns.clustermap(full_matrix, annot=True, method="average", col_cluster=False,
                   cmap='YlOrRd', fmt = "g", annot_kws={"size": 150/full_matrix.shape[0]}, figsize=(13, 10))
    plt.title(title)
    
    # write file
    plot_name = f'{outdir}/{title}.png'
    if os.path.exists(plot_name):
        os.remove(plot_name)
    plt.savefig(plot_name)
    
    
def handle_excel_file(inputfile, outdir):
    excel_data = pd.read_excel(inputfile, sheet_name=None, header=0, index_col=0)
    for sheet_name, sheet_data in excel_data.items():
        heatmap(sheet_data, title=sheet_name, rownames=sheet_data.index, colnames=sheet_data.columns, outdir=outdir)


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

