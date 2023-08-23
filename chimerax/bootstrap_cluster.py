import os
import random
import numpy as np
import pandas as pd
import time
from functools import wraps
from pvclust import PvClust
import argparse


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
        print(f"Step {progress_func.__name__} took {execution_time:.3f} seconds")
        return result
    return wrapper

@timer
def pvclustpy(matrix_file, output_dir, clusters = None):
    excel_data = pd.read_excel(matrix_file, sheet_name=None, header=0, index_col=0)
    
    random.seed(0) 
    for sheet_name, sheet_data in excel_data.items():
        full_matrix = pd.DataFrame(np.tril(sheet_data) + np.tril(sheet_data, -1).T)
        full_matrix.index = sheet_data.index
        full_matrix.columns = sheet_data.columns  
        pv = PvClust(full_matrix, method="ward", metric="euclidean", nboot=10000, parallel=False)
        pv.plot(f"{output_dir}/Dendrogram ({sheet_name}).pdf", labels=full_matrix.index,
                param_display='AU', sig_level = 95, orientation = "left")
        print(f"\n...{sheet_name}...")
        pv.print_result(digits=5)
        pv.seplot(f"{output_dir}/SEplot ({sheet_name}).pdf", annotate=True)
    

def main():
    parser = argparse.ArgumentParser(description='A simple CLI tool that builds phylogenetic trees from precomputed distance matrices'
                                     ' using `pvclust` library implemented in Python')
    parser.add_argument('--input', '-i', required=True, help="Input distance matrices in Excel format")
    parser.add_argument('--output', '-o', required=True, help="Output folder")
    parser.add_argument('--metric', '-m', default='ward', help='Metric for hierarchical clustering available at'
                        ' `scipy.cluster.hierarchy.linkage`')
    parser.add_argument('--distance', '-d', default='euclidean', help='Distance for hierarchical clustering available at'
                        ' `scipy.spatial.distance.pdist`')
    parser.add_argument('--sig_level', '-p', default=95, help='Significant level')
    args = parser.parse_args()
    
    # cluster file
    clusters_file = f'{args.input}/clusters.csv'
    if os.path.exists(clusters_file):
        clusters = pd.read_csv(clusters_file, header=0, index_col=0).iloc[:,-1]
    else:
        clusters = pd.Series()
    
    # output folder
    if not os.path.exists(args.output):
        os.mkdir(args.output)
        
    
    
    
def testcolor():
    # Color mapping
    from scipy.cluster.hierarchy import (dendrogram, set_link_color_palette,
                                     leaves_list, to_tree, linkage)
    from scipy.spatial import distance
    from sklearn.datasets import load_diabetes
    import matplotlib.pyplot as plt
    A_data = load_diabetes().data
    DF_diabetes = pd.DataFrame(A_data, columns = ["attr_%d" % j for j in range(A_data.shape[1])])

    # Absolute value of correlation matrix, then subtract from 1 for disimilarity
    DF_dism = 1 - np.abs(DF_diabetes.corr())

    # Compute average linkage
    A_dist = distance.squareform(DF_dism.values)
    Z = linkage(A_dist,method="average")
    dflt_col = "#808080"   # Unclustered gray
    D_leaf_colors = {"attr_1": dflt_col,

                    "attr_4": "#B061FF", # Cluster 1 indigo
                    "attr_5": "#61ffff",
                    "attr_2": "#B061FF",
                    "attr_8": "#B061FF",
                    "attr_6": "#B061FF",
                    "attr_7": "#B061FF",

                    "attr_0": "#61ffff", # Cluster 2 cyan
                    "attr_3": "#61ffff",
                    "attr_9": "#61ffff",
                    }

    # notes:
    # * rows in Z correspond to "inverted U" links that connect clusters
    # * rows are ordered by increasing distance
    # * if the colors of the connected clusters match, use that color for link
    link_cols = {}
    for i, i12 in enumerate(Z[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(Z) else D_leaf_colors["attr_%d"%x] for x in i12)
        link_cols[i+1+len(Z)] = c1 if c1 == c2 else dflt_col

    # Dendrogram
    D = dendrogram(Z=Z, labels=DF_dism.index, color_threshold=None,
    leaf_font_size=12, leaf_rotation=45, link_color_func=lambda x: link_cols[x])
    plt.show()

if __name__ == "__main__":
    testcolor()