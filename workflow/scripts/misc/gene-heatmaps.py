#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import numpy as np


def load_data(args):
    normcounts = pd.read_csv(args[1], sep="\t")
    gene_list = pd.read_csv(args[2], sep="\t", header=None).iloc[:,0].to_list()
    out_file = args[3]
    return(normcounts, gene_list, out_file)

def plot_heatmap(normcounts, gene_list, out_file):
    counts_data = normcounts.query("GeneID in @gene_list").set_index("GeneID")
    
    cg = sns.clustermap(data=np.log2(counts_data))
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    print(f"Saving heatmap to {out_file}")
    plt.savefig(out_file)
    return(cg)


# ### Plot clustered heatmap 
normCounts, gene_list, out_file = load_data(sys.argv)
cg = plot_heatmap(normCounts, gene_list, out_file)