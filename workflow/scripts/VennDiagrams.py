#!/usr/bin/env python3
# coding: utf-8

"""
A script to get the intersections of Differential expression results. Draws Venn diagrams. 
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import itertools
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3


direction = snakemake.wildcards['dir_']
pval = snakemake.params['padj_threshold']
dataset = snakemake.params['dataset']
comparisons = snakemake.params['comparisons'] #['BusiaParental_BusiaSelected', 'Kisumu_BusiaParental', 'Kisumu_BusiaSelected'] 

assert len(comparisons) > 1, "Only one differential expression comparison is specified, cannot run venn analysis. Please disable venn in config.yaml or specify more comparisons"

de_data = {}
de_genes = {}
for comp in comparisons:
    n_comps = len(comparisons)
    
    de_data[comp] = pd.read_csv(f"results/genediff/{comp}.csv")
    de_genes[comp] = de_data[comp].query(f"padj < {pval}")
    if direction == 'up':
        de_genes[comp] = de_genes[comp].query("FC > 1")['GeneID'].to_numpy()
    elif direction == 'down':
        de_genes[comp] = de_genes[comp].query("FC < 1")['GeneID'].to_numpy()

for comp1, comp2 in itertools.combinations(comparisons, 2):
    all_de_comps = [comp1, comp2]
    all_de_comps_string = '.'.join(all_de_comps)
    all_de_comps = [a.replace("_", " v ") for a in all_de_comps]

    de1 = set(de_genes[comp1])
    de2 = set(de_genes[comp2])

    s = [de1, de2]
    venn2(s, all_de_comps)
    plt.title(f"Venn - {all_de_comps_string} - {direction}")
    plt.savefig(f"results/genediff/venn/{all_de_comps_string}-{direction}-Venn.png")
    plt.close()

if n_comps >= 3:
    for all_de_comps in itertools.combinations(comparisons, 3):
        comp1, comp2, comp3 = all_de_comps
        all_de_comps_string = '.'.join(all_de_comps)
        all_de_comps = [a.replace("_", " v ") for a in all_de_comps]

        de_genes2 = dict((k, de_genes[k]) for k in (comp1, comp2, comp3))

        s = [set(v) for v in de_genes2.values()]

        venn3(s, all_de_comps)
        plt.title(f"Venn - {all_de_comps_string} - {direction}")
        plt.savefig(f"results/genediff/venn/{all_de_comps_string}-{direction}-Venn.png")
        plt.close()