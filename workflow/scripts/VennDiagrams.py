#!/usr/bin/env python3

"""
A script to get the intersections of Differential expression results, Fst, and differential SNPs analysis.
Draws Venn diagrams and adds columns to RNA-seq-diff.xlsx, whether the gene has high Fst/PBS/diffsnps. 
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib_venn import *
import pandas as pd
import numpy as np
from pathlib import Path

def plotvenn2(name, group1, group2, nboth,stat="DE_PBS", group1name='Significant up DE genes', group2name='High PBS'):
    
    print(f"There are {group2} high Fst genes in {name}") 
    print(f"There are {nboth} shared in {name}")
    
    venn2(subsets = (group1, group2, nboth), set_labels = (group1name, group2name), 
          set_colors=('r', 'g'), 
          alpha = 0.5);
    venn2_circles(subsets = (group1, group2, nboth))
    plt.title(f"{name}")
    plt.savefig(f"results/venn/{name}_{stat}.venn.png")
    plt.close()
    
def intersect2(one, two, df, write=True, path=None):
    inter = [x for x in list(one.GeneID) if x in list(two.GeneID)]
    length = len(inter)
    intersected_df = df[df.GeneID.isin(inter)]
    intersected_df.to_csv(f"{path}", sep="\t")
   
    return(length, intersected_df)


def add_columns_xlsx(name, de, fst, highfst, diffsnps, diffsnpsdf=None):
    
    rnaxlsx = pd.read_excel("results/genediff/RNA-Seq_diff.xlsx", 
                       sheet_name=name)
    
    highfst_bool = de.GeneID.isin(highfst.GeneID).astype(str)
    
    rnaxlsx['HighFst'] = highfst_bool
    
    if diffsnps:
        diffsnps_bool = de.GeneID.isin(diffsnpsdf.GeneID).astype(str)
        rnaxlsx['DiffSNPs'] = diffsnps_bool
    
    # add column of number of SNPs
    merged = pd.merge(de, fst, how="outer")
    rnaxlsx['nSNPs'] = merged['nSNPs']
    
    return(rnaxlsx)


#### Main ####
# Read contrasts in and other snakemake params
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons = [list(row) for i,row in comparisons.iterrows()]

percentile = snakemake.params['percentile']
diffsnps = snakemake.params['diffsnps']

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('results/RNA-Seq-full.xlsx', engine='xlsxwriter')

#### Differential expression v Fst venn diagram
for comp1,comp2 in comparisons:
    name = comp1 + "_" + comp2
    print(f"\n-------------- Venn Diagram for {name} --------------")
    de = pd.read_csv(f"results/genediff/{name}.csv")
    fst = pd.read_csv("results/variantAnalysis/selection/FstPerGene.tsv", sep="\t")
    #compare sig DE genes and top 5% fst genes?
    #get sig up and down diffexp genes
    sigde = de[de['padj'] < pval_threshold]
    sigde_up = sigde[sigde['FC'] > upper_fc]
    sigde_down = sigde[sigde['FC'] < lower_fc]

    #take top percentile of fst genes
    highfst = fst.nlargest(int(fst.shape[0]*percentile),f"{name}_zFst")

    #how many fst? how many sig de up and down?
    nfst = highfst.shape[0]
    nde_up = sigde_up.shape[0]
    nde_down = sigde_down.shape[0]

    print(f"There are {nde_up} significantly upregulated genes in {name}") 
    print(f"There are {nde_down} significantly downregulated genes in {name}")

    nboth, _ = intersect2(sigde_up, 
               highfst, 
               de, 
               write=True, 
               path=f"results/venn/{name}.DE.Fst.intersection.tsv")
    

    ###### XLSX file ######
    if diffsnps:
        diffsnpsDE = pd.read_csv("results/diffsnps/{name}.sig.kissDE.tsv", sep="\t")
        sheet = add_columns_xlsx(name, de, fst, highfst, diffsnps, diffsnpsDE)
    else:
        sheet = add_columns_xlsx(name, de, fst, highfst, diffsnps, diffsnpsDE=None)

    # Write each dataframe to a different worksheet.
    sheet.to_excel(writer, sheet_name=name)

# Close the Pandas Excel writer and output the Excel file.
writer.save()