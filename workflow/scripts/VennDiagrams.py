#!/usr/bin/env python3

"""
A script to get the intersections of Differential expression results, Fst, PBS and differential SNPs analysis.
Draws Venn diagrams and adds columns to RNA-seq-diff.xlsx, whether the gene has high Fst/PBS/diffsnps. 
"""
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


def add_columns_xlsx(name, fst, highfst, diffsnps):
    
    rnaxlsx = pd.read_excel("results/genediff/RNA-Seq_diff.xlsx", 
                       sheet_name=name)
    
    highfst_bool = de.GeneID.isin(highfst.GeneID).astype(str)
    diffsnps_bool = de.GeneID.isin(diffsnps.GeneID).astype(str)
    
    rnaxlsx['HighFst'] = highfst_bool
    rnaxlsx['DiffSNPs'] = diffsnps_bool
    
    # add column of number of SNPs
    merged = pd.merge(de, fst, how="outer")
    rnaxlsx['nSNPs'] = merged['nSNPs']
    
    return(rnaxlsx)


#### Main ####
# Read contrasts in and other snakemake params
comparisons = pd.read_csv(snakemake.input['DEcontrasts'])
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons = [list(row) for i,row in comparisons.iterrows()]

pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
percentile = snakemake.params['percentile']


# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('results/RNA-Seq-full.xlsx', engine='xlsxwriter')

#### Differential expression v Fst venn diagram
for comp1,comp2 in comparisons:
    name = comp1 + "_" + comp2
    print(f"\n-------------- Venn Diagram for {name} --------------")
    de = pd.read_csv(f"results/genediff/{name}.csv")
    fst = pd.read_csv("results/variants/fst.tsv", sep="\t")
    diffsnps = pd.read_csv(f"results/variants/diffsnps/{name}.sig.kissDE.tsv", sep="\t")
    #compare sig DE genes and top 5% fst genes?
    #get sig up and down diffexp genes
    sigde = de[de['padj'] < 0.05]
    sigde_up = sigde[sigde['FC'] > 1]
    sigde_down = sigde[sigde['FC'] < 1]

    #take top 5% of fst genes
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
    
    plotvenn2(name, nde_up, nfst, nboth,stat="DE.Fst", group1name='Significant up DE genes', group2name='High Fst')
    
    ###### XLSX file ######
    sheet = add_columns_xlsx(name, fst, highfst, diffsnps)
    # Write each dataframe to a different worksheet.
    sheet.to_excel(writer, sheet_name=name)

# Close the Pandas Excel writer and output the Excel file.
writer.save()


if pbs is True:
    for pbscomp in pbscomps:
        pop1, pop2, outpop = pbscomp.split("_")
        name2 = pop2 + "_" + pop1

        DEfile = Path(f"results/genediff/{name2}.csv")
        if DEfile.is_file() is False:
            continue
            
        de = pd.read_csv(DEfile)
        fst = pd.read_csv("results/variants/fst.tsv", sep="\t")
        pbs = pd.read_csv("results/variants/pbs.tsv", sep="\t")

        # compare sig DE genes and top 5% fst/pbs genes?
        # get sig up and down diffexp genes
        sigde = de[de['padj'] < 0.05]
        sigde_up = sigde[sigde['FC'] > 1]
        sigde_down = sigde[sigde['FC'] < 1]

        # take top 5% of fst genes
        highpbs = pbs.nlargest(int(pbs.shape[0]*percentile),f"{pbscomp}PBS")
        highfst = fst.nlargest(int(fst.shape[0]*percentile),f"{name2}_zFst")

        # how many fst? how many sig de up and down?
        nfst = highfst.shape[0]
        npbs = highpbs.shape[0]
        nde_up = sigde_up.shape[0]
        nde_down = sigde_down.shape[0]

        print(f"There are {nde_up} significantly upregulated genes in {name2}") 
        print(f"There are {nde_down} significantly downregulated genes in {name2}")

        # intersect of fst and pbs, write out
        nbothdefst, _ = intersect2(sigde_up, highfst, de, write=False)
        nbothdepbs, _ = intersect2(sigde_up, highpbs, de, write=False)
        nbothfstpbs, _ = intersect2(highfst, highpbs, de, write=True, 
                                    path=f"results/venn/{name2}.Fst.PBS.intersection.tsv")
        
        # plot venns, DEvPBS and FSTvPBS
        plotvenn2(pbscomp, nde_up, npbs, nbothdepbs,stat="DEup.PBS", group1name='Significant up DE genes', group2name='High PBS')
        plotvenn2(pbscomp, nfst, npbs, nbothfstpbs,stat="Fst.PBS", group1name='High Fst', group2name='High PBS')

        # get intersection of all three (FST PBS and DE) results 
        threesets = list(set(sigde_up.GeneID).intersection(highfst.GeneID, highpbs.GeneID))
        de_intersected3 = de[de.GeneID.isin(threesets)]
        nall = de_intersected3.shape[0]
        de_intersected3.to_csv(f"results/venn/{name2}_DE.Fst.PBS.intersection.tsv", sep="\t")

        # three way venn diagram
        venn3(subsets = (nde_up, nfst,nbothdefst , npbs, nbothdepbs , nbothfstpbs, nall), 
              set_labels = ('Significant DE', 'High Fst', 'High PBS'), 
              alpha = 0.5)
        venn3_circles(subsets = (nde_up, nfst,nbothdefst , npbs, nbothdepbs , nbothfstpbs, nall))
        plt.title(f"{pbscomp}")
        plt.savefig(f"results/venn/{name}_DE_Fst_PBS.venn.png")
        plt.close()
