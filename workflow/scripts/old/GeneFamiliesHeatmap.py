#!/usr/bin/env python3

"""
Construct heatmap of different gene families
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import rnaseqpoptools as rnaseqpop
import pandas as pd
import numpy as np
import scipy
from matplotlib.collections import LineCollection
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])
metadata = pd.read_csv("config/samples.tsv", sep="\t")

# Read in .csv file containing selection signals and associated metadata
pfam_df = pd.read_csv(snakemake.input['pfam'], sep="\s+", header=None)
go_df = pd.read_csv(snakemake.input['eggnog'], sep="\t", header=None)

pfam_df.columns = ["transcript", "pstart", "pend", "pfamid", "domain", "domseq"]
go_df.columns = ['transcript', 'GO_terms']

gene_annot_df = pfam_df.merge(go_df)
gene_annot_df.loc[:, 'gene_id'] = gene_annot_df.loc[:, 'transcript'].str.replace("Anogam_", "").str.replace("-R[A-Z]", "")

# load normalised counts 
counts = pd.read_csv(snakemake.input['normcounts'], sep="\t")


## define functions ## 

def plot_rotated_dendro(dendro, ax, linewidth=2):
    coords = zip(dendro['dcoord'], dendro['icoord'])
    lines = LineCollection([list(zip(x, y)) for x, y in coords], colors='lightgray', linewidth=linewidth)

    ax.add_collection(lines)
    number_of_leaves = len(order)
    max_dependent_coord = max(map(max, dendro['dcoord']))

    ax.yaxis.set_ticks_position('right')

    # Constants 10 and 1.05 come from
    # `scipy.cluster.hierarchy._plot_dendrogram`
    ax.set_ylim(0, number_of_leaves * 10)
    ax.set_xlim(0, max_dependent_coord * 1.05)

    ax.invert_xaxis()
    ax.invert_yaxis()

    sns.despine(ax=ax, bottom=True, left=True, right=True)
    ax.tick_params(top=False,
                   bottom=False,
                   left=False,
                   right=False,
                   labelleft=False,
                   labelright=False,
                   labelbottom=False)
    return(ax)

def gene_ids_from_domain(gene_annot_df, domain):
    gene_list = []
    if isinstance(domain, list):
        for dom in domain:
            ids = gene_annot_df.query("domain == @domain")['gene_id']
            gene_list.append(ids)
            return(np.unique(gene_list))
    else:
        return(gene_annot_df.query("domain == @domain")['gene_id'].to_numpy())





### MAIN

# a dict with gene families and their respective Pfam domain for extracting
gene_fams = {'COE': 'COesterase', 
             'OBP':'PBP_GOBP', 
             'OR':'7tm_6', 
             'UGT':'UDPGT',
             'Io':['Lig_chan','7tm_1'],
             'Gr': '7tm_7',
             'P450':'p450', 
             'Gst':['GST_N', 'GST_N_3', 'GST_C'], 
             'ABC':['ABC_membrane', 'ABC_tran'],
            'FAS':'ketoacyl-synt',
            'ELO':'ELO',
            'FAD':'FA_desaturase',
            'FAR':'NAD_binding_4'}


### read in diff exp data
genediff = {}
counts_dict = {}
for comp in comparisons['contrast']:
    df = pd.read_csv(f"results/genediff/{comp}.csv", sep=",")
    names = df[['GeneID', 'GeneName']]#, 'GeneDescription']]
    genediff[comp] = df
    
    #loop through gene fams 
    fam_dict = {}
    for fam, domain in gene_fams.items():
        gene_ids = gene_ids_from_domain(gene_annot_df, domain)
        fam_de = df.query("GeneID in @gene_ids").copy()
        
        fam_de.loc[:, 'de'] = np.where(
             np.logical_and(fam_de['padj'].between(0, 0.05, inclusive="both"), fam_de['FC'].ge(1)), 
            'Upregulated', 
             np.where(
                     np.logical_and(fam_de['padj'].between(0, 0.05, inclusive="both"), fam_de['FC'].le(1)), 'Downregulated', 'No Diff. expression'
             )
        )
        
        fam_dict[fam] = fam_de
        counts_dict[fam] = counts.query("GeneID in @gene_ids")
        counts_dict[fam] = counts_dict[fam].merge(names)
        counts_dict[fam].loc[:, 'Label'] = [id_ + " | " + name if name != "" else id_ for id_, name in zip(counts_dict[fam]['GeneID'].fillna(""), counts_dict[fam]['GeneName'].fillna(""))]
        
    genediff[comp] = dict(fam_dict)

colors = ["#d9d9d9", '#33cc33', '#cc66ff']
pal = sns.color_palette(colors)

def remap(x):
    if x == 'Upregulated':
        return(1)
    elif x == 'Downregulated':
        return(2)
    elif x == 'No Diff. expression':
        return(0)
    else:
        assert x == 'Nan', "whats wrong with x"
    

def check_fcs(diffexp, fam, comparisons):
    
    idx = np.random.choice(diffexp.shape[0])
    gene = diffexp['GeneID'][idx]
    
    vals = diffexp.set_index('GeneID').query("GeneID == @gene").to_numpy()[0]
    
    comp_vals = []
    for comp in comparisons:
        row = genediff[comp][fam].query("GeneID == @gene")
        if row['padj'].to_numpy() < 0.05:
            if row['FC'].to_numpy() > 1:
                comp_vals.append(1)
            elif row['FC'].to_numpy() < 1:
                comp_vals.append(2)
        else:
            comp_vals.append(0)
    
    assert all(np.array(comp_vals) == vals), "The values do not match! ffs"

with PdfPages("results/genediff/GeneFamiliesHeatmap.pdf") as pdf:

    for fam, domain in gene_fams.items():

        ids = counts_dict[fam]['GeneID'].to_numpy()
        diffexp = counts_dict[fam]['GeneID'].to_frame().copy()

        for comp in comparisons['contrast']:
            de_data = genediff[comp][fam].set_index("GeneID").reindex(ids).reset_index()
            assert all(diffexp['GeneID'].to_numpy() == de_data['GeneID'].to_numpy()), "wrong!!"
            diffexp.loc[:, comp] = de_data['de']

        diffexp_remap = diffexp.set_index("GeneID").reindex(index=ids)
        diffexp_remap = diffexp_remap.applymap(remap)
        check_fcs(diffexp_remap.reset_index(), fam, comparisons['contrast'])

        df = counts_dict[fam].set_index("Label").drop(columns=['GeneID', 'GeneName'])
        size = len(df)

        g = sns.clustermap(data=np.log2(df+1), 
                        cmap="Blues", 
                        vmax=10, 
                        cbar=False,
                        linewidths=4,
                        col_cluster=False,
                        linecolor="white",
                        yticklabels=False, 
                        xticklabels=False,
                        tree_kws=dict(linewidths=3, colors="lightgray"),
                        figsize=[1, 1])
        ax = g.ax_heatmap
        dendro = g.ax_row_dendrogram
        g.ax_cbar.set_visible(False)
        g.ax_row_dendrogram.set_visible(False) #suppress row dendrogram
        g.ax_col_dendrogram.set_visible(False) #suppress column dendrogram

        order = g.dendrogram_row.reordered_ind

        diffexp_remap = diffexp_remap.iloc[order, :]
        assert all(diffexp_remap.index == df.iloc[order, :].index.str[:10].to_list()), "rows do not match"

        df = df.T.assign(sampleID=df.T.index)
        df.loc[:, 'treatment'] = df['sampleID'].str.replace("[\d+]", "", regex=True)
        n_treatments = len(df.loc[:, 'treatment'].unique())

        fig, ax = plt.subplots(1, n_treatments+2, figsize=[10,(size/3)], gridspec_kw={'width_ratios':[3, 3] + list(np.repeat(1, n_treatments))})
        fig.suptitle(f"{fam}, pfam domains = {domain}", fontsize=18, fontweight='bold')
        dendro = scipy.cluster.hierarchy.dendrogram(g.dendrogram_row.linkage, no_plot=True,
                                                        above_threshold_color='lightgray', 
                                                        color_threshold=0, 
                                                        no_labels=True)
        ax[0] = plot_rotated_dendro(dendro, ax[0])    
        sns.heatmap(ax=ax[1], data=diffexp_remap, cmap=pal, cbar=False, linewidths=4, linecolor='white')
        ax[1].set_xticklabels(comparisons['contrast'], fontsize=14, fontweight='bold')
        ax[1].set_ylabel("")

        for idx, group in enumerate(metadata.treatment.unique()):
            idx += 2
            df2 = df.query("treatment == @group").drop(columns=['treatment', 'sampleID'])
            df2 = df2.T.iloc[order, :]

            sns.heatmap(ax=ax[idx], data=np.log10(df2+1), cmap="Blues", cbar=False, linewidths=4, linecolor='white')

            ax[idx].set_xlabel(group, rotation=90, fontsize=14, fontweight='bold', ha='right', rotation_mode='anchor')
            ax[idx].set_xticklabels([])
            ax[idx].set_ylabel("")
            ax[idx].tick_params(top=False,
                                bottom=False,
                                left=False,
                                right=False,
                                labelleft=False,
                                labelright=True,
                                labelbottom=True)

            sns.despine(left=True, bottom=True, right=True)
            ax[idx].set_yticklabels(labels=df2.index, rotation=0,fontsize=14)
        for axes in ax[1:-1]:
            axes.set_yticklabels([])
            axes.set_yticks([])

        fig.savefig(f"results/genediff/{fam}_expression.tiff", dpi=200)
        pdf.savefig(bbox_inches = 'tight')