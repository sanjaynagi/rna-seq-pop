#!/usr/bin/env python3

"""
A script that checks if DE genes lie underneath known selective sweeeps in the Ag1000g
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import numpy as np
from collections import defaultdict

# Read in parameters 
pval_threshold = snakemake.params['pval']
upper_fc = snakemake.params['fc']
lower_fc = 1/upper_fc # if someone wants a FC threshold of 2, need to have lower threshold of 0.5.
# Read in list of contrasts
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])

# Read in .csv file containing selection signals and associated metadata
signals = pd.read_csv("resources/signals.csv")

for comp in comparisons['contrast']:
    
    # Load DE results
    DEgenes = pd.read_csv(f"results/genediff/{comp}.csv")
    # Filter to overexpressed genes
    sigup = DEgenes[np.logical_and(DEgenes['padj'] < pval_threshold, np.logical_or(DEgenes['FC'] > upper_fc, DEgenes['FC'] < lower_fc))]
    
    sweep = {}
    nswept = {}

    # Loop through each signal, recording if any DE genes are found in one
    for i, cols in signals.iterrows():

        if pd.isnull(cols['overlapping_genes']): # Skip signal if no overlapping genes
            continue

        sweptgenes = np.array(cols['overlapping_genes'].split(" "))

        # Get boolean array - if list of swept genes isin our DE genes
        overlap = np.isin(sweptgenes, sigup['GeneID'])

        sweep[cols['uid']] = sweptgenes[overlap]
        nswept[cols['uid']] = sweptgenes
        
    genes = np.concatenate(list(sweep.values()))
    swept = sigup[np.isin(sigup['GeneID'], genes)]
    
    for k,v in sweep.items():
        sweep[k] = ' '.join(v)

    # Build dataframe and add sweep metadata columns
    sweptDE = pd.DataFrame.from_dict(sweep, orient='index', columns=['overlapping_DE_genes'])
    sweptDE = sweptDE.reset_index().rename(columns={'index': 'Ag1000g_sweep'})
    sweptDE['overlapping_genes'] = signals['overlapping_genes'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['chromosome'] = signals['peak_end_seqid'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['epicenter'] = signals['epicenter_coord'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['known_loci'] = signals['known_loci'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)

    wheresweep = defaultdict(dict)
    whatsweep = defaultdict(list)

    # Now loop through each swept gene to find name of sweeps it lies under
    # And location of sweep
    for gene in genes:

        for i, cols in sweptDE.iterrows():

            sweptgenes = np.array(cols['overlapping_DE_genes'].split(" "))

            if np.isin(sweptgenes, gene).any():
                wheresweep[gene]['contig'] = cols['chromosome']
                wheresweep[gene]['epicenter'] = cols['epicenter']
                wheresweep[gene]['known_loci'] = cols['known_loci']

                whatsweep[gene].append(cols['Ag1000g_sweep'])

    # Join name of sweeps column into a single string, so it fits in one column of data frame
    for k,v in whatsweep.items():
        whatsweep[k] = ' '.join(v)
        
    dfwhere = pd.DataFrame.from_dict(wheresweep, orient='index')
    dfwhat = pd.DataFrame.from_dict(whatsweep, orient='index', columns=['Ag1000g_sweeps'])

    df = pd.concat([dfwhat, dfwhere], axis=1)
    df = df.reset_index().rename(columns={'index': 'GeneID'})

    swept = swept.merge(df)
    swept.to_csv(f"results/genediff/ag1000gSweeps/{comp}_swept.tsv", sep="\t")