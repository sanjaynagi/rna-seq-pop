#!/usr/bin/env python3

"""
A script that checks if DE genes lie underneath known selective sweeeps in the Ag
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *


pval_threshold = snakemake.params['pval']
fc_threshold = snakemake.params['fc']
# Read in list of contrasts
comparisons = pd.read_csv(snakemake.input['DEcontrasts'])

signals = pd.read_csv("resources/signals.csv")

for comp in comparisons['contrast']:
    
    DEgenes = pd.read_csv(f"results/genediff/{comp}.csv")
    sigup = DEgenes[np.logical_and(DEgenes['padj'] < pval_threshold, DEgenes['FC'] > fc_threshold)]
    
    sweep = {}
    nswept = {}

    for i, cols in signals.iterrows():

        if pd.isnull(cols['overlapping_genes']):
            continue

        sweptgenes = np.array(cols['overlapping_genes'].split(" "))


        overlap = np.isin(sweptgenes, sigup['GeneID'])

        sweep[cols['uid']] = sweptgenes[overlap]
        nswept[cols['uid']] = sweptgenes
        
    genes = np.concatenate(list(sweep.values()))
    swept = sigup[np.isin(sigup['GeneID'], genes)]
    
    for k,v in sweep.items():
        sweep[k] = ' '.join(v)

    sweptDE = pd.DataFrame.from_dict(sweep, orient='index', columns=['overlapping_DE_genes'])
    sweptDE = sweptDE.reset_index().rename(columns={'index': 'Ag1000g_sweep'})
    sweptDE['overlapping_genes'] = signals['overlapping_genes'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['chromosome'] = signals['peak_end_seqid'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['epicenter'] = signals['epicenter_coord'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)
    sweptDE['known_loci'] = signals['known_loci'][~pd.isnull(signals['overlapping_genes'])].reset_index(drop=True)

    wheresweep = defaultdict(dict)
    whatsweep = defaultdict(list)

    for gene in genes:

        for i, cols in sweptDE.iterrows():

            sweptgenes = np.array(cols['overlapping_DE_genes'].split(" "))

            if np.isin(sweptgenes, gene).any():
                wheresweep[gene]['chrom'] = cols['chromosome']
                wheresweep[gene]['epicenter'] = cols['epicenter']
                wheresweep[gene]['known_loci'] = cols['known_loci']

                whatsweep[gene].append(cols['Ag1000g_sweep'])

    for k,v in whatsweep.items():
        whatsweep[k] = ' '.join(v)
        
    dfwhere = pd.DataFrame.from_dict(wheresweep, orient='index')
    dfwhat = pd.DataFrame.from_dict(whatsweep, orient='index', columns=['Ag1000g_sweeps'])

    df = pd.concat([dfwhat, dfwhere], axis=1)
    df = df.reset_index().rename(columns={'index': 'GeneID'})

    swept = swept.merge(df)
    swept.to_csv(f"results/genediff/ag1000gSweeps/{comp}_swept.tsv", sep="\t")