#!/usr/bin/env python3

"""
A script to perform PCA on the genotype data
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import rnaseqpoptools as rnaseqpop
import pandas as pd 
import numpy as np
import allel
from adjustText import adjust_text


# Read in parameters from snakemake
dataset = snakemake.params['dataset']
metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
contigs = snakemake.params['contigs']
ploidy = snakemake.params['ploidy']
numbers = rnaseqpop.get_numbers_dict(ploidy)
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']


for i, contig in enumerate(contigs):
    
    # Read in and Filter VCF
    path = f"results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = rnaseqpop.readAndFilterVcf(path=path,
                                                           contig=contig,
                                                           samples=metadata,
                                                           numbers=numbers,
                                                           ploidy=ploidy,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)
    

    #### Principal Components Analysis (PCA) ####
    # Set up dict to store indices for colours
    d={}
    for name, inds in subpops.items():
        for n in range(len(inds)):
            p = inds[n]
            d[p] = name

    # Store dict as a dataframe and get colours 
    treatment_indices = pd.DataFrame.from_dict(d, orient='index').reset_index()
    treatment_indices = treatment_indices.rename(columns = {'index':'sample_index', 0:"name"})
    pop_colours = rnaseqpop.get_colour_dict(treatment_indices['name'], "viridis")
    
    # Run PCA function defined in tools.py
    print(f"Performing PCA on {dataset} chromosome {contig}")
    rnaseqpop.pca(geno, contig, ploidy, dataset, populations, metadata, pop_colours, prune=True, scaler=None)