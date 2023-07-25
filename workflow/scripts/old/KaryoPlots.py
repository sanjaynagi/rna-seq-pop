#!/usr/bin/env python3

"""
Plots Karyotype data
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import rnaseqpoptools as rnaseqpop
import pandas as pd 
import matplotlib.pyplot as plt

# Read in parameters from snakemake
ploidy = snakemake.params['ploidy']
invs = snakemake.params['inversions']
metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
dataset = snakemake.params['dataset']


karyo = {}
for inv in invs:
    df = pd.read_csv(f"results/karyotype/{inv}.{dataset}.karyo.txt", sep="\s+", header=None)
    df = df.rename(columns={0:'sampleID', 1:'KaryoScore', 2:'n_SNPtags'})
    df[inv] = df['KaryoScore']/ploidy
    df.rename(columns={inv:f'{inv} frequency'}).to_csv(f"results/karyotype/{inv}.{dataset}.karyo.txt", sep="\t")
    
    karyo[inv] = df[['sampleID', inv]]

# concat all the dfs in the dict and remove duplicate cols
karyo = pd.concat(karyo.values(), axis=1).T.drop_duplicates().T.set_index("sampleID")

## transpose and round to 2 decimals
karyo = karyo.T.astype("float64").round(2)
rnaseqpop.plotRectangular(karyo, path="results/karyotype/karyoFreqs.svg" , cmap='mako_r', ylab='Inversion', figsize=[10,5])

# Produce for average karyos per treatment
df = karyo.T.reset_index()
df = df.merge(metadata[['sampleID', 'treatment']])
df = df.groupby("treatment").agg('mean').T.astype("float64").round(2)
rnaseqpop.plotRectangular(df, path="results/karyotype/karyoOverallFreqs.svg", ylab='Inversion', cbar=False, figsize=[8,4])