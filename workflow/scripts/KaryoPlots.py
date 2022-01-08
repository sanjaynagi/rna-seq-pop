#!/usr/bin/env python3

"""
Plots Karyotype data
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *

# Read in parameters from snakemake
ploidy = snakemake.params['ploidy']
invs = snakemake.params['inversions']
metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')

##
karyo = {}
for inv in invs:
    df = pd.read_csv(f"../../results/karyotype/{inv}.karyo.txt", sep="\s+", header=None)
    df = df.rename(columns={0:'Sample', 1:'KaryoScore', 2:'n_SNPtags'})
    df[inv] = df['KaryoScore']/ploidy
    df.rename({inv:f'{inv} frequency'}).to_csv(f"../../results/karyotype/{inv}.karyo.txt", sep="\t")
    
    karyo[inv] = df[['Sample', inv]]

# concat all the dfs in the dict and remove duplicate cols
karyo = pd.concat(karyo.values(), axis=1).T.drop_duplicates().T.set_index("Sample")

## transpose and round to 2 decimals
karyo = karyo.T.astype("float16").round(2)

plotRectangular(karyo, path="results/karyotype/karyoFreqs.png" , cmap='mako_r', ylab='Inversion', figsize=[10,5])

# Produce for average karyos per treatment
df = karyo.T.reset_index()
df['Sample'] = df['Sample'].str.strip('12345678910')
df = df.groupby("Sample").agg('mean').T
plotRectangular(df, path="../../results/karyotype/karyoOverallFreqs.png", ylab='Inversion', cbar=False, figsize=[2,len(invs)])