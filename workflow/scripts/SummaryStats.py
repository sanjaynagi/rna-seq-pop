#!/usr/bin/env python3

"""
A script to calculate various SNP statistics, windowed population genetic statistics (Fst, PBS) and PCA.
Currently not modularised further to reduce the repetition of loading and filtering VCFs (which is slow). 
"""

import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
import allel
from collections import defaultdict
import rnaseqpoptools as rnaseqpop

# Read in parameters from snakemake
dataset = snakemake.params['dataset']
metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
contigs = snakemake.params['contigs']
ploidy = snakemake.params['ploidy']
numbers = rnaseqpop.get_numbers_dict(ploidy)
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']


# Initialise dicts to store genetic diversity statistic
pi = {}
theta = {}
coefdictchrom= {}

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


    # Genome-wide statistics (Pi, Wattersons Theta, inbreeding coefficient)
    pi[contig] = rnaseqpop.windowedDiversity(geno=geno, pos=pos, subpops=subpops, statistic='pi', window_size=20_000)
    theta[contig] = rnaseqpop.windowedDiversity(geno=geno, pos=pos, subpops=subpops, statistic='theta', window_size=20_000)    
    
    coefdict= {}
    allcoef = defaultdict(list)
    for pop in metadata['treatment'].unique():
        # Inbreeding coefficient
        if ploidy > 1:
            gn = geno.take(subpops[pop], axis=1)
            coef = allel.moving_statistic(gn, statistic=allel.inbreeding_coefficient, size=1000, step=100)
            coef = np.nanmean(coef, axis=1)
            coefdict[pop] = np.mean(coef)
            allcoef[pop].append(np.array(coef))

        if ploidy > 1: print(f"{pop} | {contig} | Inbreeding Coef =", np.mean(coef), "\n")
    if ploidy > 1: coefdictchrom[contig] = dict(coefdict)

# Concat contigs, get CIs for Pi and Theta and save to file
pi_df = rnaseqpop.diversity_ci_table(div_dict=pi, statistic='pi').to_csv("results/variantAnalysis/diversity/SequenceDiversity.tsv", sep="\t", index=True)
theta_df = rnaseqpop.diversity_ci_table(div_dict=theta, statistic='theta').to_csv("results/variantAnalysis/diversity/WattersonsTheta.tsv", sep="\t", index=True)

if ploidy > 1: coefdictchrom = rnaseqpop.flip_dict(coefdictchrom)
if ploidy > 1: pd.DataFrame.from_dict(coefdictchrom).to_csv("results/variantAnalysis/diversity/inbreedingCoef.tsv", sep="\t", index=True)
# Get genome wide average stats
if ploidy > 1:
    for pop in allcoef.keys():
        allcoef[pop] = np.nanmean(allcoef[pop])

    coefdf = pd.DataFrame.from_dict(allcoef, orient='index', columns=['InbreedingCoefficient'])
    coefdf.to_csv(f"results/variantAnalysis/diversity/inbreedingCoef.mean.tsv", sep="\t", index=True)