#!/usr/bin/env python3

"""
A script to plot the frequencies of variants of interest as a heatmap
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *

### AIMS ###
voiData = snakemake.input['VariantsOfInterest']

## Read VOI data
muts = pd.read_csv(voiData, sep="\t")

## separate chrom and pos data and sort 
muts['chrom'] = muts['Location'].str.split(":").str.get(0)
muts['pos'] = muts['Location'].str.split(":").str.get(1).str.split("-").str.get(0)
muts = muts.sort_values(['chrom', 'pos'])


def getAlleleFreqTable(muts, Path, var="sample", mean_=False):
    """
    This function loops through the mutations, reads in its csv file containing frequencies
    and stores as a rectangular matrix for heatmap plotting.
    """
    Dict = {}
    for mut in muts['Name']:
        df = pd.read_csv(Path.format(mut=mut))
        if mean_:
            df['gene'] = muts[muts.Name == mut]['Gene'].iloc[0]
        df['name'] = df['chrom'] + ":"+ df['pos'].astype(str) + "  " + df['gene'] + " | " + df['mutation']
        df['frequency'] = df.filter(like="proportion").sum(axis=1)
        Dict[mut] = df[['name', var, 'frequency']]

    voiData = pd.concat(Dict)
    # Make rectangular
    voiFreqTable = voiData.pivot(index="name", columns=var).round(2).droplevel(0, axis=1) 
    return(voiFreqTable)



## Run for all samples
df = getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/{mut}_alleleBalance.csv", var="sample")
plotRectangular(df, path="results/variantAnalysis/variantsOfInterest/csvs/{mut}_allele_balance.csv")


## Run for avarage frequencies across treatments
df = getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/mean_{mut}_allele_balance.csv", var="treatment", mean_=True)
plotRectangular(df, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.png", xlab="strain")

    