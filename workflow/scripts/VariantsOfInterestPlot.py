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


## Run for all samples
df, annot = getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/{mut}_alleleBalance.csv", var="sample")
plotRectangular(df, path="results/variantAnalysis/variantsOfInterest/csvs/{mut}_allele_balance.csv")


## Run for avarage frequencies across treatments
df2, annot2 = getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/mean_{mut}_allele_balance.csv", var="treatment", mean_=True)
plotRectangular(df, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.png", xlab="strain")


plotTwoRectangular(df, annot, df2, annot2, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapBothPlots.png", ratio='atuo')