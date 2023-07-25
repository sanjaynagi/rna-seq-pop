#!/usr/bin/env python3

"""
A script to plot the frequencies of variants of interest as a heatmap
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import matplotlib
import rnaseqpoptools as rnaseqpop


### Variants of Interest patjh ###
voiPath = snakemake.input['VariantsOfInterest']
## Read VOI data
muts = pd.read_csv(voiPath, sep="\t")

## separate contig and pos data and sort 
muts['chrom'] = muts['Location'].str.split(":").str.get(0)
muts['pos'] = muts['Location'].str.split(":").str.get(1).str.split("-").str.get(0)
muts = muts.sort_values(['chrom', 'pos'])


## Run for all samples
df, annot = rnaseqpop.getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/{mut}_alleleBalance.csv", var="sample")
rnaseqpop.plotRectangular(df, annot=annot, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerSample.svg")
rnaseqpop.plotRectangular(df, annot=annot, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerSample.pdf")


## Run for avarage frequencies across treatments
df2, annot2 = rnaseqpop.getAlleleFreqTable(muts, "results/variantAnalysis/variantsOfInterest/csvs/mean_{mut}_alleleBalance.csv", var="treatment", mean_=True)
rnaseqpop.plotRectangular(df2, annot=annot2, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.svg", xlab="strain")
rnaseqpop.plotRectangular(df2, annot=annot2, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapPerTreatment.pdf", xlab="strain")

# Join both plots
rnaseqpop.plotTwoRectangular(df, annot, df2, annot2, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapBothPlots.svg", ratio='auto')
rnaseqpop.plotTwoRectangular(df, annot, df2, annot2, path="results/variantAnalysis/variantsOfInterest/VOI.heatmapBothPlots.pdf", ratio='auto')