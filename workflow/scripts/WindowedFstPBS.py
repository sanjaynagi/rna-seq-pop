#!/usr/bin/env python3

"""
A script to calculate windowed population genetic statistics (Fst, PBS).
"""

import sys
sys.stderr = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np
import allel
import rnaseqpoptools as rnaseqpop
from pathlib import Path

dataset = snakemake.params['dataset']
metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
contigs = snakemake.params['contigs']
ploidy = snakemake.params['ploidy']
numbers = rnaseqpop.get_numbers_dict(ploidy)
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']

#Fst/PBS window size
windownames =  ['1000snp_window', '2000snp_window', '5000snp_window']
windowsizes =  [1000, 2000, 5000]
windowsteps = [500, 1000, 1000]

# Read in list of contrasts
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']
comparisons = [list(row) for i,row in comparisons.iterrows()]

for i, contig in enumerate(contigs):

    path = f"results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = rnaseqpop.readAndFilterVcf(path=path,
                                                           contig=contig,
                                                           samples=metadata,
                                                           numbers=numbers,
                                                           ploidy=ploidy,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)

    #### Fst in windows #### 
    for sus, res in comparisons:
        name = sus + "_" + res
        cohortText = f"{sus} v {res}"
        print(f"Calculating Fst values in sliding windows for {name}\n")

        for wname, size, step in zip(windownames, windowsizes, windowsteps):
            
            if geno.shape[0] < size:
                print(f"Skipping {wname} for {name} because there are not enough SNPs in {contig}.")
                print("Touching file to prevent snakemake from erroring out.")
                Path(f"results/variantAnalysis/selection/fst/{wname}/{name}.Fst.{contig}.svg").touch()
            else:
                FstArray = allel.moving_hudson_fst(acsubpops[sus], 
                                acsubpops[res], 
                                size=size, step=step)
                midpoint = allel.moving_statistic(pos, np.median, size=size, step=step)
                
                cohortNoSpaceText = wname + "/" + name 
                rnaseqpop.plotWindowed(statName="Fst",
                            cohortText=cohortText,
                            cohortNoSpaceText=cohortNoSpaceText,
                            values=FstArray, 
                            midpoints=midpoint,
                            colour='dodgerblue',
                            prefix="results/variantAnalysis/selection/fst", 
                            contig=contig, 
                            ylim=1, 
                            save=True)

        
    #### Population Branch Statistic (PBS) in windows ####
    if pbs:
        for pbscomp in pbscomps:
            pop1, pop2, outpop = pbscomp.split("_")
            cohortText = f"(({pop1}, {pop2}), {outpop})"
            print(f"Calculating PBS values in sliding window for {pbscomp}\n")
        
            for wname, size, step in zip(windownames, windowsizes, windowsteps):

                if geno.shape[0] < size:
                    print(f"Skipping {wname} for {pbscomp} because there are not enough SNPs in {contig}.")
                    print("Touching file to prevent snakemake from erroring out.")
                    Path(f"results/variantAnalysis/selection/pbs/{wname}/{pbscomp}.PBS.{contig}.svg").touch()
                else:
                    pbsArray = allel.pbs(acsubpops[pop1], 
                                    acsubpops[pop2], 
                                    acsubpops[outpop], 
                                    window_size=size, window_step=step, normed=True)
                    midpoint = allel.moving_statistic(pos, np.median, size=size, step=step)

                    cohortNoSpaceText =  wname + "/" + pbscomp
                    rnaseqpop.plotWindowed(statName="PBS", 
                                cohortText=cohortText,
                                cohortNoSpaceText=cohortNoSpaceText,
                                values=pbsArray, 
                                midpoints=midpoint, 
                                colour='dodgerblue',
                                prefix="results/variantAnalysis/selection/pbs",
                                contig=contig, 
                                ylim=0.5, 
                                save=True)