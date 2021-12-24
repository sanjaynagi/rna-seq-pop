#!/usr/bin/env python3

"""
A script to calculate windowed population genetic statistics (Fst, PBS).
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *

metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
chroms = snakemake.params['chroms']
ploidy = snakemake.params['ploidy']
numbers = get_numbers_dict(ploidy)
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']

#Fst/PBS window size
windowsizes = snakemake.params['windowsizes']
windowsteps = snakemake.params['windowsteps']
windownames = snakemake.params['windownames']

# Read in list of contrasts
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']
comparisons = [list(row) for i,row in comparisons.iterrows()]

for i, chrom in enumerate(chroms):

    path = f"results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=metadata,
                                                           numbers=numbers,
                                                           ploidy=ploidy,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)

    #### Fst in windows #### 
    for sus, res in comparisons:
        name = sus + "_" + res
        print(f"Calculating Fst values in sliding windows for {name}\n")

        for wname, size, step in zip(windownames, windowsizes, windowsteps):
            FstArray = allel.moving_hudson_fst(acsubpops[sus], 
                            acsubpops[res], 
                            size=size, step=step)
            midpoint = allel.moving_statistic(pos, np.mean, size=size, step=step)


            saveAndPlot("Fst", 
                        FstArray, 
                        midpoint, 
                        name=name,
                        prefix="results/variantAnalysis/selection/Fst/{name}.{chrom}.fst.{wname}.png", 
                        species=species, 
                        chrom=chrom, 
                        ylim=0.5, 
                        save=True)

            
        
    #### Population Branch Statistic (PBS) in windows ####
    if pbs:
        for pbscomp in pbscomps:
            pop1, pop2, outpop = pbscomp.split("_")

            print(f"Calculating PBS values in sliding window for {pbscomp}\n")
        
            for wname, size, step in zip(windownames, windowsizes, windowsteps):
                pbsArray = allel.pbs(acsubpops[pop1], 
                                acsubpops[pop2], 
                                acsubpops[outpop], 
                                window_size=size, window_step=step, normed=True)
                midpoint = allel.moving_statistic(pos, np.mean, size=size, step=step)

                plt.figure(figsize=[20,8])
                sns.lineplot(midpoint, pbsArray)
                plt.title(f"PBS {chrom} {pbscomp}")
                plt.savefig(f"results/variantAnalysis/selection/pbs/{pbscomp}.{chrom}.pbs.{wname}.png")
                plt.close()