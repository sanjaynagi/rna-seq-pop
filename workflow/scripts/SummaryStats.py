#!/usr/bin/env python3

"""
A script to calculate various SNP statistics, windowed population genetic statistics (Fst, PBS) and PCA.
Currently not modularised further to reduce the repetition of loading and filtering VCFs (which is slow). 
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *

# Read in parameters from snakemake
dataset = snakemake.params['dataset']
metadata = pd.read_csv(snakemake.input['metadata'], sep="\t")
metadata = metadata.sort_values(by='species')
chroms = snakemake.params['chroms']
ploidy = snakemake.params['ploidy']
numbers = get_numbers_dict(ploidy)
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']

# define functions
def isnotmissing(gn):
    """
    Function to check if missing values are present at a SNP
    """
    return((gn != -1).all())

#initialise dicts
total_snps_per_chrom = {}
snps_per_gene_allchroms = {}
snpeffdict = {}

seqdivdictchrom = {}
thetadictchrom = {}
coefdictchrom= {}

for i, chrom in enumerate(chroms):
    
    # Read in and Filter VCF
    path = f"results/variantAnalysis/vcfs/{dataset}.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=metadata,
                                                           numbers=numbers,
                                                           ploidy=ploidy,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)


    #### Genome-wide statistics (seqDiv, Wattersons Theta, LD, inbreeding coefficient) ####
    seqdivdict = {}
    thetadict = {}
    coefdict= {}
    allcoef = defaultdict(list)

    for pop in metadata['treatment'].unique():

        # Sequence diversity 
        seqdivdict[pop] = allel.sequence_diversity(pos, acsubpops[pop])
        
        # Wattersons theta
        thetadict[pop] = allel.watterson_theta(pos, acsubpops[pop])

        # Inbreeding coefficient
        if ploidy > 1:
            gn = geno.take(subpops[pop], axis=1)
            coef = allel.moving_statistic(gn, statistic=allel.inbreeding_coefficient, 
                                                size=1000, step=100)
            coef = np.nanmean(coef, axis=1)
            coefdict[pop] = np.mean(coef)
            allcoef[pop].append(np.array(coef))

        print(f"{pop},{chrom}, Sequence Diversity = ", seqdivdict[pop])
        print(f"{pop},{chrom}, Wattersons Theta = ", thetadict[pop])
        if ploidy > 1: print(f"{pop},{chrom}, Inbreeding Coef = ", np.mean(coef), "\n")

    seqdivdictchrom[chrom] = dict(seqdivdict)
    thetadictchrom[chrom] = dict(thetadict)
    if ploidy > 1: coefdictchrom[chrom] = dict(coefdict)

seqdivdictchrom= flip_dict(seqdivdictchrom)
thetadictchrom = flip_dict(thetadictchrom)
if ploidy > 1: coefdictchrom = flip_dict(coefdictchrom)

# Get stats per chromosome and plot heatmap
pidf = pd.DataFrame.from_dict(seqdivdictchrom)
pidf.to_csv("results/variantAnalysis/diversity/SequenceDiversity.tsv", sep="\t", index=True)
plotRectangular(pidf, path="results/variantAnalysis/diversity/piPerChrom.png", ylab="Chromosome", xlab="Treatment", figsize=[5,5], title=r'$\pi$')
thetadf = pd.DataFrame.from_dict(thetadictchrom)
plotRectangular(thetadf, path="results/variantAnalysis/diversity/thetaPerChrom.png", ylab="Chromosome", xlab="Treatment", figsize=[5,5], title=r'$\theta$')
thetadf.to_csv("results/variantAnalysis/diversity/WattersonsTheta.tsv", sep="\t", index=True)

thetamean = thetadf.apply(np.mean, axis=0)
pimean = pidf.apply(np.mean, axis=0)
summaryStats = pd.DataFrame({r'$\theta$':thetamean, r'$\pi$':pimean})
plotRectangular(summaryStats, path="results/variantAnalysis/diversity/pi_theta.overall.png", ylab="Treatment", xlab="Statistic", figsize=[5,5], rotate=False)

if ploidy > 1: pd.DataFrame.from_dict(coefdictchrom).to_csv("results/variantAnalysis/diversity/inbreedingCoef.tsv", sep="\t", index=True)

# Get genome wide average stats
if ploidy > 1:
    for pop in allcoef.keys():
        allcoef[pop] = np.nanmean(allcoef[pop])

    coefdf = pd.DataFrame.from_dict(allcoef, orient='index', columns=['InbreedingCoefficient'])
    coefdf.to_csv(f"results/variantAnalysis/diversity/inbreedingCoef.mean.tsv", sep="\t", index=True)