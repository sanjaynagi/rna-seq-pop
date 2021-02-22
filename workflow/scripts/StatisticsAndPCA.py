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
samples = pd.read_csv(snakemake.input['samples'], sep="\t")
samples = samples.sort_values(by='species')
chroms = snakemake.params['chroms']
ploidy = snakemake.params['ploidy']
numbers = get_numbers_dict(ploidy)
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']
gffpath = snakemake.input['gff']

# load gff
features = allel.gff3_to_dataframe(gffpath,
                        attributes=["ID", "description"])
gff = features[features.type == 'gene']

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
    path = f"results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=samples,
                                                           numbers=numbers,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)
    # Store total SNPs per chromosome
    # And summaries from SnpEff
    total_snps_per_chrom[chrom] = geno.shape[0]
    snpeffdict[chrom] = snpeff[1].value_counts(normalize=True)

    ######## SNP counts per gene ########
    # Subset GFF to appropriate chromosome
    genes = gff[gff.seqid == f"{chrom}"].sort_values('start').reset_index(drop=True)

    snpsnotmissing = {}
    missing_array = {}
    snps_per_gene = {}
    snps_per_sample = {}

    for sample in samples['samples']:

        bool_ = sample == populations
        gn = geno.compress(bool_, axis=1)

        res = map(isnotmissing, gn[:,0])
        missing_array[sample] = np.fromiter(res, dtype=bool)
        snpsnotmissing[sample] = missing_array[sample].sum()

    for i, gene in genes.iterrows():
        ID = gene['ID']
        # locate_ranges() to get a boolean, needed as locate_range() will throw errors if no snps found in gene
        gene_bool = pos.locate_ranges([gene['start']], [gene['end']], strict=False)

        for sample in samples['samples']:

            presentSNPs = missing_array[sample].compress(gene_bool, axis=0)
            snps_per_sample[sample] = presentSNPs.sum()

        snps_per_gene[ID] = dict(snps_per_sample)

    snps_per_gene_allchroms[chrom] = pd.DataFrame.from_dict(flip_dict(snps_per_gene))
    
    
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
    pop_colours = get_colour_dict(treatment_indices['name'], "viridis")
    
    # Run PCA function defined in tools.py
    print(f"Performing PCA on {dataset} chromosome {chrom}")
    pca(geno, chrom, dataset, populations, samples, pop_colours, prune=True, scaler=None)

    ######## Plot variant density over genome (defined in tools.py) ########
    plot_density(pos, window_size=100000, title=f"Variant Density chromosome {chrom}", path=f"results/variants/plots/{dataset}_SNPdensity_{chrom}.png")

    #### Genome-wide statistics (seqDiv, Wattersons Theta, LD, inbreeding coefficient) ####
    seqdivdict = {}
    thetadict = {}
    coefdict= {}
    allcoef = defaultdict(list)

    for pop in samples.treatment.unique():

        # Sequence diversity 
        seqdivdict[pop] = allel.sequence_diversity(pos, acsubpops[pop])
        
        # Wattersons theta
        thetadict[pop] = allel.watterson_theta(pos, acsubpops[pop])

        # Inbreeding coefficient
        coef = allel.moving_statistic(gn, statistic=allel.inbreeding_coefficient, 
                                            size=1000, step=100)
        coef = np.nanmean(coef, axis=1)
        coefdict[pop] = np.mean(coef)
        allcoef[pop].append(np.array(coef))

        print(f"{pop},{chrom}, Sequence Diversity = ", seqdivdict[pop])
        print(f"{pop},{chrom}, Wattersons Theta = ", thetadict[pop])
        print(f"{pop},{chrom}, Inbreeding Coef = ", np.mean(coef), "\n")

    seqdivdictchrom[chrom] = dict(seqdivdict)
    thetadictchrom[chrom] = dict(thetadict)
    coefdictchrom[chrom] = dict(coefdict)

seqdivdictchrom= flip_dict(seqdivdictchrom)
thetadictchrom = flip_dict(thetadictchrom)
coefdictchrom = flip_dict(coefdictchrom)

# Get stats per chromosome
pd.DataFrame.from_dict(seqdivdictchrom).to_csv("results/variants/stats/SequenceDiversity.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(thetadictchrom).to_csv("results/variants/stats/WattersonsTheta.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(coefdictchrom).to_csv("results/variants/stats/inbreedingCoef.tsv", sep="\t", index=True)

# Get genome wide average stats
for pop in allcoef.keys():
    allcoef[pop] = np.nanmean(allcoef[pop])

coefdf = pd.DataFrame.from_dict(allcoef, orient='index', columns=['InbreedingCoefficient'])
coefdf.to_csv(f"results/variants/stats/inbreedingCoef.mean.tsv", sep="\t", index=True)

# Total SNPs per chrom, all samples
totalsnpsdf = pd.DataFrame.from_dict(total_snps_per_chrom, orient='index', columns=['Total_SNPs'])
totalsnpsdf.to_csv(f"results/variants/stats/totalSNPs.tsv", sep="\t", index=True)

# SNPcount per gene
snpcountsdf = pd.concat(snps_per_gene_allchroms)  
snpcountsdf.to_csv("results/variants/stats/nSNPsPerGene.tsv", sep="\t", index=True)
genesNoSNPs = pd.DataFrame((snpcountsdf == 0).sum(axis=0), columns=['Genes with zero SNPs'])
genesNoSNPs.to_csv("results/variants/stats/nGenesZeroSNPs.tsv", sep="\t", index=True)

# snpEff
snpeffdf = pd.concat(snpeffdict)
snpeffdf.to_csv("results/variants/stats/snpEffProportions.tsv", sep="\t", index=True)