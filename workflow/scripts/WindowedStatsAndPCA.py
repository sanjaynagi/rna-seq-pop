#!/usr/bin/env python3

"""
A script to calculate various windowed population genetic statistics and PCA 
"""

import pandas as pd
import numpy as np
import allel
from scipy import stats
import seaborn as sns
import zarr
import matplotlib
import matplotlib.pyplot as plt
from tools import *

dataset = snakemake.params['dataset']
samples = pd.read_csv(snakemake.input['samples'], sep="\t")
samples = samples.sort_values(by='species')
chroms = snakemake.params['chroms']
comparisons_path = snakemake.input['comparisons']
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']
gffpath = snakemake.input['gff']


# Read in list of contrasts
comparisons = pd.read_csv(comparisons_path)
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']
comparisons = [list(row) for i,row in comparisons.iterrows()]

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
snpeff = {}

for chrom in chroms:

    path = f"results/variants/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=samples,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)
    
    total_snps_per_chrom[chrom] = geno.shape[0]
    snpeff[chrom] = snpeff[1].value_counts(normalize=True)

    ######## SNP counts per gene ########
    # subset gff to appropriate chrom
    genes = gff[gff.seqid == f"{chrom}"].sort_values('start').reset_index(drop=True)

    snpsnotmissing = {}
    missing_array = {}
    snps_per_gene = {}
    snps_per_sample = {}

    for sample in samples['samples']:

        bool_ = sample == pops
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
    
    ######## pbs windowed ########
    if pbs:
        for pbscomp in pbscomps:
            name = pbscomp[0] + "_" + pbscomp[1] + "_" + pbscomp[2]

            print(f"Calculating PBS values in sliding window for {name}\n")

            pbs = allel.pbs(acsubpops[pbscomp[0]], 
                            acsubpops[pbscomp[1]], 
                            acsubpops[pbscomp[2]], 
                            window_size=1000, window_step=100, normed=True)
            midpoint = allel.moving_statistic(pos, np.mean, 1000, step=100)

            plt.figure(figsize=[20,8])
            sns.lineplot(midpoint, pbs)
            plt.title(f"PBS {chrom} {name}")
            plt.savefig(f"results/variants/plots/PBS_{name}.{chrom}.line.png")
            plt.close()
            plt.figure()
            sns.scatterplot(midpoint, pbs)
            plt.title(f"PBS {chrom} {name}")
            plt.savefig(f"results/variants/plots/PBS_{name}.{chrom}.scatter.png")


    #### pattersons f3 statistic ####
    #pattersonf3(acsubpops['gambiaeCont'], acsubpops['coluzziiCont'], acsubpops['coluzziiDelta'], pos, f"gambcont_{chrom}")
    #pattersonf3(acsubpops['gambiaeDelta'], acsubpops['coluzziiCont'], acsubpops['coluzziiDelta'], pos, f"gambdelta_{chrom}")
    #pattersonf3(acsubpops['coluzziiDelta'], acsubpops['gambiaeCont'], acsubpops['gambiaeDelta'], pos, f"coludelta_{chrom}")
    

    #### PCA #####
    d={}
    for name, inds in subpops.items():
        for n in range(4):
            p = inds[n]
            d[p] = name

    treatment_indices = pd.DataFrame.from_dict(d, orient='index')
    treatment_indices.reset_index(inplace=True)
    treatment_indices = treatment_indices.rename(columns = {'index':'sample_index', 0:"name"})

    pop_colours = get_colour_dict(treatment_indices['name'], "viridis")

    print(f"\n Performing PCA on {dataset} chromosome {chrom}")
    pca(geno, chrom, dataset, treatment_indices, prune=True, scaler=None)

    ######## Variant density ####
    plot_density(pos, window_size=100000, title=f"Variant Density chromosome {chrom}", path=f"results/variants/plots/{dataset}_SNPdensity_{chrom}.png")

    #### genome-wide mean statistics (seqDiv, LD, inbreeding coefficient) ####
    coefdict= {}
    seqdivdict = {}
    ldict = {}

    coefdictchrom= {}
    seqdivdictchrom = {}
    ldictchrom = {}
    allcoef = defaultdict(list)
    allld = defaultdict(list)

    for pop in samples.treatment.unique():
        
        gn = geno.take(subpops[pop], axis=1)
        bial_ = acsubpops[pop].is_biallelic()
        gnalt = gn.compress(bial_, axis=0).to_n_alt()
        
        # inbreeding coefficient
        coef = allel.moving_statistic(gn,statistic=allel.inbreeding_coefficient, 
                                            size=1000, step=100)
        coef = np.nanmean(coef, axis=1)
        coefdict[pop] = np.mean(coef)
        allcoef[pop].append(np.array(coef))

        # linkage
        ld = allel.rogers_huff_r(gnalt)
        ldict[pop] = np.mean(ld)
        allld[pop].append(ld)

        # sequence diversity 
        seqdiv = allel.sequence_diversity(pos, acsubpops[pop])
        seqdivdict[pop] = seqdiv


        print(f"\n{pop}, {chrom}, inbreeding coef = ", np.mean(coef))
        print(f"{pop},{chrom}, ld (rogers huff r2) = ", ld, "\n")
        print(f"{pop},{chrom}, sequence diversity = ", seqdiv, "\n")

    coefdictchrom[chrom] = dict(coefdict)
    seqdivdictchrom[chrom] = dict(seqdivdict)
    ldictchrom[chrom] = dict(ldict)

coefdictchrom = flip_dict(coefdictchrom)
seqdivdictchrom= flip_dict(seqdivdictchrom)
ldictchrom = flip_dict(ldictchrom)

#get AIM fractions per chromosome
pd.DataFrame.from_dict(coefdictchrom).to_csv("results/variants/stats/inbreedingCoef.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(seqdivdictchrom).to_csv("results/variants/stats/SequenceDiversity.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(ldictchrom).to_csv("results/variants/stats/LD.tsv", sep="\t", index=True)

# get genome wide average AIM fractions
for k in allcoef.keys():
    allld[k] = np.nanmean(allld[k])
    allcoef[k] = np.nanmean(allcoef[k])

df1 = pd.DataFrame.from_dict(allld, orient='index',columns=['LinkageDisequilibrium'])
df2 = pd.DataFrame.from_dict(allcoef, orient='index', columns=['InbreedingCoefficient'])

df1.to_csv(f"results/variants/stats/LD.mean.tsv", sep="\t", index=True)
df2.to_csv(f"results/variants/stats/inbreedingCoef.mean.tsv", sep="\t", index=True)

### Total SNPs per chrom, all samples
total_snps_df = pd.DataFrame.from_dict(total_snps_per_chrom, orient='index', columns=['Total_SNPs'])
total_snps_df.to_csv(f"results/variants/stats/totalSNPs.tsv", sep="\t", index=True)

# SNPcount per gene
snpcountsdf = pd.concat(snps_per_gene_allchroms)  
snpcountsdf.to_csv("results/variants/stats/nSNPsPerGene.tsv", sep="\t", index=True)
genesNoSNPs = pd.DataFrame((snpcountsdf == 0).sum(axis=0), columns=['Genes with zero SNPs'])
genesNoSNPs.to_csv("results/variants/stats/nGenesZeroSNPs.tsv", sep="\t", index=True)

# snpEff
snpeffdf = pd.concat(snpeff)
snpeffdf.to_csv("results/variants/stats/snpEffProportions.tsv", sep="\t", index=True)