#!/usr/bin/env python3

"""
A script to calculate various SNP statistics, windowed population genetic statistics (Fst, PBS) and PCA.
Currently not modularised to reduce the repetition of loading and filtering VCFs (which is slow). 
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

from tools import *

dataset = snakemake.params['dataset']
samples = pd.read_csv(snakemake.input['samples'], sep="\t")
samples = samples.sort_values(by='species')
chroms = snakemake.params['chroms']
ploidy = snakemake.params['ploidy']
numbers = get_numbers_dict(ploidy)
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']
gffpath = snakemake.input['gff']
linkage = snakemake.params['linkage']

#Fst/PBS window size
windowsizes = snakemake.params['windowsizes']
windowsteps = snakemake.params['windowsteps']
windownames = snakemake.params['windownames']

# Read in list of contrasts
comparisons = pd.read_csv(snakemake.input['contrasts'])
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
snpeffdict = {}

seqdivdictchrom = {}
thetadictchrom = {}
coefdictchrom= {}
ldictchrom = {}

for i, chrom in enumerate(chroms):

    path = f"results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=samples,
                                                           numbers=numbers,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)

    total_snps_per_chrom[chrom] = geno.shape[0]
    snpeffdict[chrom] = snpeff[1].value_counts(normalize=True)

    ######## SNP counts per gene ########
    # subset gff to appropriate chrom
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

    print("comp len:", len(comparisons))
    ######## Fst in windows ########## 
    for sus, res in comparisons:
        name = sus + "_" + res
        print(f"Calculating Fst values in sliding windows for {name}\n")

        for wname, size, step in zip(windownames, windowsizes, windowsteps):
            FstArray = allel.moving_hudson_fst(acsubpops[sus], 
                            acsubpops[res], 
                            size=size, step=step)
            midpoint = allel.moving_statistic(pos, np.mean, size=size, step=step)

            plt.figure(figsize=[20,8])
            sns.lineplot(midpoint, FstArray)
            plt.title(f"Fst {chrom} {name}")
            plt.savefig(f"results/variants/plots/fst/{name}.{chrom}.fst.{wname}.png")
            plt.close()
            
        

    ######## Population Branch Statistic (PBS) in windows ########
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
                plt.savefig(f"results/variants/plots/pbs/{pbscomp}.{chrom}.pbs.{wname}.png")
                plt.close()


    ######## Principal Components Analysis (PCA) ########
    d={}
    for name, inds in subpops.items():
        for n in range(len(list(subpops.values())[0])):
            p = inds[n]
            d[p] = name

    treatment_indices = pd.DataFrame.from_dict(d, orient='index')
    treatment_indices.reset_index(inplace=True)
    treatment_indices = treatment_indices.rename(columns = {'index':'sample_index', 0:"name"})

    pop_colours = get_colour_dict(treatment_indices['name'], "viridis")
    print(f"Performing PCA on {dataset} chromosome {chrom}")
    pca(geno, chrom, dataset, populations, samples, pop_colours, prune=True, scaler=None)

    ######## Variant density over genome ########
    plot_density(pos, window_size=100000, title=f"Variant Density chromosome {chrom}", path=f"results/variants/plots/{dataset}_SNPdensity_{chrom}.png")

    #### Genome-wide statistics (seqDiv, Wattersons Theta, LD, inbreeding coefficient) ####
    seqdivdict = {}
    thetadict = {}
    coefdict= {}
    ldict = {}
    allcoef = defaultdict(list)
    allld = defaultdict(list)

    for pop in samples.treatment.unique():

        # subset to each population and filter to biallelic markers for LD calculation
        gn = geno.take(subpops[pop], axis=1)
        bial_ = acsubpops[pop].is_biallelic()
        gnalt = gn.compress(bial_, axis=0).to_n_alt()

        # sequence diversity 
        seqdivdict[pop] = allel.sequence_diversity(pos, acsubpops[pop])
        
        # wattersons theta
        thetadict[pop] = allel.watterson_theta(pos, acsubpops[pop])

        # inbreeding coefficient
        coef = allel.moving_statistic(gn, statistic=allel.inbreeding_coefficient, 
                                            size=1000, step=100)
        coef = np.nanmean(coef, axis=1)
        coefdict[pop] = np.mean(coef)
        allcoef[pop].append(np.array(coef))

        # linkage
        if linkage:
            ld = allel.rogers_huff_r(gnalt)
            allld[pop].append(ld)
	    # remove nan and infs to calculate average LD
            ld = ld[~np.logical_or(np.isinf(ld), np.isnan(ld))]
            ldict[pop] = np.nanmean(ld)

        print(f"{pop},{chrom}, Sequence Diversity = ", seqdivdict[pop])
        print(f"{pop},{chrom}, Wattersons Theta = ", thetadict[pop])
        print(f"{pop},{chrom}, Inbreeding Coef = ", np.mean(coef), "\n")
        if linkage is True: print(f"{pop},{chrom}, LD (rogers huff r2) = ", np.nanmean(ld))

    seqdivdictchrom[chrom] = dict(seqdivdict)
    thetadictchrom[chrom] = dict(thetadict)
    coefdictchrom[chrom] = dict(coefdict)
    if linkage: ldictchrom[chrom] = dict(ldict)

seqdivdictchrom= flip_dict(seqdivdictchrom)
thetadictchrom = flip_dict(thetadictchrom)
coefdictchrom = flip_dict(coefdictchrom)
if linkage: ldictchrom = flip_dict(ldictchrom)

#get stats per chromosome
pd.DataFrame.from_dict(seqdivdictchrom).to_csv("results/variants/stats/SequenceDiversity.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(thetadictchrom).to_csv("results/variants/stats/WattersonsTheta.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(coefdictchrom).to_csv("results/variants/stats/inbreedingCoef.tsv", sep="\t", index=True)
if linkage: pd.DataFrame.from_dict(ldictchrom).to_csv("results/variants/stats/LD.tsv", sep="\t", index=True)

# get genome wide average stats
for pop in allcoef.keys():
    allld[pop] = np.nanmean(allld[pop])
    allcoef[pop] = np.nanmean(allcoef[pop])

if linkage:
    df1 = pd.DataFrame.from_dict(allld, orient='index',columns=['LinkageDisequilibrium'])
    df1.to_csv(f"results/variants/stats/LD.mean.tsv", sep="\t", index=True)

df2 = pd.DataFrame.from_dict(allcoef, orient='index', columns=['InbreedingCoefficient'])
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
snpeffdf = pd.concat(snpeffdict)
snpeffdf.to_csv("results/variants/stats/snpEffProportions.tsv", sep="\t", index=True)
