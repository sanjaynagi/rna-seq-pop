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
gffpath = snakemake.input['gff']

# load gff
features = allel.gff3_to_dataframe(gffpath,
                        attributes=["ID", "Parent", "description"])

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
    path = f"results/variantAnalysis/vcfs/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations = readAndFilterVcf(path=path,
                                                           chrom=chrom,
                                                           samples=metadata,
                                                           numbers=numbers,
                                                           ploidy=ploidy,
                                                           qualflt=qualflt,
                                                           missingfltprop=missingprop)
    # Store total SNPs per chromosome
    # And summaries from SnpEff
    total_snps_per_chrom[chrom] = geno.shape[0]
    snpeffdict[chrom] = snpeff[1].value_counts(normalize=True)

    ######## SNP counts per gene ########
    # Subset GFF to appropriate chromosome
    gff = features.query("seqid == @chrom").sort_values('start').reset_index(drop=True)
    genes = gff.query("type == 'gene'")
    exons = features.query("type == 'exon'").reset_index(drop=True)
    
    ## Proportion SNPs per GFF feature
    snpsPerGff = getSNPGFFstats(gff,  pos)
    snpsPerGff.to_csv(snakemake.output['snpsPerGenomicFeature'], sep="\t", index=True)

    ## Calculate missing SNPs per sample, SNPs per gene etc
    snpsnotmissing = {}
    missing_array = {}
    snps_per_gene = {}
    snps_per_sample = {}

    # For each sample in metadata, compress the genotype array and check is data is missing
    for sample in metadata['sampleID']:

        bool_ = sample == populations
        gn = geno.compress(bool_, axis=1)

        res = map(isnotmissing, gn[:,0])
        missing_array[sample] = np.fromiter(res, dtype=bool)
        snpsnotmissing[sample] = missing_array[sample].sum()

    # For each gene, find out how many SNPs per gene 
    for i, gene in genes.iterrows():
        # locate_ranges() to get a boolean, needed as locate_range() will throw errors if no snps found in gene
        gene_bool = pos.locate_ranges([gene['start']], [gene['end']], strict=False)

        for sample in metadata['sampleID']:
            presentSNPs = missing_array[sample].compress(gene_bool, axis=0) # subset to each gene
            snps_per_sample[sample] = presentSNPs.sum()

        snps_per_gene[gene['ID']] = dict(snps_per_sample)

    snps_per_gene_allchroms[chrom] = pd.DataFrame.from_dict(flip_dict(snps_per_gene))

# Total SNPs per chrom, all samples
totalsnpsdf = pd.DataFrame.from_dict(total_snps_per_chrom, orient='index', columns=['Total_SNPs'])
totalsnpsdf.to_csv(f"results/variantAnalysis/stats/totalSNPs.tsv", sep="\t", index=True)

# SNPcount per gene
snpcountsdf = pd.concat(snps_per_gene_allchroms)  
snpcountsdf.to_csv("results/variantAnalysis/stats/nSNPsPerGene.tsv", sep="\t", index=True)
genesNoSNPs = pd.DataFrame((snpcountsdf == 0).sum(axis=0), columns=['Genes with zero SNPs'])
genesNoSNPs.to_csv("results/variantAnalysis/stats/nGenesZeroSNPs.tsv", sep="\t", index=True)

# snpEff summary
snpeffdf = pd.concat(snpeffdict)
snpeffdf.to_csv("results/variantAnalysis/stats/snpEffProportions.tsv", sep="\t", index=True)