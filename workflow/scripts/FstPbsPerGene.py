#!/usr/bin/env python3

"""
A script to perform Fst and population branch statistic analysis by gene on genotype data
"""

from tools import *
from scipy import stats
import warnings
warnings.filterwarnings('ignore') # suppress numpy runtime warnings, this is a bit dangerous, should be removed for release or resolve source of warnings

# snakemake inputs and params
metadata_path = snakemake.input['samples']
samples = pd.read_csv(metadata_path, sep="\s+")
gffpath = snakemake.input['gff']
comparisons_path = snakemake.input['DEcontrasts']
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
chroms = snakemake.params['chroms']
ploidy = snakemake.params['ploidy']
numbers = get_numbers_dict(ploidy)
missingprop = snakemake.params['missingprop']

# gff
features = allel.gff3_to_dataframe(gffpath,
                       attributes=["ID", "description"])
gff = features[features.type == 'gene']
# gene names file, rename
gene_names = pd.read_csv(snakemake.input['geneNames'], sep="\t")
gene_names.columns = ['GeneID' if x=='Gene_stable_ID' else x for x in gene_names.columns]

### main ####
# Read in list of contrasts
comparisons = pd.read_csv(comparisons_path)
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']
comparisons = [list(row) for i,row in comparisons.iterrows()]

print(f"The pairwise comparisons for Fst are {comparisons}")

fstbychrom={}
if pbs: pbsbychrom={}
tajdbychrom={}
gdivbychrom = {}

for chrom in chroms:

    # path to vcf
    path = f"results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
    # function to read in vcfs and associated SNP data
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops =  readAndFilterVcf(path=path,
                                                               chrom=chrom, samples=samples,
                                                               numbers=numbers,
                                                               qualflt=30,
                                                               missingfltprop=missingprop,
                                                               plot=False)
    # subset gff to appropriate chrom
    genes = gff[gff.seqid == chrom].sort_values('start').reset_index(drop=True)

    ### Average Fst, pbs, tajima d for each gene
    fst_per_comp = {}
    fst_per_gene = {}
    pbs_per_gene = {}
    pbs_per_comp = {}
    tajd_per_pop = {}
    tajd_per_gene = {}
    gdiv_per_pop = {}
    gdiv_per_gene = {}
    se_per_comp = {}
    se_per_gene = {}
    pos_dict = {}
    n_dict = {}

    # loop through each gene and calculate fst, pbs, tajimas d, or sequence diversity for each comparison
    for i, gene in genes.iterrows():
        ID = gene.ID
        # locate_ranges() to get a boolean, needed as locate_range() will throw errors if no snps found in gene
        gene_bool = pos.locate_ranges([gene['start']], [gene['end']], strict=False)
        nsnps = gene_bool.sum()

        # if there are less than 3 snps in this gene then skip to next in loop
        if nsnps < 2:
            continue

        # store number of snps per gene
        n_dict[ID] = nsnps
        # store midpoint positions of gene
        pos_dict[ID] = (gene['start'] + gene['end'])/2

        # fst per gene between each comparison
        for comp1,comp2 in comparisons:
            name = comp1 + "_" + comp2
            ac1 = acsubpops[comp1].compress(gene_bool, axis=0)
            ac2 = acsubpops[comp2].compress(gene_bool, axis=0)

            fst_per_comp[name], se_per_comp[name],_,_= allel.average_hudson_fst(ac1, ac2, blen=1)
         #  fst_per_comp[name] = stats.zscore(fst_per_comp[name], nan_policy='omit') Need to implement z-score transformation

        # tajimas d and sequence diversity per gene for each subpop(i.e treatment)
        for subpop in subpops:
            ac = acsubpops[subpop].compress(gene_bool)
            genepos = pos[gene_bool]
            tajd_per_pop[subpop] = allel.tajima_d(ac=ac, pos=genepos)
            gdiv_per_pop[subpop] = allel.sequence_diversity(ac=ac, pos=genepos)

        # pbs for each gene for each pbc comparison as defined in config.yaml
        if pbs is True:
            for pbscomp in pbscomps:
                name = pbscomp[0] + "_" + pbscomp[1] + "_" + pbscomp[2]
                pbs_per_comp[name],se,_,_ = meanPBS(acsubpops[pbscomp[0]].compress(gene_bool, axis=0),
                                          acsubpops[pbscomp[1]].compress(gene_bool, axis=0),
                                          acsubpops[pbscomp[2]].compress(gene_bool, axis=0),
                                                     window_size=1,
                                                    normalise=True)
        # store inner dict in outer dicts
        fst_per_gene[ID] = dict(fst_per_comp)
        se_per_gene[ID] = dict(se_per_comp)
        if pbs is True : pbs_per_gene[ID] = dict(pbs_per_comp)
        tajd_per_gene[ID] = dict(tajd_per_pop)
        gdiv_per_gene[ID] = dict(gdiv_per_pop)

    #reverse the dicts so the comparisons/subpops are on the outer dict
    fst_per_gene = flip_dict(fst_per_gene)
    se_per_gene = flip_dict(se_per_gene)
    if pbs is True : pbs_per_gene = flip_dict(pbs_per_gene)
    tajd_per_gene = flip_dict(tajd_per_gene)
    gdiv_per_gene = flip_dict(gdiv_per_gene)

    print(f"Chromosome {chrom} complete...\n")
    for comp1,comp2 in comparisons:
        name = comp1 + "_" + comp2
        a = np.array(list(fst_per_gene[name].values()))
        print(f"Overall Fst for chromosome {chrom} between {name} is {np.nanmean(a)}")

    #make dataframe of number of snps per gene (that pass quality and missingness filters)
    ndf = pd.DataFrame.from_dict(n_dict, orient='index').reset_index(drop=False)
    ndf.columns = ['GeneID', 'nSNPs']

    #make dataframe of midpoints of each gene
    posdf = pd.DataFrame.from_dict(pos_dict, orient='index').reset_index(drop=False)
    posdf.columns = ['GeneID', 'Gene_midpoint']

    #make dataframe of fst for each comparison
    fst_dfs = {}
    se_dfs = {}
    for comp1,comp2 in comparisons:
        name = comp1 + "_" + comp2
        fst_df = pd.DataFrame.from_dict(fst_per_gene[name], orient='index').reset_index(drop=False)
        fst_df.columns = ['GeneID', (name + '_zFst')]
        fst_dfs[name] = fst_df

        se_df = pd.DataFrame.from_dict(se_per_gene[name], orient='index').reset_index(drop=False)
        se_df.columns = ['GeneID', (name + '_SE')]
        se_dfs[name] = se_df

    my_reduce = partial(pd.merge, on='GeneID', how='outer')
    fst_allcomparisons = reduce(my_reduce, fst_dfs.values())
    se_allcomparisons = reduce(my_reduce, se_dfs.values())
    fst_allcomparisons = reduce(my_reduce, [fst_allcomparisons, se_allcomparisons])
    fst_allcomparisons['chrom'] = chrom

    tajd_dfs = {}
    gdiv_dfs = {}
    #store sequence diversityt and tajimas d for each gene and each subpop
    for subpop in subpops:
        tajd_df = pd.DataFrame.from_dict(tajd_per_gene[subpop], orient='index').reset_index(drop=False)
        tajd_df.columns = ['GeneID', (subpop+"_Tajima_d")]
        tajd_dfs[subpop] = tajd_df
        gdiv_df = pd.DataFrame.from_dict(gdiv_per_gene[subpop], orient='index').reset_index(drop=False)
        gdiv_df.columns = ['GeneID', (subpop+"_SeqDiv")]
        gdiv_dfs[subpop] = gdiv_df

    #combine tajimas d and sequence diversity for each sample
    tajdall = reduce(my_reduce, tajd_dfs.values())
    gdivall = reduce(my_reduce, gdiv_dfs.values())
    tajdall['chrom'] = chrom
    gdivall['chrom'] = chrom

    if pbs is True:
        #pbs store as dataframes 
        pbs_dfs = {}
        for pbscomp in pbscomps:
            name = pbscomp[0] + "_" + pbscomp[1] + "_" + pbscomp[2]
            pbs_df = pd.DataFrame.from_dict(pbs_per_gene[name], orient='index').reset_index(drop=False)
            pbs_df.columns = ['GeneID', (name+"PBS")]
            pbs_dfs[name] = pbs_df

        pbs_allcomparisons = reduce(my_reduce, pbs_dfs.values())
        pbs_allcomparisons['chrom'] = chrom

    fstbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), [fst_allcomparisons, gene_names, ndf, posdf])
    tajdbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), [tajdall, gene_names, ndf,posdf])
    gdivbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), [gdivall, gene_names, ndf, posdf])
    if pbs is True:
        pbsbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                    how='inner'), [pbs_allcomparisons, gene_names, ndf, posdf])
    



fstall = pd.concat(fstbychrom.values(), ignore_index=True)
tajdall = pd.concat(tajdbychrom.values(), ignore_index=True)
gdivall = pd.concat(gdivbychrom.values(), ignore_index=True)

#write to csv
fstall.to_csv(f"results/variants/Fst.tsv", index=False, sep="\t")
tajdall.to_csv(f"results/variants/TajimasD.tsv", index=False, sep="\t")
gdivall.to_csv(f"results/variants/SequenceDiv.tsv", index=False, sep="\t")

if pbs is True:
    pbsall = pd.concat(pbsbychrom.values(), ignore_index=True)
    pbsall.to_csv(f"results/variants/PBS.tsv", index=False, sep="\t")
