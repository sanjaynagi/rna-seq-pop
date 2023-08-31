#!/usr/bin/env python3

"""
A script to perform Fst and population branch statistic analysis by gene on genotype data
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import rnaseqpoptools as rnaseqpop
import pandas as pd
import numpy as np
import allel
from scipy import stats
from functools import partial, reduce
import warnings
warnings.filterwarnings('ignore') # suppress numpy runtime warnings, this is a bit dangerous, should be removed for release or resolve source of warnings

# snakemake inputs and params
dataset = snakemake.params['dataset']
metadata_path = snakemake.input['metadata']
metadata = rnaseqpop.load_metadata(metadata_path)
gffpath = snakemake.input['gff']
pbs = snakemake.params['pbs']
pbscomps = snakemake.params['pbscomps']
contigs = snakemake.params['contigs']
ploidy = snakemake.params['ploidy']
numbers = rnaseqpop.get_numbers_dict(ploidy)
missingprop = snakemake.params['missingprop']

# gff
features = allel.gff3_to_dataframe(gffpath,
                       attributes=["ID", "description"])
gff = features[features.type == 'gene']
# gene names file, rename
gene_names = pd.read_csv(snakemake.input['geneNames'], sep="\t").drop(columns=['TranscriptID']).drop_duplicates()

### main ####
# Read in list of contrasts
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']
comparisons = [list(row) for i,row in comparisons.iterrows()]

print(f"The pairwise comparisons for Fst are {comparisons}")

fstbychrom={}
if pbs: pbsbychrom={}
tajdbychrom={}
gdivbychrom = {}
dxybychrom = {}

for contig in contigs:

    # path to vcf
    path = f"results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz"
    # function to read in vcfs and associated SNP data
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops = rnaseqpop.readAndFilterVcf(path=path,
                                                               contig=contig, 
                                                               samples=metadata,
                                                               numbers=numbers,
                                                               ploidy=ploidy,
                                                               qualflt=30,
                                                               missingfltprop=missingprop)
    # subset gff to appropriate contig
    genes = gff[gff.seqid == contig].sort_values('start').reset_index(drop=True)

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
    dxy_per_comp = {}
    dxy_per_gene = {}
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

        # fst and dxy per gene between each comparison
        for comp1,comp2 in comparisons:
            name = comp1 + "_" + comp2
            ac1 = acsubpops[comp1].compress(gene_bool, axis=0)
            ac2 = acsubpops[comp2].compress(gene_bool, axis=0)

            fst_per_comp[name], se_per_comp[name],_,_= allel.average_hudson_fst(ac1, ac2, blen=1)
            
            dxy_per_comp[name] = allel.sequence_divergence(pos[gene_bool], ac1, ac2)

        # tajimas d and sequence diversity per gene for each subpop(i.e treatment)
        for subpop in subpops:
            ac = acsubpops[subpop].compress(gene_bool)
            genepos = pos[gene_bool]
            tajd_per_pop[subpop] = allel.tajima_d(ac=ac, pos=genepos)
            gdiv_per_pop[subpop] = allel.sequence_diversity(ac=ac, pos=genepos)

        # pbs for each gene for each pbc comparison as defined in config.yaml
        if pbs is True:
            for pbscomp in pbscomps:
                pop1, pop2, outpop = pbscomp.split("_")
                pbs_per_comp[pbscomp],se,_,_ = rnaseqpop.meanPBS(acsubpops[pop1].compress(gene_bool, axis=0),
                                          acsubpops[pop2].compress(gene_bool, axis=0),
                                          acsubpops[outpop].compress(gene_bool, axis=0),
                                                     window_size=1,
                                                    normalise=True)
        # store inner dict in outer dicts
        fst_per_gene[ID] = dict(fst_per_comp)
        se_per_gene[ID] = dict(se_per_comp)
        if pbs is True : pbs_per_gene[ID] = dict(pbs_per_comp)
        tajd_per_gene[ID] = dict(tajd_per_pop)
        gdiv_per_gene[ID] = dict(gdiv_per_pop)
        dxy_per_gene[ID] = dict(dxy_per_comp)

    #reverse the dicts so the comparisons/subpops are on the outer dict
    fst_per_gene = rnaseqpop.flip_dict(fst_per_gene)
    se_per_gene = rnaseqpop.flip_dict(se_per_gene)
    if pbs is True : pbs_per_gene = rnaseqpop.flip_dict(pbs_per_gene)
    tajd_per_gene = rnaseqpop.flip_dict(tajd_per_gene)
    gdiv_per_gene = rnaseqpop.flip_dict(gdiv_per_gene)
    dxy_per_gene = rnaseqpop.flip_dict(dxy_per_gene)

    print(f"Chromosome {contig} complete...\n")
    for comp1,comp2 in comparisons:
        name = comp1 + "_" + comp2
        a = np.array(list(fst_per_gene[name].values()))
        print(f"Overall Fst (averaged across genes) for chromosome {contig} between {name} is {np.nanmean(a)}")

    # Make dataframe of number of snps per gene (that pass quality and missingness filters)
    ndf = pd.DataFrame.from_dict(n_dict, orient='index').reset_index(drop=False)
    ndf.columns = ['GeneID', 'nSNPs']

    # Make dataframe of midpoints of each gene
    posdf = pd.DataFrame.from_dict(pos_dict, orient='index').reset_index(drop=False)
    posdf.columns = ['GeneID', 'Gene_midpoint']

    # Make dataframe of fst for each comparison
    fst_dfs = {}
    se_dfs = {}
    dxy_dfs = {}
    for comp1,comp2 in comparisons:
        name = comp1 + "_" + comp2
        fst_df = pd.DataFrame.from_dict(fst_per_gene[name], orient='index').reset_index(drop=False)
        fst_df.columns = ['GeneID', (name + '_zFst')]
        fst_dfs[name] = fst_df

        se_df = pd.DataFrame.from_dict(se_per_gene[name], orient='index').reset_index(drop=False)
        se_df.columns = ['GeneID', (name + '_SE')]
        se_dfs[name] = se_df

        dxy_df = pd.DataFrame.from_dict(dxy_per_gene[name], orient='index').reset_index(drop=False)
        dxy_df.columns = ['GeneID', (name + '_dxy')]
        dxy_dfs[name] = dxy_df


    my_reduce = partial(pd.merge, on='GeneID', how='outer')
    fst_allcomparisons = reduce(my_reduce, fst_dfs.values())
    se_allcomparisons = reduce(my_reduce, se_dfs.values())
    fst_allcomparisons = reduce(my_reduce, [fst_allcomparisons, se_allcomparisons])
    fst_allcomparisons['contig'] = contig
    dxy_allcomparisons = reduce(my_reduce, dxy_dfs.values())
    dxy_allcomparisons['contig'] = contig

    tajd_dfs = {}
    gdiv_dfs = {}
    # Store sequence diversity and tajimas d for each gene and each subpop
    for subpop in subpops:
        tajd_df = pd.DataFrame.from_dict(tajd_per_gene[subpop], orient='index').reset_index(drop=False)
        tajd_df.columns = ['GeneID', (subpop+"_Tajima_d")]
        tajd_dfs[subpop] = tajd_df
        gdiv_df = pd.DataFrame.from_dict(gdiv_per_gene[subpop], orient='index').reset_index(drop=False)
        gdiv_df.columns = ['GeneID', (subpop+"_SeqDiv")]
        gdiv_dfs[subpop] = gdiv_df

    # Combine tajimas d and sequence diversity for each sample
    tajdall = reduce(my_reduce, tajd_dfs.values())
    gdivall = reduce(my_reduce, gdiv_dfs.values())
    tajdall['contig'] = contig
    gdivall['contig'] = contig

    if pbs is True:
        # PBS store as dataframes 
        pbs_dfs = {}
        for pbscomp in pbscomps:
            pbs_df = pd.DataFrame.from_dict(pbs_per_gene[pbscomp], orient='index').reset_index(drop=False)
            pbs_df.columns = ['GeneID', (pbscomp+"PBS")]
            pbs_dfs[pbscomp] = pbs_df

        pbs_allcomparisons = reduce(my_reduce, pbs_dfs.values())
        pbs_allcomparisons['contig'] = contig

    fstbychrom[contig] = reduce(lambda left, right: pd.merge(left,right, on=['GeneID'],
                                                how='inner'), [fst_allcomparisons, gene_names, ndf, posdf])
    tajdbychrom[contig] = reduce(lambda left, right: pd.merge(left,right, on=['GeneID'],
                                                how='inner'), [tajdall, gene_names, ndf,posdf])
    gdivbychrom[contig] = reduce(lambda left, right: pd.merge(left,right, on=['GeneID'],
                                                how='inner'), [gdivall, gene_names, ndf, posdf])
    dxybychrom[contig] = reduce(lambda left, right: pd.merge(left, right, on=['GeneID'],
                                                how='inner'), [dxy_allcomparisons, gene_names, ndf, posdf])
    if pbs is True:
        pbsbychrom[contig] = reduce(lambda left,right: pd.merge(left,right,on=['GeneID'],
                                                    how='inner'), [pbs_allcomparisons, gene_names, ndf, posdf])
    

## Concatenate chromosome dfs to one big data frame and remove duplicates
fstall = pd.concat(fstbychrom.values(), ignore_index=True).drop_duplicates()
tajdall = pd.concat(tajdbychrom.values(), ignore_index=True).drop_duplicates()
gdivall = pd.concat(gdivbychrom.values(), ignore_index=True).drop_duplicates()
dxyall = pd.concat(dxybychrom.values(), ignore_index=True).drop_duplicates()

# Write to csv
fstall.to_csv(f"results/variantAnalysis/selection/FstPerGene.tsv", index=False, sep="\t")
tajdall.to_csv(f"results/variantAnalysis/selection/TajimasDPerGene.tsv", index=False, sep="\t")
gdivall.to_csv(f"results/variantAnalysis/diversity/SequenceDivPerGene.tsv", index=False, sep="\t")
dxyall.to_csv(f"results/variantAnalysis/diversity/DxyPerGene.tsv", index=False, sep="\t")

if pbs is True:
    pbsall = pd.concat(pbsbychrom.values(), ignore_index=True).drop_duplicates()
    pbsall.to_csv(f"results/variantAnalysis/selection/PbsPerGene.tsv", index=False, sep="\t")
