#!/usr/bin/env python3

"""
A script to perform Fst and population branch statistic analysis by gene on genotype data
"""

from tools import *
from tqdm import tqdm

# snakemake inputs and params
metadata_path = snakemake.input[0]
samples = pd.read_csv(metadata_path, sep="\s+")
gffpath = snakemake.input[1]
comparisons_path = snakemake.input[2]
pbs = snakemake.params[0]
pbscomps = snakemake.params[1]
chroms = snakemake.params[2]
ploidy = snakemake.params[3]
missingprop = snakemake.params[4]
gffchromprefix= snakemake.params[5]

#gff
features = allel.gff3_to_dataframe(gffpath,
                       attributes=["ID", "description"])
gff = features[features.type == 'gene']
# gene names file, rename 
gene_names = pd.read_csv(snakemake.input[3], sep="\t")
gene_names.columns = ['GeneID' if x=='Gene_stable_ID' else x for x in gene_names.columns]

### main ####
# Read in list of contrasts
comparisons = pd.read_csv(comparisons_path, header=None)
comparisons.columns = ['contrast']
comparisons = comparisons.contrast.str.split("_", expand=True)
comparisons.columns = ['sus', 'res']

comparisons = [list(row) for i,row in comparisons.iterrows()]
print(f"The pairwise comparisons for Fst are {comparisons}")


fstpbsbychrom=dict()
tajdbychrom=dict()
gdivbychrom = dict()

for chrom in chroms:
    
    #path to vcf
    path = f"results/variants/annot.variants.{chrom}.vcf.gz"
    #function to read in vcfs and associated SNP data
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops =  readAndFilterVcf(path=path, 
                                                               chrom=chrom, 
                                                               qualflt=30, 
                                                               plot=False)
    #subset gff to appropriate chrom
    genes = gff[gff.seqid == f"{gffchromprefix}{chrom}"].sort_values('start').reset_index(drop=True)

    ### Average Fst, pbs, tajima d for each gene 
    fst_per_comp = dict()
    fst_per_gene = dict()
    pbs_per_gene = dict()
    pbs_per_comp = dict()
    tajd_per_pop = dict()
    tajd_per_gene = dict()
    gdiv_per_pop = dict()
    gdiv_per_gene = dict()
    se_per_comp = dict()
    se_per_gene = dict()
    pos_dict = dict()
    n_dict = dict()

    #loop through each gene and calculate fst, pbs, tajimas d, or sequence diversity for each comparison
    for i, gene in tqdm(genes.iterrows()):
        ID = gene.ID
        #locate_ranges() to get a boolean, needed as locate_range() will throw errors if no snps found in gene
        gene_bool = pos.locate_ranges([gene['start']], [gene['end']], strict=False)
        nsnps = gene_bool.sum()

        #if there are less than 3 snps in this gene then skip to next in loop
        if nsnps < 3:
            continue

        #store number of snps per gene
        n_dict[ID] = nsnps
        #store midpoint positions of gene
        pos_dict[ID] = (gene['start'] + gene['end'])/2

        #fst per gene between each comparison
        for comp1,comp2 in comparisons:
            name = comp1 + "_" + comp2
            ac1 = acsubpops[comp1].compress(gene_bool, axis=0)
            ac2 = acsubpops[comp2].compress(gene_bool, axis=0)

            fst_per_comp[name], se_per_comp[name],_,_= allel.average_hudson_fst(ac1, ac2, blen=1)   

        #tajimas d and sequence diversity per gene for each subpop(i.e treatment)
        for subpop in subpops:
            ac = acsubpops[subpop].compress(gene_bool)
            genepos = pos[gene_bool]
            tajd_per_pop[subpop] = allel.tajima_d(ac=ac, pos=genepos)
            gdiv_per_pop[subpop] = allel.sequence_diversity(ac=ac, pos=genepos)

        #pbs for each gene for each pbc comparison as defined in config.yaml
        if pbs is True:
            for pbscomp in pbscomps:
                name = pbscomp[0] + "_" + pbscomp[1] + "_" + pbscomp[2]
                pbs_per_comp[name],se,_,_ = meanPBS(acsubpops[pbscomp[0]].compress(gene_bool, axis=0),
                                          acsubpops[pbscomp[1]].compress(gene_bool, axis=0), 
                                          acsubpops[pbscomp[2]].compress(gene_bool, axis=0), 
                                                     window_size=1,
                                                    normalise=True)
        #store inner dict in outer dicts
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
    fst_dfs = dict()
    se_dfs = dict()
    for comp1,comp2 in comparisons:
        name = comp1 + "_" + comp2
        fst_df = pd.DataFrame.from_dict(fst_per_gene[name], orient='index').reset_index(drop=False)
        fst_df.columns = ['GeneID', (name + '_Fst')]
        fst_dfs[name] = fst_df

        se_df = pd.DataFrame.from_dict(se_per_gene[name], orient='index').reset_index(drop=False)
        se_df.columns = ['GeneID', (name + '_SE')]
        se_dfs[name] = se_df
    
    my_reduce = partial(pd.merge, on='GeneID', how='outer')
    fst_allcomparisons = reduce(my_reduce, fst_dfs.values())
    se_allcomparisons = reduce(my_reduce, se_dfs.values())
    fst_allcomparisons = reduce(my_reduce, [fst_allcomparisons, se_allcomparisons])
    fst_allcomparisons['chrom'] = chrom

    tajd_dfs = dict()
    gdiv_dfs = dict()
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
        pbs_dfs = dict()
        for pbscomp in pbscomps:
            name = pbscomp[0] + "_" + pbscomp[1] + "_" + pbscomp[2]
            pbs_df = pd.DataFrame.from_dict(pbs_per_gene[name], orient='index').reset_index(drop=False)
            pbs_df.columns = ['GeneID', (name+"PBS")]
            pbs_dfs[name] = pbs_df

        pbs_allcomparisons = reduce(my_reduce, pbs_dfs.values())  
        dfs = [fst_allcomparisons,pbs_allcomparisons, gene_names, ndf, posdf]

    dfs = [fst_allcomparisons, gene_names, ndf, posdf]
    
    tajdbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), [tajdall, gene_names, ndf,posdf])
    gdivbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), [gdivall, gene_names, ndf, posdf])
    fstpbsbychrom[chrom] = reduce(lambda  left,right: pd.merge(left,right,on=['GeneID'],
                                                how='inner'), dfs)

fstpbsall = pd.concat(fstpbsbychrom.values(), ignore_index=True)
tajdall = pd.concat(tajdbychrom.values(), ignore_index=True)
gdivall = pd.concat(gdivbychrom.values(), ignore_index=True)

#write to csv
if pbs is True:
    fstpbsall.to_csv(f"results/variants/Fst_PBS.tsv", index=False, sep="\t")
else:
    fstpbsall.to_csv(f"results/variants/Fst.tsv", index=False, sep="\t")

tajdall.to_csv(f"results/variants/tajimas_d.tsv", index=False, sep="\t")
gdivall.to_csv(f"results/variants/sequence_div.tsv", index=False, sep="\t")