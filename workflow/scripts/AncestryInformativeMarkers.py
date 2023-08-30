#!/usr/bin/env python3

"""
A script to determine the proportion of Anopheles gambiae/coluzzii alleles at Ancestry Informative Markers across the genome
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import rnaseqpoptools as rnaseqpop
import pandas as pd
import numpy as np
import zarr
import allel
from collections import defaultdict

### AIMS ###
dataset = snakemake.params['dataset']
metadata = rnaseqpop.load_metadata(snakemake.input['metadata'])
metadata = metadata.sort_values(by='species').reset_index(drop=True)
contigs = snakemake.params['contigs']
ploidy = snakemake.params['ploidy']
numbers = rnaseqpop.get_numbers_dict(ploidy)
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']

# read AIMs
aims = zarr.open(snakemake.input['aims_zarr_gambcolu'], mode='r')

## initialize dicts
ancestryPerAim = {}
aims_chrom_gamb = {}
aims_chrom_colu = {}
all_gamb = defaultdict(list)
all_colu = defaultdict(list)
n_aims_per_chrom = {}

for contig in contigs:

    # read in and filter data
    path = f"results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops =  rnaseqpop.readAndFilterVcf(path=path,
                                                               contig=contig,
                                                               samples=metadata,
                                                               numbers=numbers,
                                                               ploidy=ploidy,
                                                               qualflt=qualflt,
                                                               missingfltprop=missingprop)
    aimspos = aims[contig]['POS'][:]

    # get intersection of aims and our SNPs
    aims_pos_mask, aims_mask_2 = pos.locate_intersection(aimspos)
    our_aims = pos[aims_pos_mask]
    print(f"\n In the data, across all samples there are {our_aims.shape[0]} Ancestry Informative markers on Chromosome {contig}")

    # get gamb and colu alleles, and subset to aims that we have in the rna-seq data 
    aimscolu = aims[contig]['colu_allele'][:][aims_mask_2]
    aimsgamb = aims[contig]['gamb_allele'][:][aims_mask_2]

    # get mask that was used in readAndFilterVcf()
    mask = pos.locate_intersection(vcf['variants/POS'])[1]
    ref  = vcf['variants/REF'][mask][aims_pos_mask]
    alt = vcf['variants/ALT'][mask][aims_pos_mask]

    # filter geno array to set of aims
    geno_aims = geno.compress(aims_pos_mask, axis=0)

    totalgambscore = {}
    totalcoluscore = {}

    for aim in our_aims:

        gambscore = {}
        coluscore = {}

        # filter arrays 
        mask = our_aims == aim
        ref_ = ref[mask]
        alt_ = alt[mask]
        aimscolu_ = aimscolu[mask]
        aimsgamb_ = aimsgamb[mask]

        gn_aim = geno_aims.compress(mask, axis=0)

        # convert genotypes to nucleotides
        gn2nucleotide = {0:ref_[0],
                        1:alt_[0][0],
                         2:alt_[0][1],
                         3:alt_[0][2],
                        -1:float("nan")}
        gn = rnaseqpop.replace_with_dict2_generic(gn_aim, gn2nucleotide)

        # for each sample, get proportion of gambiae/coluzzii alleles
        # alleles that are different to both will be missed here
        for sample in metadata.treatment.unique():
            alleles = gn.take(subpops[sample], axis=1).flatten()
            
            # at each AIM, do we have gamb or colu alleles
            gamb = alleles[alleles != 'nan'] == aimsgamb_
            colu = alleles[alleles != 'nan'] == aimscolu_

            # get proportion of gamb v colu alleles at each locus
            gambscore[sample] = np.mean(gamb)
            coluscore[sample] = np.mean(colu)

        totalgambscore[aim] = dict(gambscore)
        totalcoluscore[aim] = dict(coluscore)

        gambscores = rnaseqpop.flip_dict(totalgambscore)
        coluscores = rnaseqpop.flip_dict(totalcoluscore)

        prop_gambiae = {}
        prop_colu = {}
        n_aims_per_sample = {}

        for sample in metadata.treatment.unique():

            prop_gambiae[sample] = np.nanmean(np.array(list(gambscores[sample].values())))
            all_gamb[sample].append(np.nanmean(np.array(list(gambscores[sample].values()))))
            prop_colu[sample] = np.nanmean(np.array(list(coluscores[sample].values())))
            all_colu[sample].append(np.nanmean(np.array(list(coluscores[sample].values()))))
            
            arr = np.array(list(gambscores[sample].values()))
            dim = arr.shape[0]
            n_aims_per_sample[sample] = dim-np.sum(np.isnan(arr))
            
    # store AIM fractions for each chromosome in outer dict 
    aims_chrom_gamb[contig] = dict(prop_gambiae)
    aims_chrom_colu[contig] = dict(prop_colu)
    n_aims_per_chrom[contig] = dict(n_aims_per_sample)

    # Store ancestry score per aim
    ancestryPerAim[contig] = pd.concat([pd.DataFrame(gambscores).add_suffix("_gamb"), pd.DataFrame(coluscores).add_suffix("_colu")], axis=1)
    ancestryPerAim[contig]['contig'] = contig

    # plot and store for each chromosome
    coludf = pd.DataFrame.from_dict(prop_colu, orient='index', columns=['AIM_fraction_coluzzii'])
    gambdf = pd.DataFrame.from_dict(prop_gambiae, orient='index', columns=['AIM_fraction_gambiae'])
    perchromdf = gambdf.merge(coludf, left_index=True, right_index=True)
    aimsperchromdf = pd.DataFrame.from_dict(n_aims_per_sample, orient='index', columns=['n_AIMs'])

    perchromdf.to_csv(f"results/variantAnalysis/ancestry/AIM_fraction_{contig}.tsv", sep="\t", index=True)
    rnaseqpop.plot_aims(perchromdf, aimsperchromdf, species1="coluzzii", species2="gambiae", figtitle=f"AIM_fraction_{contig}", total=False)


aims_chrom_gamb = rnaseqpop.flip_dict(aims_chrom_gamb)
aims_chrom_colu = rnaseqpop.flip_dict(aims_chrom_colu)
n_aims_per_chrom = rnaseqpop.flip_dict(n_aims_per_chrom)

# get ancestry per aim for later plotting on chromosome
ancestryPerAim = pd.concat(ancestryPerAim, axis=0)
ancestryPerAim.to_csv("results/variantAnalysis/ancestry/ancestryScorePerAim.tsv", sep="\t")

# get genome wide average AIM fractions
for k in all_gamb:
    all_gamb[k] = np.nanmean(all_gamb[k])
    all_colu[k] = np.nanmean(all_colu[k])

df1 = pd.DataFrame.from_dict(all_gamb, orient='index',columns=['AIM_fraction_gambiae'])
df2 = pd.DataFrame.from_dict(all_colu, orient='index', columns=['AIM_fraction_coluzzii'])
n_aimsdf = pd.DataFrame.from_dict(n_aims_per_chrom)
n_aimsdf.to_csv(f"results/variantAnalysis/ancestry/n_AIMS_per_chrom.tsv", sep="\t", index=True)

df = df1.merge(df2, left_index=True, right_index=True)
df.to_csv(f"results/variantAnalysis/ancestry/AIMs_summary.tsv", sep="\t", index=True)

rnaseqpop.plot_aims(df, n_aimsdf, species1="coluzzii", species2="gambiae", figtitle="AIM_fraction_whole_genome", total=True)




########### The same but for arabiensis v gambiae/coluzzii 
# script should be modularised but no time atm, isnt necessary

if metadata['species'].isin(['arabiensis']).any():

    aims = zarr.open(snakemake.input['aims_zarr_arab'], mode='r')

    aims_chrom_gamb = {}
    aims_chrom_arab = {}
    all_gamb = defaultdict(list)
    all_arab = defaultdict(list)
    n_aims_per_chrom = {}

    for contig in contigs:

        # read in and filter data
        path = f"results/variantAnalysis/vcfs/{dataset}.{contig}.vcf.gz"
        vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops = rnaseqpop.readAndFilterVcf(path=path,
                                                                contig=contig,
                                                                samples=metadata,
                                                                numbers=numbers,
                                                                qualflt=qualflt,
                                                                missingfltprop=missingprop)
        aimspos = aims[contig]['POS'][:]

        # get intersection of aims and our SNPs
        aims_pos_mask, aims_mask_2 = pos.locate_intersection(aimspos)
        our_aims = pos[aims_pos_mask]
        print(f"\n In the data, across all samples there are {our_aims.shape[0]} gamb-arab Ancestry Informative markers on Chromosome {contig}")

        # get gamb and colu alleles, and subset to aims that we have in the rna-seq data 
        aimsgamb = aims[contig]['gambcolu_allele'][:][aims_mask_2]
        aimsarab = aims[contig]['arab_allele'][:][aims_mask_2]

        # get mask that was used in readAndFilterVcf()
        mask = pos.locate_intersection(vcf['variants/POS'])[1]
        ref  = vcf['variants/REF'][mask][aims_pos_mask]
        alt = vcf['variants/ALT'][mask][aims_pos_mask]

        # filter geno array to set of aims
        geno_aims = geno.compress(aims_pos_mask, axis=0)

        totalgambscore = {}
        totalarabscore = {}

        for aim in our_aims:

            gambscore = {}
            arabscore = {}

            # filter arrays 
            mask = our_aims == aim
            ref_ = ref[mask]
            alt_ = alt[mask]
            aimsarab_ = aimsarab[mask]
            aimsgamb_ = aimsgamb[mask]

            gn_aim = geno_aims.compress(mask, axis=0)

            #convert genotypes to nucleotides
            gn2nucleotide = {0:ref_[0],
                            1:alt_[0][0],
                            2:alt_[0][1],
                            3:alt_[0][2],
                            -1:float("nan")}
            gn = rnaseqpop.replace_with_dict2_generic(gn_aim, gn2nucleotide)

            # for each sample, get proportion of gambiae/arabiensis alleles
            # alleles that are different to both will be missed here
            for sample in metadata.treatment.unique():
                alleles = gn.take(subpops[sample], axis=1).flatten()
                
                # at each AIM, do we have gamb or arab alleles
                gamb = alleles[alleles != 'nan'] == aimsgamb_
                arab = alleles[alleles != 'nan'] == aimsarab_

                # get proportion of gamb v arab alleles at each locus
                gambscore[sample] = np.mean(gamb)
                arabscore[sample] = np.mean(arab)

            totalgambscore[aim] = dict(gambscore)
            totalarabscore[aim] = dict(arabscore)

            gambscores = rnaseqpop.flip_dict(totalgambscore)
            arabscores = rnaseqpop.flip_dict(totalarabscore)

            prop_gambiae = {}
            prop_arab = {}
            n_aims_per_sample = {}

            for sample in metadata.treatment.unique():

                prop_gambiae[sample] = np.nanmean(np.array(list(gambscores[sample].values())))
                all_gamb[sample].append(np.nanmean(np.array(list(gambscores[sample].values()))))
                prop_arab[sample] = np.nanmean(np.array(list(arabscores[sample].values())))
                all_arab[sample].append(np.nanmean(np.array(list(arabscores[sample].values()))))

                arr = np.array(list(gambscores[sample].values()))
                dim = arr.shape[0]
                n_aims_per_sample[sample] = dim-np.sum(np.isnan(arr))

        aims_chrom_gamb[contig] = dict(prop_gambiae)
        aims_chrom_arab[contig] = dict(prop_arab)
        n_aims_per_chrom[contig] = dict(n_aims_per_sample)

        # Store ancestry score per aim
        ancestryPerAim[contig] = pd.concat([pd.DataFrame(gambscores).add_suffix("_gamb"), pd.DataFrame(coluscores).add_suffix("_colu")], axis=1)
        ancestryPerAim[contig]['contig'] = contig

        # plot and store for each chromosome
        gambdf = pd.DataFrame.from_dict(prop_gambiae, orient='index', columns=['AIM_fraction_gambiae'])
        arabdf = pd.DataFrame.from_dict(prop_arab, orient='index', columns=['AIM_fraction_arabiensis'])
        perchromdf = gambdf.merge(arabdf, left_index=True, right_index=True)
        aimsperchromdf = pd.DataFrame.from_dict(n_aims_per_sample, orient='index', columns=['n_AIMs'])

        perchromdf.to_csv(f"results/variantAnalysis/ancestry/AIM_fraction_{contig}.arab.tsv", sep="\t", index=True)
        rnaseqpop.plot_aims(perchromdf, aimsperchromdf, species1="arabiensis", species2="gambiae", figtitle=f"AIM_fraction_arab_{contig}", total=False)

    aims_chrom_gamb = rnaseqpop.flip_dict(aims_chrom_gamb)
    aims_chrom_arab = rnaseqpop.flip_dict(aims_chrom_arab)
    n_aims_per_chrom = rnaseqpop.flip_dict(n_aims_per_chrom)

    # get ancestry per aim for later plotting on chromosome
    ancestryPerAim = pd.concat(ancestryPerAim, axis=0)
    ancestryPerAim.to_csv("results/variantAnalysis/ancestry/ancestryScorePerAim.tsv", sep="\t")

    # get genome wide average AIM fractions
    for k in all_gamb:
        all_gamb[k] = np.nanmean(all_gamb[k])
        all_arab[k] = np.nanmean(all_arab[k])

    df1 = pd.DataFrame.from_dict(all_gamb, orient='index', columns=['AIM_fraction_gambiae'])
    df2 = pd.DataFrame.from_dict(all_arab, orient='index', columns=['AIM_fraction_arabiensis'])
    n_aimsdf = pd.DataFrame.from_dict(n_aims_per_chrom)
    n_aimsdf.to_csv(f"results/variantAnalysis/ancestry/n_AIMS_per_chrom_arab.tsv", sep="\t", index=True)

    df = df1.merge(df2, left_index=True, right_index=True)
    df.to_csv(f"results/variantAnalysis/ancestry/AIMs_summary_arab.tsv", sep="\t", index=True)

    rnaseqpop.plot_aims(df, n_aimsdf, species1="arabiensis", species2="gambiae", figtitle="AIM_fraction_whole_genome", total=True)