#!/usr/bin/env python3

"""
A script to determine the proportion of Anopheles gambiae/coluzzii alleles at Ancestry Informative Markers across the genome
"""

from tools import *

### AIMS ###
samples = pd.read_csv(snakemake.input['samples'], sep="\t")
samples = samples.sort_values(by='species').reset_index(drop=True)
chroms = snakemake.params['chroms']
qualflt = snakemake.params['qualflt']
missingprop = snakemake.params['missingprop']

#read AIMs
aims = zarr.open(snakemake.input['aims_zarr_gambcolu'], mode='r')

## initialize dicts
aims_chrom_gamb = {}
aims_chrom_colu = {}
all_gamb = defaultdict(list)
all_colu = defaultdict(list)

for chrom in chroms:

    # read in and filter data
    path = f"results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops =  readAndFilterVcf(path=path,
                                                               chrom=chrom,
                                                               samples=samples,
                                                               qualflt=qualflt,
                                                               missingfltprop=missingprop)
    aimspos = aims[chrom]['POS'][:]

    # get intersection of aims and our SNPs
    aims_pos_mask, aims_mask_2 = pos.locate_intersection(aimspos)
    our_aims = pos[aims_pos_mask]
    print(f"\n In the data, across all samples there are {our_aims.shape[0]} Ancestry Informative markers on Chromosome {chrom}")

    # get gamb and colu alleles, and subset to aims that we have in the rna-seq data 
    aimscolu = aims[chrom]['colu_allele'][:][aims_mask_2]
    aimsgamb = aims[chrom]['gamb_allele'][:][aims_mask_2]

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

        #convert genotypes to nucleotides
        gn2nucleotide = {0:ref_[0],
                        1:alt_[0][0],
                         2:alt_[0][1],
                         3:alt_[0][2],
                        -1:float("nan")}
        gn = replace_with_dict2_generic(gn_aim, gn2nucleotide)

        # for each sample, get proportion of gambiae/coluzzii alleles
        # alleles that are different to both will be missed here
        for sample in samples.treatment.unique():
            alleles = gn.take(subpops[sample], axis=1).flatten()
            
            # at each AIM, do we have gamb or colu alleles
            gamb = alleles[alleles != 'nan'] == aimsgamb_
            colu = alleles[alleles != 'nan'] == aimscolu_

            # get proportion of gamb v colu alleles at each locus
            gambscore[sample] = np.mean(gamb)
            coluscore[sample] = np.mean(colu)

        totalgambscore[aim] = dict(gambscore)
        totalcoluscore[aim] = dict(coluscore)

        gambscores = flip_dict(totalgambscore)
        coluscores = flip_dict(totalcoluscore)

        prop_gambiae = {}
        prop_colu = {}

        for sample in samples.treatment.unique():

            prop_gambiae[sample] = np.nanmean(np.array(list(gambscores[sample].values())))
            all_gamb[sample].append(np.nanmean(np.array(list(gambscores[sample].values()))))
            prop_colu[sample] = np.nanmean(np.array(list(coluscores[sample].values())))
            all_colu[sample].append(np.nanmean(np.array(list(coluscores[sample].values()))))

        aims_chrom_gamb[chrom] = dict(prop_gambiae)
        aims_chrom_colu[chrom] = dict(prop_colu)

aims_chrom_gamb = flip_dict(aims_chrom_gamb)
aims_chrom_colu = flip_dict(aims_chrom_colu)

#get AIM fractions per chromosome
pd.DataFrame.from_dict(aims_chrom_gamb).to_csv("results/variants/AIMs/AIMs_gambiae.tsv", sep="\t", index=True)
pd.DataFrame.from_dict(aims_chrom_colu).to_csv("results/variants/AIMs/AIMs_coluzzii.tsv", sep="\t", index=True)

# get genome wide average AIM fractions
for k in all_gamb:
    all_gamb[k] = np.nanmean(all_gamb[k])
    all_colu[k] = np.nanmean(all_colu[k])

df1 = pd.DataFrame.from_dict(all_gamb, orient='index',columns=['AIM_fraction_gambiae'])
df2 = pd.DataFrame.from_dict(all_colu, orient='index', columns=['AIM_fraction_coluzzii'])

df = df1.merge(df2, left_index=True, right_index=True)
df.to_csv(f"results/variants/AIMs/AIMs_summary.tsv", sep="\t", index=True)




#### Seaborn stacked barplots, overall AIM fraction ####
sns.set_style("white")
sns.set_context({"figure.figsize":(24,10)})

total = df['AIM_fraction_coluzzii'] + df['AIM_fraction_gambiae']

sns.barplot(x = df.index, y=total, color="Red")
bottom_plot = sns.barplot(x = df.index, y =df['AIM_fraction_gambiae'], color = "#0000A3")

topbar = plt.Rectangle((0,0),1,1,fc="red", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor = 'none')
l = plt.legend([bottombar, topbar], ['An. gambiae', 'An. coluzzii'], loc=1, ncol = 2, prop={'size':16})
l.draw_frame(True)

#Optional code - Make plot look nicer
bottom_plot.set_ylabel("AIM fraction")
bottom_plot.set_xlabel("Sample")

#Set fonts to consistent 16pt size
for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
             bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
    item.set_fontsize(16)

plt.savefig(f"results/variants/AIMs/AIM_fraction_overall.png")




########### The same but for arabiensis v gambiae/coluzzii 
# script should be modularised but no time atm, isnt necessary

if samples['species'].isin(['arabiensis'].any():

    aims = zarr.open(snakemake.input['aims_zarr_arab'], mode='r')

    aims_chrom_gamb = {}
    aims_chrom_arab = {}
    all_gamb = defaultdict(list)
    all_arab = defaultdict(list)

    for chrom in chroms:

        # read in and filter data
        path = f"results/variants/vcfs/annot.variants.{chrom}.vcf.gz"
        vcf, geno, acsubpops, pos, depth, snpeff, subpops, pops =  readAndFilterVcf(path=path,
                                                                chrom=chrom,
                                                                samples=samples,
                                                                qualflt=qualflt,
                                                                missingfltprop=missingprop)
        aimspos = aims[chrom]['POS'][:]

        # get intersection of aims and our SNPs
        aims_pos_mask, aims_mask_2 = pos.locate_intersection(aimspos)
        our_aims = pos[aims_pos_mask]
        print(f"\n In the data, across all samples there are {our_aims.shape[0]} gamb-arab Ancestry Informative markers on Chromosome {chrom}")

        # get gamb and colu alleles, and subset to aims that we have in the rna-seq data 
        aimsgamb = aims[chrom]['gambcolu_allele'][:][aims_mask_2]
        aimsarab = aims[chrom]['arab_allele'][:][aims_mask_2]

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
            gn = replace_with_dict2_generic(gn_aim, gn2nucleotide)

            # for each sample, get proportion of gambiae/arabiensis alleles
            # alleles that are different to both will be missed here
            for sample in samples.treatment.unique():
                alleles = gn.take(subpops[sample], axis=1).flatten()
                
                # at each AIM, do we have gamb or arab alleles
                gamb = alleles[alleles != 'nan'] == aimsgamb_
                arab = alleles[alleles != 'nan'] == aimsarab_

                # get proportion of gamb v arab alleles at each locus
                gambscore[sample] = np.mean(gamb)
                arabscore[sample] = np.mean(arab)

            totalgambscore[aim] = dict(gambscore)
            totalarabscore[aim] = dict(arabscore)

            gambscores = flip_dict(totalgambscore)
            arabscores = flip_dict(totalarabscore)

            prop_gambiae = {}
            prop_arab = {}

            for sample in samples.treatment.unique():

                prop_gambiae[sample] = np.nanmean(np.array(list(gambscores[sample].values())))
                all_gamb[sample].append(np.nanmean(np.array(list(gambscores[sample].values()))))
                prop_arab[sample] = np.nanmean(np.array(list(arabscores[sample].values())))
                all_arab[sample].append(np.nanmean(np.array(list(arabscores[sample].values()))))

            aims_chrom_gamb[chrom] = dict(prop_gambiae)
            aims_chrom_arab[chrom] = dict(prop_arab)

    aims_chrom_gamb = flip_dict(aims_chrom_gamb)
    aims_chrom_arab = flip_dict(aims_chrom_arab)

    #get AIM fractions per chromosome
    pd.DataFrame.from_dict(aims_chrom_gamb).to_csv("results/variants/AIMs/AIMs_gamb_arab.tsv", sep="\t", index=True)
    pd.DataFrame.from_dict(aims_chrom_arab).to_csv("results/variants/AIMs/AIMs_arab_gamb.tsv", sep="\t", index=True)

    # get genome wide average AIM fractions
    for k in all_gamb:
        all_gamb[k] = np.nanmean(all_gamb[k])
        all_arab[k] = np.nanmean(all_arab[k])

    df1 = pd.DataFrame.from_dict(all_gamb, orient='index', columns=['AIM_fraction_gambiae'])
    df2 = pd.DataFrame.from_dict(all_arab, orient='index', columns=['AIM_fraction_arabiensis'])

    df = df1.merge(df2, left_index=True, right_index=True)
    df.to_csv(f"results/variants/AIMs/AIMs_summary_arab.tsv", sep="\t", index=True)


    #### Seaborn stacked barplots, overall AIM fraction ####
    sns.set_style("white")
    sns.set_context({"figure.figsize":(24,10)})

    total = df['AIM_fraction_arabiensis'] + df['AIM_fraction_gambiae']

    sns.barplot(x = df.index, y=total, color="Red")
    bottom_plot = sns.barplot(x = df.index, y =df['AIM_fraction_gambiae'], color = "#0000A3")

    topbar = plt.Rectangle((0,0),1,1,fc="red", edgecolor = 'none')
    bottombar = plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor = 'none')
    l = plt.legend([bottombar, topbar], ['An. gamb/colu', 'An. arabiensis'], loc=1, ncol = 2, prop={'size':16})
    l.draw_frame(True)

    #Optional code - Make plot look nicer
    bottom_plot.set_ylabel("AIM fraction")
    bottom_plot.set_xlabel("Sample")

    #Set fonts to consistent 16pt size
    for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
                bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
        item.set_fontsize(16)

    plt.savefig(f"results/variants/AIMs/AIM_fraction_Arab_overall.png")
