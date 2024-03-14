import allel
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def load_metadata(metadata_path):
    # load panel metadata
    if metadata_path.endswith('.xlsx'):
        metadata = pd.read_excel(metadata_path, engine='openpyxl')
    elif metadata_path.endswith('.tsv'):
        metadata = pd.read_csv(metadata_path, sep="\t")
    elif metadata_path.endswith('.csv'):
        metadata = pd.read_csv(metadata_path, sep=",")
    else:
        raise ValueError("Metadata file must be .xlsx or .csv")
    return metadata


def plotWindowed(statName, cohortText, cohortNoSpaceText, values, midpoints, prefix, contig, ylim, colour, save=True):

    """
    Saves to .tsv and plots windowed statistics
    """
    
    assert midpoints.shape == values.shape, f"arrays not same shape, midpoints shape - {midpoints.shape}, value shape - {values.shape}"
    
    if save:
        # store windowed statistics as .tsv 
        df = pd.DataFrame({'midpoint':midpoints, statName:values})
        df.to_csv(f"{prefix}/{cohortNoSpaceText}.{statName}.{contig}.tsv", sep="\t", index=False)

    xtick = np.arange(0, midpoints.max(), 2000000)
    ylim = np.max([ylim, values.max()])
    plt.figure(figsize=[16,8])
    sns.lineplot(x=midpoints, y=values, color=colour, linewidth=3)
    plt.xlim(0, midpoints.max()+1000)
    plt.ylim(0, ylim)
    plt.yticks(fontsize=14)
    plt.xticks(xtick, rotation=45, ha='right', fontsize=14)
    plt.ticklabel_format(style='plain', axis='x')
    plt.title(f"{statName} | {cohortText} | Chromosome {contig}", fontdict={'fontsize':20})
    if save: plt.savefig(f"{prefix}/{cohortNoSpaceText}.{statName}.{contig}.svg",format="svg", dpi=300)
    if save: plt.savefig(f"{prefix}/{cohortNoSpaceText}.{statName}.{contig}.pdf",format="pdf", dpi=300)

    

def plotRectangular(voiFreqTable, path, annot=True, xlab="Sample", ylab="Variant Of Interest", title=None, figsize=[10,10], cbar=True, vmax=None, rotate=True, cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True), dpi=300):
    plt.figure(figsize=figsize)
    #voiFreqTable = (voiFreqTable*100).astype(int)
    sns.heatmap(voiFreqTable, cmap=cmap, vmax=vmax, cbar=cbar,
                   linewidths=0.8,linecolor="white",annot=annot, fmt = '',annot_kws={"size": 18})
    if title != None: plt.title(title, pad=10)
    
    if rotate:
        plt.xticks(fontsize=16, rotation=45, ha='right',rotation_mode="anchor")#, labels=['Busia Parental', 'Busia Selected', 'Kisumu'], ticks=[0.5,1.5,2.5])
    else:
        plt.xticks(fontsize=13)
        
    plt.yticks(fontsize=14)
    plt.xlabel(xlab, fontdict={'fontsize':14}, labelpad=20)
    plt.ylabel(ylab, fontdict={'fontsize':14})
    plt.savefig(path, bbox_inches='tight', dpi=dpi)
    
def getAlleleFreqTable(muts, Path, var="sample", mean_=False, lowCov = 10):
    freqDict = {}
    covDict = {}
    cov_var = "cov" if mean_== False else "cov_mean"

    for mut in muts['Name']:
        df = pd.read_csv(Path.format(mut=mut))
        if mean_:
            df['gene'] = muts[muts.Name == mut]['Gene'].iloc[0]
        df['name'] = df['chrom'].astype(str) + ":"+ df['pos'].astype(str) + "  " + df['gene'].astype(str) + " | " + df['mutation'].astype(str)
        df['frequency'] = df.filter(like="proportion").sum(axis=1)
        freqDict[mut] = df[['name', var, 'frequency']]
        covDict[mut] = df[['name', var, cov_var]]

    voiData = pd.concat(freqDict)
    covData = pd.concat(covDict)
    voiFreqTable = voiData.pivot(index="name", columns=var).round(2).droplevel(0, axis=1)
    voiCovTable = covData.pivot(index="name", columns=var).round(2).droplevel(0, axis=1)

    #annotTable = (voiFreqTable*100).astype(int).astype(str) + "%" ## percentages
    annotTable = voiFreqTable.astype(str).apply(lambda x: x.str.strip("0")).applymap(addZeros)  ## decimals
    ## adding asterisks if low Cov 
    asteriskTable = voiCovTable.applymap(lambda x: "*" if x < lowCov else "")
    annotTable = annotTable + asteriskTable
    #annotTable = annotTable.applymap(lambda x: "" if x == "0%*" else x)
    #annotTable = annotTable.applymap(lambda x: "" if x == "0%" else x)
    return(voiFreqTable, annotTable)

def addZeros(x):
    if x == "1.":
        return("1")
    elif x == ".":
        return("")
    elif len(x) < 3:
        return(x + "0")
    else: 
        return(x)

def plotTwoRectangular(FreqTable1, annotdf1, FreqTable2, annotdf2, path, ylab="Variant Of Interest", annotFontsize=50, ylabfontsize=28 ,ytickfontsize=18, title1=None, title2=None, figsize=[20,10], ratio='auto', vmax=None, rotate=True, cmap=sns.cubehelix_palette(start=.5, rot=-.75, as_cmap=True), dpi=100):
    
    if ratio=='auto':
        ratio=[FreqTable1.shape[1],FreqTable2.shape[1]]
    else:
        ratio=[2,1]
        
    ## Load subplots
    fig, ax = plt.subplots(1,2, figsize=figsize, constrained_layout=True, gridspec_kw={'width_ratios': ratio})
    if title1 != None: ax[0].set_title(title1, fontsize=28)
    ## First heatmap
    sns.heatmap(ax=ax[0],
                data=FreqTable1, 
                cmap=cmap, 
                vmax=vmax, 
                cbar=False,
                linewidths=0.8,
                linecolor="white",
                annot=annotdf1,
                fmt = '', 
                annot_kws={"size": annotFontsize / np.sqrt(len(FreqTable1))})
    ax[0].set(xlabel="")
    plt.setp(ax[0].get_xticklabels(),fontsize=18, rotation=45, ha='right',rotation_mode="anchor")
    ax[0].set_ylabel(ylab, fontsize=ylabfontsize)
    plt.xlabel(None)     
    plt.setp(ax[0].get_yticklabels(),fontsize=ytickfontsize)
    
    ## Second heatmap
    sns.heatmap(ax=ax[1],
                data=FreqTable2, 
                cmap=cmap, 
                vmax=vmax, 
                cbar=True,
                linewidths=0.8,
                linecolor="white",
                annot=annotdf2, 
                yticklabels=False, 
                fmt = '', 
                annot_kws={"size": annotFontsize / np.sqrt(len(FreqTable2))})
    
    plt.setp(ax[1].get_xticklabels(),fontsize=18, rotation=45, ha='right',rotation_mode="anchor")
    ax[1].set(xlabel="", ylabel="")
    if title2 != None: plt.title(title2, fontsize=28)
    if path != None: plt.savefig(path, bbox_inches='tight', dpi=dpi)

    

# get indices of duplicate names
def list_duplicates(seq):
    """
    This function will find duplicate records and record their index 
    """
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return (tally)

def replace_with_dict2_generic(ar, dic, assume_all_present=False):
    """
    This function replaces values in a numpy array depending on a supplied dictionary
    """
    # Extract out keys and values
    k = np.array(list(dic.keys()))
    v = np.array(list(dic.values()))

    # Get argsort indices
    sidx = k.argsort()

    ks = k[sidx]
    vs = v[sidx]
    idx = np.searchsorted(ks,ar)

    if assume_all_present==0:
        idx[idx==len(vs)] = 0
        mask = ks[idx] == ar
        return np.where(mask, vs[idx], ar)
    else:
        return vs[idx]
        
def get_colour_dict(populations, palette="Set1"):
    """
    This function creates a colour palette for a provided list, returning a dict
    """

    cmap = plt.get_cmap(palette, len(np.unique(populations)))    # PiYG
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    pop_colours = {A: B for A, B in zip(np.unique(populations), colors)}
    return(pop_colours)


def get_numbers_dict(ploidy):
    numbers = {
        'samples':1,
        'variants/CHROM': 1,
        'variants/POS': 1,
        'variants/ID': 1,
        'variants/REF': 1,
        'variants/ALT': 'A',
        'variants/QUAL': 1,
        'variants/DP': 1,
        'variants/AN': 1,
        'variants/AC': 'A',
        'variants/AF': 'A',
        'variants/MQ': 1,
        'variants/ANN': 1,
        'calldata/DP': 1,
        'calldata/GT': ploidy,
        'calldata/GQ': 1,
        'calldata/HQ': 2,
        'calldata/AD': 'R',
        'calldata/MQ0': 1,
        'calldata/MQ': 1,
        }
    return(numbers)


def readAndFilterVcf(path, contig, samples, numbers, ploidy, qualflt=30, missingfltprop=0.6, verbose=False):

    """
    This function reads a VCF file, and filters it to a given quality and missingness proportion
    """
    
    print(f"\n-------------- Reading VCF for chromosome {contig} --------------")
    vcf = allel.read_vcf(path, 
                      numbers=numbers,
                     fields=['calldata/*', 'variants/*', 'samples'])
    
    #get sample names and indices
    samplenames = vcf['samples']
    ind = defaultdict(list)
    for s,names in enumerate(samplenames):
        idx = np.where(np.isin(samples['sampleID'],names))[0][0]
        t = samples.treatment[idx]
        ind[t].append(s)
        subpops = dict(ind)

    if verbose: print(subpops, "\n")
    
    print(f"------- Filtering VCF at QUAL={qualflt} and missingness proportion of {missingfltprop} -------")
    #apply quality filters
    qual = vcf['variants/QUAL']
    passfilter = qual >= qualflt
    print(f"QUAL filter will retain {passfilter.sum()} SNPs retained out of {passfilter.shape[0]} for chromosome {contig}")

       #missingness filters 
    ac = allel.GenotypeArray(vcf['calldata/GT']).count_alleles()
    snpcounts = ac.sum(axis=1)
    missingflt = snpcounts.max()*missingfltprop # must have at least 1/p alleles present
    missingness_flt = snpcounts >= missingflt
    print(f"Missingness filter will retain {missingness_flt.sum()} SNPs out of {missingness_flt.shape[0]} for chromosome {contig}")

    passfilter = np.logical_and(passfilter, missingness_flt)
    print(f"The combined filter will retain {passfilter.sum()} SNPs out of {passfilter.shape[0]} for chromosome {contig}")
    
    if ploidy == 1:
        geno = allel.HaplotypeArray(vcf['calldata/GT'].compress(passfilter, axis=0))
    else: 
        geno = allel.GenotypeArray(vcf['calldata/GT'].compress(passfilter, axis=0))

    pos = allel.SortedIndex(vcf['variants/POS'].compress(passfilter, axis=0))     
    depth = vcf['variants/DP'].compress(passfilter, axis=0)
    #extract snpeff info and filter   
    snpeff = pd.DataFrame(vcf['variants/ANN'])[0].str.split("|", expand=True)[passfilter]
    
    ac_subpops = geno.count_alleles_subpops(subpops)
    
    return(vcf, geno, ac_subpops, pos, depth, snpeff, subpops, samplenames)
    

def meanPBS(ac1, ac2, ac3, window_size, normalise):
    """
    This function calculate PBS on allele counts arrays and then takes the mean of all pbs values.
    """
    #pbs per variant
    pbs = allel.pbs(ac1, ac2, ac3, window_size=window_size, normed=normalise)
    #get average of all pbs values (will be per gene)
    meanpbs = np.nanmean(pbs)
    
    _, se, stats = allel.stats.misc.jackknife(pbs, statistic=lambda n: np.mean(n))
    
    return(meanpbs, se, pbs, stats)


def windowedDiversity(geno, pos, subpops, statistic='pi', window_size=20000):
    ### Estimate in windows separately
    pi_dict = {}
    for pop, idx in subpops.items():
        ac = geno.take(idx, axis=1).count_alleles()
        if statistic == 'pi':
            pi, windows, d, f = allel.windowed_diversity(pos, ac, size=window_size)
        elif statistic == 'theta':
            pi, windows, d, f = allel.windowed_watterson_theta(pos, ac, size=window_size)
        else:
            assert "statistic is neither pi or theta"
        pi_dict[pop] = pd.DataFrame(pi).rename(columns={0:statistic})

    pi = pd.concat(pi_dict).reset_index().rename(columns={'level_0':'treatment'})
    return(pi)

def diversity_ci_table(div_dict, statistic='pi'):
    import math
    div_stat = pd.concat(div_dict)
    stats = div_stat.groupby(['treatment'])[statistic].agg(['mean', 'count', 'std'])

    ci95_hi = []
    ci95_lo = []
    for i in stats.index:
        m, c, s = stats.loc[i]
        ci95_hi.append(m + 1.96*s/math.sqrt(c))
        ci95_lo.append(m - 1.96*s/math.sqrt(c))

    stats['ci95_hi'] = ci95_hi
    stats['ci95_lo'] = ci95_lo
    return(stats)













def flip_dict(dict_):
    """
    Inverts a nested dictionary
    """
    flipped = defaultdict(dict)
    for key, val in dict_.items():
        for subkey, subval in val.items():
            flipped[subkey][key] = subval
    return(flipped)

def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    """
    Performs LD pruning, originally from Alistair Miles' blog. 
    """
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn

def plot_ld(gn, title):
    m = allel.rogers_huff_r(gn) ** 2
    ax = allel.plot_pairwise_ld(m)
    ax.set_title(title)

def plot_density(pos, window_size, title, path):
    
    fig, ax = plt.subplots(figsize=(30, 10))
    sns.despine(ax=ax, offset=5)
    y, windows = allel.windowed_count(pos, size=window_size)
    x = np.mean(windows, axis=1)
    ax.plot(x, y/window_size)
    ax.set_ylabel('Density (bp$^{-1}$)')
    ax.set_xlabel('Position (bp)')
    if title:
        ax.set_title(title)
    fig.savefig(path)

def meanPBS(ac1, ac2, ac3, window_size, normalise):
    # pbs per variant
    pbs = allel.pbs(ac1, ac2, ac3, window_size=window_size, normed=normalise)
    # get average of all pbs values (will be per gene)
    meanpbs = np.nanmean(pbs)
    
    _, se, stats = allel.stats.misc.jackknife(pbs, statistic=lambda n: np.mean(n))
    
    return(meanpbs, se, pbs, stats)

def legend_without_duplicate_labels(ax, **kwargs):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), **kwargs)

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population, samples, pop_colours):
        sns.despine(ax=ax, offset=5)
        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in sample_population:
            treatment = samples[samples['sampleID'] == pop]['treatment'].values[0]
            flt = (sample_population == pop)
            ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[treatment], 
                    label=treatment, markersize=14, mec='k', mew=.5)

        #texts = [plt.text(x[i], y[i], pop, fontsize='small', ha='center', va='center') for i, pop in enumerate(sample_population)]
        #adjust_text(texts)

        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title, path, samples, pop_colours,sample_population=None):
        if sample_population is None:
            sample_population = samples['sampleID'].values
        # plot coords for PCs 1 vs 2, 3 vs 4
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 2, 1)
        plot_pca_coords(coords, model, 0, 1, ax, sample_population, samples, pop_colours)
        ax = fig.add_subplot(1, 2, 2)
        plot_pca_coords(coords, model, 2, 3, ax, sample_population, samples, pop_colours)
        legend_without_duplicate_labels(ax)
        fig.suptitle(title, y=1.02)
        
        fig.savefig(f"{path}.svg", format='svg', bbox_inches='tight', dpi=300)
        fig.savefig(f"{path}.pdf", format='pdf', bbox_inches='tight', dpi=300)


def pca(geno, contig, ploidy, dataset, populations, samples, pop_colours, prune=True, scaler=None):
    if prune is True:
        if ploidy > 1:
            geno = geno.to_n_alt()
        geno = ld_prune(geno, size=500, step=250,threshold=0.01)
    else:
        if ploidy > 1:
            geno = geno.to_n_alt()
        
    n_components = 10
    coords1, model1 = allel.pca(geno, n_components=n_components, scaler=scaler)
    coords = pd.DataFrame(coords1, columns=[f"PC{i}" for i in range(1, n_components)])

    fig_pca(coords1, model1, f"PCA {contig} {dataset}", f"results/variantAnalysis/pca/PCA-{contig}-{dataset}", samples, pop_colours, sample_population=populations)
    return coords, model1






### AIMs plotting ###

def plot_aims(df, n_aims, species1="coluzzii", species2="gambiae", figtitle="AIM_fraction_overall", total=True):
    
    # if we are plotting the total genome wide AIM fraction, sum the number of obs for each contig
    if total:
        n_aims = n_aims.sum(axis=0)
    else:
        n_aims = n_aims['n_AIMs']
    totalaims = ["n=" + str(t) for t in n_aims]

    #### Seaborn stacked barplots, overall AIM fraction ####
    sns.set_style("white")
    sns.set_context({"figure.figsize":(16,10)})

    total = df[f'AIM_fraction_{species1}'] + df[f'AIM_fraction_{species2}']

    sns.barplot(x = df.index, y=total, color="#e84c3d")
    bottom_plot = sns.barplot(x = df.index, y = df[f'AIM_fraction_{species2}'], color = "#3598db")

    box = bottom_plot.get_position()
    bottom_plot.set_position([box.x0, box.y0, box.width * 0.95, box.height])

    # Legend
    topbar = plt.Rectangle((0,0),1,1,fc="#e84c3d", edgecolor = 'none')
    bottombar = plt.Rectangle((0,0),1,1,fc='#3598db',  edgecolor = 'none')
    l = plt.legend([bottombar, topbar], 
                   [f'An. {species2}', f'An. {species1}'],
                   loc='best',# bbox_to_anchor=(1, 0.95),
              ncol=1, fancybox=True, shadow=True, prop={'size':20}, framealpha=0.8)
    l.draw_frame(True)

    # label the bars with n_obs
    pos = range(len(totalaims))
    for tick,label in zip(pos, bottom_plot.get_xticklabels()):
        bottom_plot.text(pos[tick], 0.75, totalaims[tick], 
                horizontalalignment='center', size=15, color="k")

    # Label axes
    bottom_plot.set_ylabel("AIM fraction")
    bottom_plot.set_xlabel("Population", labelpad=10)

    #Set fonts to consistent 16pt size
    for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
             bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
        item.set_fontsize(22)
    
    # add title and save figure
    plt.title(f"{figtitle}", fontsize=22, pad=20)
    plt.savefig(f"results/variantAnalysis/ancestry/{figtitle}.svg", dpi=300)
    plt.savefig(f"results/variantAnalysis/ancestry/{figtitle}.pdf", dpi=300)
    plt.close()    

def getSNPGffstats(gff, pos):
    """
    Calculates number of sites found that intersect with a GFF feature and the proportion % 
    """
    
    refBases = allel.SortedIndex(np.arange(1, gff['end'].max()+1))

    exons = gff.query("type == 'exon'")
    genes = gff.query("type == 'gene'")
    
    ## Get intergenic and intronic SNP numbers
    snpsInGenesBool = pos.locate_intersection_ranges(genes['start'], genes['end'])[0] # Locate genic snps
    intergenicSNPs = pos[~snpsInGenesBool]                 # get complement of genic snps - intergenic
    genicSNPs = pos[snpsInGenesBool]             
    snpsInExonsBool = genicSNPs.locate_intersection_ranges(exons['start'], exons['end'])[0]
    exonicSNPs = genicSNPs[snpsInExonsBool]
    intronicSNPs = genicSNPs[~snpsInExonsBool]

    ## Get intergenic and intronic SNP numbers in Reference genome
    snpsInGenesBool = refBases.locate_intersection_ranges(genes['start'], genes['end'])[0] # Locate genic snps
    RefIntergenicSNPs = refBases[~snpsInGenesBool]                 # get complement of genic snps - intergenic
    RefGenicSNPs = refBases[snpsInGenesBool]             
    snpsInExonsBool = RefGenicSNPs.locate_intersection_ranges(exons['start'], exons['end'])[0]
    RefExonicSNPs = RefGenicSNPs[snpsInExonsBool]
    RefIntronicSNPs = RefGenicSNPs[~snpsInExonsBool]
    
    RefDict = {}
    SamplesDict = {}
    for feature in ['chromosome', 'gene', 'exon', 'intron', 'three_prime_UTR', 'five_prime_UTR']:
                RefDict[feature] = getSNPsinType(gff, refBases, feature, verbose=False)
                SamplesDict[feature] = getSNPsinType(gff, pos, feature, verbose=False)
    
    RefDict['intergenic'] = [RefIntergenicSNPs.shape[0], refBases.shape[0], RefIntergenicSNPs.shape[0]/refBases.shape[0]]
    SamplesDict['intergenic'] = [intergenicSNPs.shape[0], pos.shape[0], intergenicSNPs.shape[0]/pos.shape[0]]
    RefDict['intron'] = [RefIntronicSNPs.shape[0], refBases.shape[0], RefIntronicSNPs.shape[0]/refBases.shape[0]]
    SamplesDict['intron'] = [intronicSNPs.shape[0], pos.shape[0], intronicSNPs.shape[0]/pos.shape[0]]
    
    ref = pd.DataFrame.from_dict(RefDict).T.rename(columns={2:'ReferenceGenomeProportion'})[['ReferenceGenomeProportion']]
    res = pd.DataFrame.from_dict(SamplesDict).T.rename(columns={0:'called', 1:'total', 2:'proportion'})
    res[['total', 'called']] = res[['total', 'called']].astype(int)
    ref = ref.merge(res, left_index=True, right_index=True)
    return(ref)



def getSNPsinType(gff, pos, feature, verbose=False):
    gf = gff.query("type == @feature")

    snps = pos.intersect_ranges(gf['start'], gf['end'])
    total = pos.shape[0]
    called = snps.shape[0]
    prop = called/total
    return(called, total, prop)
