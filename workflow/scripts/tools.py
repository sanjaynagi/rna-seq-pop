# tools.py
import allel
import zarr
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial, reduce
from collections import defaultdict
from adjustText import adjust_text

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


def readAndFilterVcf(path, chrom, samples, numbers, ploidy, qualflt=30, missingfltprop=0.6, plot=True, verbose=False):

    """
    This function reads a VCF file, and filters it to a given quality and missingness proportion
    """
    
    print(f"\n-------------- Reading VCF for chromosome {chrom} --------------")
    vcf = allel.read_vcf(path, 
                      numbers=numbers,
                     fields=['calldata/*', 'variants/*', 'samples'])
    
    #get sample names and indices
    samplenames = vcf['samples']
    ind = defaultdict(list)

    for s,names in enumerate(samplenames):
        idx = np.where(np.isin(samples['samples'],names))[0][0]
        t = samples.treatment[idx]
        ind[t].append(s)
        subpops = dict(ind)

    if verbose: print(subpops, "\n")
    
    print(f"------- Filtering VCF at QUAL={qualflt} and missingness proportion of {missingfltprop} -------")
    #apply quality filters
    qual = vcf['variants/QUAL']
    passfilter = (qual >= qualflt)

    if ploidy == 1:
        geno = allel.HaplotypeArray(vcf['calldata/GT'].compress(passfilter, axis=0))
    else: 
        geno = allel.GenotypeArray(vcf['calldata/GT'].compress(passfilter, axis=0))

    pos = allel.SortedIndex(vcf['variants/POS'].compress(passfilter, axis=0))
    depth = vcf['variants/DP'].compress(passfilter, axis=0)
    print(f"After QUAL filter, {passfilter.sum()} SNPs retained out of {passfilter.shape[0]} for chromosome {chrom}")
      
    #missingness filters 
    ac = geno.count_alleles()
    snpcounts = ac.sum(axis=1)
    missingflt = snpcounts.max()*missingfltprop # must have at least 1/p alleles present
    missingness_flt = snpcounts >= missingflt
    geno = geno.compress(missingness_flt, axis=0)
    pos = pos[missingness_flt]
    print(f"After missingness filter, {missingness_flt.sum()} SNPs retained out of {missingness_flt.shape[0]} for chromosome {chrom}")
    
    #extract snpeff info and filter   
    snpeff = pd.DataFrame(vcf['variants/ANN'])[0].str.split("|", expand=True)[passfilter][missingness_flt]
    
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
            treatment = samples[samples['samples'] == pop]['treatment'].values[0]
            flt = (sample_population == pop)
            ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[treatment], 
                    label=treatment, markersize=14, mec='k', mew=.5)

        texts = [plt.text(x[i], y[i], pop, fontsize='small', ha='center', va='center') for i, pop in enumerate(sample_population)]
        adjust_text(texts)

        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title, path, samples, pop_colours,sample_population=None):
        if sample_population is None:
            sample_population = samples.samples.values
        # plot coords for PCs 1 vs 2, 3 vs 4
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 2, 1)
        plot_pca_coords(coords, model, 0, 1, ax, sample_population, samples, pop_colours)
        ax = fig.add_subplot(1, 2, 2)
        plot_pca_coords(coords, model, 2, 3, ax, sample_population, samples, pop_colours)
        legend_without_duplicate_labels(ax)
        fig.suptitle(title, y=1.02)
        
        fig.savefig(path, bbox_inches='tight', dpi=300)

def pca(geno, chrom, ploidy, dataset, populations, samples, pop_colours, prune=True, scaler=None):
    if prune is True:
        if ploidy > 1:
            geno = geno.to_n_alt()
        geno = ld_prune(geno, size=500, step=200,threshold=0.2)
    else:
        if ploidy > 1:
            geno = geno.to_n_alt()
        
    coords1, model1 = allel.pca(geno, n_components=10, scaler=scaler)

    fig_pca(coords1, model1, f"PCA {chrom} {dataset}", f"results/variants/plots/PCA-{chrom}-{dataset}", samples, pop_colours, sample_population=populations)






### AIMs plotting ###

def plot_aims(df, n_aims, species1="coluzzii", species2="gambiae", figtitle="AIM_fraction_overall", total=True):
    
    # if we are plotting the total genome wide AIM fraction, sum the number of obs for each chrom
    if total:
        n_aims = n_aims.sum(axis=0)
    else:
        n_aims = n_aims['n_AIMs']
    totalaims = ["n=" + str(t) for t in n_aims]

    #### Seaborn stacked barplots, overall AIM fraction ####
    sns.set_style("white")
    sns.set_context({"figure.figsize":(24,10)})

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
                   loc='right', bbox_to_anchor=(1.18, 0.5),
              ncol=1, fancybox=True, shadow=True, prop={'size':20})
    l.draw_frame(True)

    # label the bars with n_obs
    pos = range(len(totalaims))
    for tick,label in zip(pos, bottom_plot.get_xticklabels()):
        bottom_plot.text(pos[tick], 0.75, totalaims[tick], 
                horizontalalignment='center', size=15, color="k")

    # Label axes
    bottom_plot.set_ylabel("AIM fraction")
    bottom_plot.set_xlabel("Population")

    #Set fonts to consistent 16pt size
    for item in ([bottom_plot.xaxis.label, bottom_plot.yaxis.label] +
             bottom_plot.get_xticklabels() + bottom_plot.get_yticklabels()):
        item.set_fontsize(22)
    
    # add title and save figure
    plt.title(f"{figtitle}", fontsize=22, pad=20)
    plt.savefig(f"results/variants/AIMs/{figtitle}.png")
    plt.close()    
