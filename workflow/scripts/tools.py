#tools.py
import allel
import numpy as np
import pandas as pd 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial, reduce
from collections import defaultdict

#get indices of duplicate names
def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return (tally)

def get_colour_dict(populations, palette="Set1"):

    cmap = plt.get_cmap(palette, len(np.unique(populations)))    # PiYG
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        colors.append(matplotlib.colors.rgb2hex(rgb))

    pop_colours = {A: B for A, B in zip(np.unique(populations), colors)}
    return(pop_colours)

ploidy=10
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

def readAndFilterVcf(path, chrom, qualflt=30, missingfltprop=0.6, plot=True):
    
    print(f"\n-------------- Reading VCF for chromosome {chrom} --------------")
    vcf = allel.read_vcf(path, 
                      numbers=numbers,
                     fields=['calldata/*', 'variants/*', 'samples'])
    
    samples = pd.read_csv("config/samples.tsv", sep="\t")
    #get sample names and indices
    samplenames = vcf['samples']
    ind = defaultdict(list)

    for s,names in enumerate(samplenames):
        idx = np.where(np.isin(samples.samples,names))[0][0]
        t = samples.treatment[idx]
        ind[t].append(s)
        subpops = dict(ind)

    print(subpops, "\n")
    
    print(f"------- Filtering VCF at QUAL={qualflt} and missingness proportion of {missingfltprop} -------")
    #apply quality filters
    qual = vcf['variants/QUAL']
    passfilter = (qual >= qualflt)
    geno = allel.GenotypeArray(vcf['calldata/GT'].compress(passfilter, axis=0))
    pos = allel.SortedIndex(vcf['variants/POS'].compress(passfilter, axis=0))
    depth = vcf['variants/DP'].compress(passfilter, axis=0)
    print(f"After QUAL filter, {passfilter.sum()} SNPs retained out of {passfilter.shape[0]} for chromosome {chrom}")
    
    #plot histogram of QUAL values
    if plot is True:
        print("Plotting histogram of QUAL scores, with extreme outliers removed")
        fltqual = qual[~is_outlier(qual)]
        plt = sns.distplot(fltqual)
        plt.figure.savefig(f"../../results/variants/qc/qual.{chrom}.hist.png")
    
    #missingness filters 
    ac = geno.count_alleles()
    snpcounts = ac.sum(axis=1)
    missingflt = snpcounts.max()*missingfltprop # must have at least 2/3rds alleles present
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
    flipped = defaultdict(dict)
    for key, val in dict_.items():
        for subkey, subval in val.items():
            flipped[subkey][key] = subval
    return(flipped)

def ld_prune(gn, size, step, threshold=.1, n_iter=1):
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
    #pbs per variant
    pbs = allel.pbs(ac1, ac2, ac3, window_size=window_size, normed=normalise)
    #get average of all pbs values (will be per gene)
    meanpbs = np.nanmean(pbs)
    
    _, se, stats = allel.stats.misc.jackknife(pbs, statistic=lambda n: np.mean(n))
    
    return(meanpbs, se, pbs, stats)
