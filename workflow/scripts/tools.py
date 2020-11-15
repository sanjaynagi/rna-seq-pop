#tools.py
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

#get indices of duplicate names
def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return (tally)

def replace_with_dict2_generic(ar, dic, assume_all_present=False):
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

def readAndFilterVcf(path, chrom, samples, qualflt=30, missingfltprop=0.6, plot=True, verbose=False):
    
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

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population, samples, pop_colours):
        sns.despine(ax=ax, offset=5)
        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in sample_population:
            treatment = samples[samples['samples'] == pop]['treatment'].values[0]
            flt = (sample_population == pop)
            ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[treatment], 
                    label=treatment, markersize=6, mec='k', mew=.5)
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
        ax.legend(bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        
        fig.savefig(path, bbox_inches='tight', dpi=300)

def pca(geno, chrom, dataset, populations, samples, pop_colours, prune=True, scaler=None):
    if prune is True:
        geno = geno.to_n_alt()
        gn = ld_prune(geno, size=500, step=200,threshold=0.2)
    else:
        gn = geno.to_n_alt()
        
    coords1, model1 = allel.pca(gn, n_components=10, scaler=scaler)

    fig_pca(coords1, model1, f"PCA-{chrom}-{dataset}", f"results/variants/plots/PCA-{chrom}-{dataset}", samples, pop_colours, sample_population=populations)



#### Garuds G12 ####
# Instead of haps.discrete_frequencies() in allel.garuds_h() which does not allow for differences between multi-locus genotypes

def cluster_G12(gnalt, cut_height=0.1, metric='euclidean'):
    
    #cluster the genotypes in the window
    dist = scipy.spatial.distance.pdist(gnalt.T, metric=metric)
    if metric in {'hamming', 'jaccard'}:
        # convert distance to number of SNPs, easier to interpret
        dist *= gnalt.shape[0]

    Z = scipy.cluster.hierarchy.linkage(dist, method='single')

    cut = scipy.cluster.hierarchy.cut_tree(Z, height=cut_height)[:, 0]
    cluster_sizes = np.bincount(cut)
    clusters = [np.nonzero(cut == i)[0] for i in range(cut.max() + 1)]
    
    #get freq of clusters and sort by largest freq
    f = cluster_sizes/gnalt.shape[1]
    f = np.sort(f)[::-1]
    
    #calculate g12
    g12 = np.sum(f[:2])**2 + np.sum(f[2:]**2)
    
    return(g12)

def garuds_G12(gnalt, pos, cut_height=10, window_size=1000, save=False, name=None, metric='euclidean'):
    
    g12 = allel.moving_statistic(gnalt, cluster_G12, size=window_size, metric=metric, cut_height=cut_height)
    midpoint = allel.moving_statistic(pos, np.median, size=window_size)
    
    plt.figure(figsize=[20,10])
    sns.scatterplot(midpoint, g12)
    plt.title("G12 clustered BFgam")
    plt.show()
    if save: plt.savefig(f"{name}.png")
