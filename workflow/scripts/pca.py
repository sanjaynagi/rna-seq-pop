#!/usr/bin/env python3

"""
A script to perform Principal components analysis on genotype data from RNA-Seq
"""

import matplotlib
import sys
from tools import *

dataset = snakemake.params[0]
chroms = snakemake.params[1]
missing_prop = snakemake.params[2]

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
        sns.despine(ax=ax, offset=5)
        x = coords[:, pc1]
        y = coords[:, pc2]
        for pop in populations:
            flt = (sample_population == pop)
            ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
                    label=pop, markersize=6, mec='k', mew=.5)
        ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
        ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title,path, sample_population=None):
        if sample_population is None:
            sample_population = samples.population.values
        # plot coords for PCs 1 vs 2, 3 vs 4
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 2, 1)
        plot_pca_coords(coords, model, 0, 1, ax, sample_population)
        ax = fig.add_subplot(1, 2, 2)
        plot_pca_coords(coords, model, 2, 3, ax, sample_population)
        ax.legend(bbox_to_anchor=(1, 1))
        fig.suptitle(title, y=1.02)
        
        fig.savefig(path, bbox_inches='tight', dpi=300)

def pca(geno, chrom, dataset, populations, prune=True, scaler=None):
    if prune is True:
        geno = geno.to_n_alt()
        gn = ld_prune(geno, size=500, step=200,threshold=0.2)
    else:
        gn = geno.to_n_alt()
        
    coords1, model1 = allel.pca(gn, n_components=10, scaler=scaler)
    
    fig_pca(coords1, model1, f"PCA-{chrom}-{dataset}", f"results/variants/PCA-{chrom}-{dataset}", sample_population=populations)

for chrom in chroms:
    path = f"results/variants/annot.variants.{chrom}.vcf.gz"
    #function to read in vcfs and associated SNP data
    vcf, geno, acsubpops, pos, depth, snpeff, subpops, populations =  readAndFilterVcf(path=path, 
                                                            chrom=chrom, 
                                                            qualflt=30,
                                                            missingfltprop=missing_prop,
                                                            plot=False)

    pop_colours = get_colour_dict(populations, "viridis")

    print(f"\n Performing PCA on {dataset} chromosome {chrom}")
    pca(geno, chrom, dataset, populations, prune=True, scaler=None)

    ### variant density 
    plot_density(pos, 500000, f"Variant Density chromosome {chrom}", title=None, path="results/variants/{dataset}_SNPdensity_{chrom}.png")