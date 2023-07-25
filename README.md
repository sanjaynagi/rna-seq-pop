<div align="center">

<h1 align="center">
  RNA-Seq-Pop
</h1>

[<img src="https://github.com/sanjaynagi/rna-seq-pop/blob/master/RNA-Seq-Pop-Logo.png?raw=True" width="400"/>](https://github.com/sanjaynagi/rna-seq-pop/blob/master/RNA-Seq-Pop-Logo.png?raw=True)   

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.11.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![C.I.](https://github.com/sanjaynagi/rna-seq-pop/workflows/rna-seq-pop/badge.svg)](https://github.com/sanjaynagi/rna-seq-pop/actions?query=workflow:"rna-seq-pop-paired-end")
[![GitHub release](https://img.shields.io/github/release/sanjaynagi/rna-seq-pop?include_prereleases=&sort=semver&color=blue)](https://github.com/sanjaynagi/rna-seq-pop/releases/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6078337.svg)](https://doi.org/10.5281/zenodo.6078337)
[![Twitter Follow](https://img.shields.io/twitter/follow/sanjay_c_nagi.svg?style=social)](https://twitter.com/sanjay_c_nagi)

</div>

**Documentation**: https://sanjaynagi.github.io/rna-seq-pop/    

**Walkthrough**: https://www.youtube.com/watch?v=5QQe7DLHO4M      

RNA-Seq-Pop is a computational pipeline to analyse Illumina RNA-Sequencing data of any organism. As well as performing standard transcriptomic analyses, such as differential expression, RNA-Seq-Pop also calls and analyses genetic polymorphisms, extracting population genomic signals. The workflow can perform the following analyses of illumina paired-end RNA-Sequencing data:

* Quality control of fastq reads with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), QC metrics integrated into a final report with [multiQC](https://multiqc.info/)
* Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) at the gene level ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and transcript level ([Sleuth](https://github.com/pachterlab/sleuth)), and gene set enrichment analyses.

* Variant calling with [freebayes](https://github.com/freebayes/freebayes)
* Fst and [Population branch statistic (PBS)](https://science.sciencemag.org/content/329/5987/75), both in windows across contigs and at the gene-level ([Scikit-allel](https://scikit-allel.readthedocs.io/en/stable/)).
* Allele frequencies at variants of interest (pre-specified loci of choice).
* Various summary statistics (Wattersons Theta, Sequence Diversity, Dxy etc)    
* Analysis of Ancestry Informative Markers (AIMs) to determine relative ancestry of *An.gambiae/coluzzii/arabiensis* (*Anopheles gambiae s.l*)
* Determines Karyotype of chromosomal inversions using [compkaryo](https://academic.oup.com/g3journal/article/9/10/3249/6026680) (*Anopheles gambiae s.l*)

The workflow is generalised, and will function with any Illumina single or paired-end RNA-sequencing. However, certain modules, such as the ancestry and karyotyping, are only appropriate for *An. gambiae s.l*. These can be activated in the configuration file (config.yaml). The workflow works with pooled samples, diploid, or haploid individuals. 

If you have any feedback on how the workflow may be improved, please do get in touch, or feel free to fork the github repo and create a pull request for any additional features you would like to implement. If you are using the workflow and would like to give feedback or troubleshoot, consider joining the discord server [here](https://discord.gg/RaXjP8APCq), otherwise email or log an issue on github. 

## Authors

* Sanjay Curtis Nagi (@sanjaynagi) - sanjay.c.nagi@gmail.com

## Citation

Sanjay C Nagi, Ambrose Oruni, David Weetman, Martin J Donnelly (2022). **RNA-Seq-Pop**: Exploiting the sequence in RNA-Seq - a Snakemake workflow reveals patterns of insecticide resistance in the malaria vector *Anopheles gambiae*. *Molecular Ecology Resources*; doi: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13759

If you use this workflow in a paper, please give credits to the author by citing the manuscript and its DOI - https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13759 

## Usage 

Please see the [documentation](https://sanjaynagi.github.io/rna-seq-pop/    
) for instructions on how to run and contribute to RNA-Seq-Pop. 

## Release notes

* 1.0.4 - config files have been slightly restructured to make them neater. FDR control in GSEA. exampleconfig now up to date.
* 1.0.3 - support for single-end reads, venn and custom heatmaps added, updated results folder structure. fastqs can be specified in sample metadata.  
* 1.0.2 - Changed Pi, theta calculations to be in windows across genome, and removed plotting from SummaryStats.py. Changed use of 'chrom' to 'contig'
* 1.0.1 - New feature to plot a heatmap of various gene families 
