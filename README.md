<h1 align="center">
  RNA-Seq-Pop
</h1>

[<img src="https://github.com/sanjaynagi/rna-seq-pop/blob/master/RNA-Seq-Pop-Logo.png?raw=True" width="200"/>](https://github.com/sanjaynagi/rna-seq-pop/blob/master/RNA-Seq-Pop-Logo.png?raw=True)   

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.11.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/sanjaynagi/rna-seq-pop.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-pop)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6078337.svg)](https://doi.org/10.5281/zenodo.6078337)

   
Sanjay C Nagi, Ambrose Oruni, David Weetman, Martin J Donnelly (2022). **RNA-Seq-Pop**: Exploiting the sequence in RNA-Seq - a Snakemake workflow reveals patterns of insecticide resistance in the malaria vector *Anopheles gambiae*. *bioRxiv* 2022.06.17.493894; doi: https://doi.org/10.1101/2022.06.17.493894 


This snakemake workflow can perform various analyses of illumina paired-end RNA-Sequencing data:

* Quality control of fastq reads with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), QC metrics integrated into one final QC report with [multiQC](https://multiqc.info/)
* Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) at the gene level ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and transcript level ([Sleuth](https://github.com/pachterlab/sleuth))
* Variant calling with [freebayes](https://github.com/freebayes/freebayes), and an Fst and [Population branch statistic (PBS)](https://science.sciencemag.org/content/329/5987/75) analysis, both in windows and at the gene-level ([Scikit-allel](https://scikit-allel.readthedocs.io/en/stable/)).
* Allele counts and frequencies at variants of interest (pre-specified loci of choice).
* Calculation of various summary statistics (Wattersons Theta, Sequence Diversity, Dxy etc)
* Gene Set Enrichment analyses.
  
* *Anopheles gambiae s.l* - Analysis of Ancestry Informative Markers (AIMs) to determine relative ancestry of *An.gambiae/coluzzii/arabiensis*. 
* *Anopheles gambiae s.l* - Reports if DE genes are found underneath known selective sweep signals in the [Ag1000g](https://www.nature.com/articles/nature24995).
* *Anopheles gambiae s.l* - Determines Karyotype of chromosomal inversions using [compkaryo](https://academic.oup.com/g3journal/article/9/10/3249/6026680) - [GH](https://github.com/sanjaynagi/compkaryo)

The workflow is generalised, and will function with any Illumina paired-end RNA-sequencing. However, certain modules, such as the ancestry and karyotyping, are only appropriate for *An. gambiae s.l*. These can be activated in the configuration file (config.yaml). The workflow works with pooled samples, diploid, or haploid individuals. 

If you have any feedback on how the workflow may be improved, please do get in touch, or feel free to fork the github repo and create a pull request for any additional features you would like to implement. If you are using the workflow and would like to give feedback or troubleshoot, consider joining the discord server [here](https://discord.gg/RaXjP8APCq), otherwise email or log an issue on github. 

## Authors

* Sanjay Curtis Nagi (@sanjaynagi) - sanjay.c.nagi@gmail.com

## Citation

Sanjay C Nagi, Ambrose Oruni, David Weetman, Martin J Donnelly (2022). **RNA-Seq-Pop**: Exploiting the sequence in RNA-Seq - a Snakemake workflow reveals patterns of insecticide resistance in the malaria vector *Anopheles gambiae*. *bioRxiv* 2022.06.17.493894; doi: https://doi.org/10.1101/2022.06.17.493894 

If you use this workflow in a paper, please give credits to the author by citing the preprint on bioRxiv and its DOI - https://doi.org/10.1101/2022.06.17.493894 

## Usage 

Please see the [wiki](https://github.com/sanjaynagi/rna-seq-pop/wiki) for instructions on how to run RNA-Seq-Pop. 

## Release notes

* 1.0.2 - Changed Pi, theta calculations to be in windows across genome, and removed plotting from SummaryStats.py. Changed use of 'chrom' to 'contig'
* 1.0.1 - New feature to plot a heatmap of various gene families 
