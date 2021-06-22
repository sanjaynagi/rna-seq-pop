# RNA-Seq-Pop

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.11.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/sanjaynagi/rna-seq-pop.svg?branch=master)](https://travis-ci.org/snakemake-workflows/rna-seq-pop)

This snakemake workflow performs various analyses of illumina paired-end RNA-Sequencing data:

* Quality control of fastq reads with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* QC metrics integrated into one final QC report with [multiQC](https://multiqc.info/)
* Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) at the gene level ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) and transcript level ([Sleuth](https://github.com/pachterlab/sleuth))
* Variant calling with [freebayes](https://github.com/freebayes/freebayes), and an Fst and [Population branch statistic (PBS)](https://science.sciencemag.org/content/329/5987/75) analysis, both in windows and at the gene-level ([Scikit-allel](https://scikit-allel.readthedocs.io/en/stable/)).
* Various summary statistics are calculated (Wattersons Theta, Sequence Diversity, Dxy etc)
* Differential SNP testing with the R package [kissDE](https://bioconductor.org/packages/release/bioc/html/kissDE.html), which accounts for allele-specific expression.
* Gene Set Enrichment analyses and Venn diagrams.
* Allele counts at pre-specified loci of choice.
* *Anopheles gambiae s.l* - Analysis of Ancestry Informative Markers (AIMs) to determine relative ancestry of *An.gambiae/coluzzii/arabiensis*. 
* *Anopheles gambiae s.l* - Reports if DE genes are found underneath known selective sweep signals in the [Ag1000g](https://www.nature.com/articles/nature24995).
* *Anopheles gambiae s.l* - Determines Karyotype of chromosome 2 inversions using [compkaryo](https://academic.oup.com/g3journal/article/9/10/3249/6026680) - [GH](https://github.com/sanjaynagi/compkaryo)

The workflow is generalised, and will function with any trimmed Illumina paired-end RNA-sequencing. However, certain modules, such as the AIMs analysis, are only appropriate for specific species. These can be activated in the configuration file (config.yaml). The workflow works with pooled samples, or diploid or haploid individuals. 

The workflow is still in construction, and not yet ready for release. If you have any feedback on how the workflow may be improved, please get in touch, or feel free to fork the github repo and create a pull request for any additional features you would like to implement. If you are using the workflow and would like to give feedback or troubleshoot, consider joining the discord server [here](https://discord.gg/RaXjP8APCq)

## Authors

* Sanjay Curtis Nagi (@sanjaynagi) - sanjay.c.nagi@gmail.com

## Usage

If you use this workflow in a paper, don't forget to give credits to the author by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.
3. As the workflow contains submodules (compkaryo & mpileup2readcounts), `git clone --recursive` should be used to clone the repo. 
4. If you want to run the differential SNPs analysis, please compile mpileup2readcounts by navigating to its folder in workflow/scripts and: `g++ -std=c++11 -O3 mpileup2readcounts.cc -o mpileup2readcounts`

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust the example `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/rna-seq-ir.git` or `git remote add -f upstream https://github.com/snakemake-workflows/rna-seq-ir.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 7: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).

