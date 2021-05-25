#!/usr/bin/env python3

"""
A python script to check the validity of input files and parameters to the snakemake workflow. 
It will check the following and raise assertion errors if they fail.
- Fastq table samples matches that in samples.tsv
- The contrasts in config.yaml exist in samples.tsv treatment column
- The genome reference fasta chromosome headers match the config yaml
- The gff file chromosomes match the config yaml
- The gene_names.tsv file has appropriate column names.
"""

import sys
sys.stderr = open(snakemake.log[0], "w")
import os
import pandas as pd
import numpy as np
import re
import allel

# Read in parameters 
metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")
ref = snakemake.input['ref']
gffpath = snakemake.params['gffpath']
chroms = snakemake.params['chroms']
autofastq = snakemake.params['fastq']
sweeps = snakemake.params['sweeps']
signalpath = "resources/signals.csv"


## Read fasta file and grep chromosome names, check if they match config
file = open(ref, "r")
refchroms = []
for line in file:
    if re.search(">", line):
        refchroms.append(line.split()[0].replace(">", ""))
assert (np.isin(chroms, refchroms)).all(), f"Not all chromosome names specified in config.yaml are found in the reference fasta file ({ref})"


# Check that the contrasts specified in config.yaml are all in the treatment column
comparisons = pd.DataFrame(snakemake.params['contrasts'], columns=['contrast'])
comparisons = comparisons.contrast.str.split("_", expand=True) # split contrasts into case v control 
comparisons = comparisons.values.ravel() # make 1d array 
assert (np.isin(comparisons, metadata['treatment']).all()), f"treatments specified in config.yaml contrasts do not exist in samples.tsv. {comparisons}"


# Check that samples in fastq.tsv match those in metadata.tsv
if autofastq is False:
    fastq = pd.read_csv(snakemake.params['table'], sep="\t")
    assert np.isin(fastq['samples'], metadata['samples']).all(), f"all samples specified in fastq table do not match those in samples.tsv, please check your fastq.tsv file"
else:
    # Check to see if fastq files match the metadata
    for sample in metadata['samples']:
        for n in [1,2]:
            fqpath = f"resources/reads/{sample}_{n}.fastq.gz"
            assert os.path.isfile(fqpath), f"all sample names in 'samples.tsv' do not match a .fastq.gz file in the {snakemake.workflow.basedir}/resources/reads/ directory, please rename or consider using the fastq table option"


# Check column names of gene_names.tsv
gene_names = pd.read_csv(snakemake.params['gene_names'], sep="\t")
colnames = ['Gene_stable_ID', 'Gene_description', 'Gene_name']
assert (gene_names.columns.isin(colnames)).all(), f"Column names of gene_names.tsv should be {colnames}"


# Check that the chromosomes are present in the gff file 
gff = allel.gff3_to_dataframe(gffpath)
assert np.isin(chroms, gff.seqid).all(), f"All provided chromosome names ({chroms}) are not present in the reference GFF file"

# Check for signals.csv
if sweeps:
    assert os.path.isfile(signalpath), f"Ag1000g sweeps is activated, but no signals file is present in {signalpath}"