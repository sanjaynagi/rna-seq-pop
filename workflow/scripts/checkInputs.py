#!/usr/bin/env python3

"""
A script to check the validity of input files and parameters to the workflow
"""

import sys
sys.stderr = open(snakemake.log[0], "w")
import os
from tools import *

metadata = pd.read_csv(snakemake.params['metadata'], sep="\t")
gffpath = snakemake.params['gffpath']
chroms = snakemake.params['chroms']

# Check to see if fastq files match the metadata
for sample in metadata.samples:
    for n in [1,2]:
        fqpath = f"resources/reads/{sample}_{n}.fastq.gz"
        assert os.path.isfile(fqpath), f"all sample names in 'samples.tsv' do not match a .fastq.gz file in the {snakemake.workflow.basedir}/resources/reads/ directory"

# Check that the chromosomes are present in the gff file 
gff = allel.gff3_to_dataframe(gffpath)
assert np.isin(chroms, gff.seqid).all(), f"All provided chromosome names ({chroms}) are not present in GFF file"

# Check column names of gene_names.tsv
gene_names = pd.read_csv(snakemake.params['gene_names'], sep="\t")
colnames = ['Gene_stable_ID', 'Gene_description', 'Gene_name']

assert (gene_names.columns == colnames).all(), f"Column names of gene_names.tsv should be {colnames}"
