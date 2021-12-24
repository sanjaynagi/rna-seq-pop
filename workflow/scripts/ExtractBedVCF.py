#!/usr/bin/env python3

import sys
sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import pandas as pd
import allel

chroms = snakemake.params['chroms']

for chrom in chroms:
    vcf = allel.read_vcf(f"results/variantAnalysis/vcfs/annot.missense.{chrom}.vcf")
        
    pos = vcf['variants/POS']
    pos1 = pos+1
    
    data = {'chrom':chrom, 
        'start':pos,
        'stop':pos1}

    bed = pd.DataFrame(data)
    bed.to_csv(f"resources/regions/missense.pos.{chrom}.bed", sep="\t", header=None, index=None)
