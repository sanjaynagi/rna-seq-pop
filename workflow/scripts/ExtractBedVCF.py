#!/usr/bin/env python3

import sys
sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import pandas as pd
import allel

contigs = snakemake.params['contigs']

for contig in contigs:
    vcf = allel.read_vcf(f"results/variantAnalysis/vcfs/annot.missense.{contig}.vcf")
        
    pos = vcf['variants/POS']
    pos1 = pos+1
    
    data = {'contig':contig, 
        'start':pos,
        'stop':pos1}

    bed = pd.DataFrame(data)
    bed.to_csv(f"resources/regions/missense.pos.{contig}.bed", sep="\t", header=None, index=None)
