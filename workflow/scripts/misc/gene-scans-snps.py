
import allel
import numpy as np
import pandas as pd
import rnaseqpoptools as rnaseqpop
from operator import itemgetter
import sys

geneID=sys.argv[1]
gffpath=sys.argv[2]
vcfpath=sys.argv[3]
metadatapath=sys.argv[4]
ploidy=int(sys.argv[5])

metadata = rnaseqpop.load_metadata(metadatapath)
gff = allel.gff3_to_dataframe(gffpath, attributes=['Parent', 'ID'])
contig, start, end = gff.query("ID == @geneID")[['seqid', 'start', 'end']].values[0]

numbers = rnaseqpop.get_numbers_dict(ploidy)
vcf, geno, ac_subpops, pos, depth, snpeff, subpops, samplenames = rnaseqpop.readAndFilterVcf(vcfpath, contig, metadata, numbers, ploidy, qualflt=0, missingfltprop=0, verbose=False)
ann = allel.read_vcf(vcfpath,  fields=['ANN'], numbers={'ANN':1})


gene_bool = pos.locate_ranges([start], [end])
gene_ann = ann['variants/ANN'].compress(gene_bool, axis=0)
gene_pos = pos[gene_bool]
print(f"Locating {geneID}")
print(f"There are {len(gene_pos)} called in {geneID}")

df = pd.Series(gene_ann).str.split("|").to_frame()
indices = range(len(df[0][0]))
df = df[0].transform({f'{i+1}': itemgetter(i) for i in indices})
df = df.rename(columns={'1':'alt', '2':'type', '3':'mod', '4':'chr', '5':'chr2', '6':'region'})
print(f"Before the missense mask, there are {len(df)} rows")
missense_mask = df['type'] == 'missense_variant'
gene_df = df.query("@missense_mask").drop(columns=['chr2', '7', '8' ,'9', 'region', '12', '13', '14', '15', '16']).reset_index()
print(f"Manipulating annotation data")

assert missense_mask.sum() > 0, f"There are no missense variants in {geneID}"


for pop,ac in ac_subpops.items():

    gene_acs = ac.compress(gene_bool, axis=0).compress(missense_mask, axis=0)
    gene_freqs = gene_acs.to_frequencies()

    gene_freqs = pd.DataFrame(gene_freqs).rename(columns={1:f'{pop}_allele_2', 2:f'{pop}_allele_3', 3:f'{pop}_allele_4'}).drop(columns=[0,4,5])
    gene_df = pd.concat([gene_df, gene_freqs] , axis=1)
    
gene_df['max_af'] = gene_df.iloc[:,7:].max(axis=1)
gene_df.to_csv(f"{geneID}.aa.frequencies.tsv", sep="\t")
print("Writing frequencies to file...")