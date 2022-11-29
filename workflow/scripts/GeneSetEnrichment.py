#!/usr/bin/env python3

"""
A script to perform GSEA on DE and FST data
"""
import sys
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import numpy as np

###### configuration - metadata and parameters ######
selection = snakemake.params['selection']
comparisons = pd.DataFrame(snakemake.params['DEcontrasts'], columns=['contrast'])

def load_go_descriptions():
    import urllib.request
    ids = []
    descriptions = []
    with urllib.request.urlopen("http://current.geneontology.org/ontology/go.obo") as url:
        for line in url:
            if line.startswith(b"id"):
                value = line.lstrip(b"id: ").rstrip(b"\n")
                if value.startswith(b"GO"):
                    ids.append(value)
                    descriptions.append(next(url, '').lstrip(b"name:").lstrip().rstrip(b"\n"))
    return(pd.DataFrame({'go_term': [go.decode('utf8') for go in ids], 'descriptions':[desc.decode('utf8') for desc in descriptions]}))

def go_hypergeometric(target_gene_list, gaf_df):
    
    # load gene annotation file 
    go_annotations = gaf_df[['go_term', 'descriptions']].rename(columns={'go_term':'annotation'}).drop_duplicates()
    gaf_df = gaf_df[['GeneID', 'go_term']].drop_duplicates()
    N = gaf_df.GeneID.unique().shape[0] #Total number of genes with some annotation 
    k = np.isin(gaf_df.loc[:, 'GeneID'].unique(), target_gene_list).sum() 
  
    hyper_geo = _hypergeometric(
      annotation_df=gaf_df, 
      column_name='go_term', 
      target_gene_list=target_gene_list,
      N=N,
      k=k)    
    hyper_geo = hyper_geo.merge(go_annotations, how='left')
    return(hyper_geo)

def _hypergeometric(annotation_df, column_name, target_gene_list, N, k):
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import fdrcorrection
    from tqdm import tqdm

    sig_list = []
    res_list = []
    unique_annots = annotation_df.loc[:, column_name].unique()
    for annot in tqdm(unique_annots):

        annot_genes = annotation_df.query("{col} == @annot".format(col=column_name))['GeneID']
        m = len(annot_genes)
        x = annot_genes.isin(target_gene_list).sum()
        # scipy hypergeom 
        res = hypergeom(M=N, 
                        n=m, 
                        N=k).sf(x-1)
        sig_list.append(annot)
        res_list.append(res)    

    hyper_geo = pd.DataFrame({'annotation': sig_list, 'pval':res_list})
    hypo, hyper_geo.loc[:, 'padj'] =  fdrcorrection(hyper_geo['pval'])    
    return(hyper_geo.sort_values(by='padj'))
    

# load gene annotation files and descriptions
gaffile = snakemake.input['gaf']
gaf_df = pd.read_csv(gaffile, sep="\t")
if len(gaf_df.columns) < 3:
    gaf_df = gaf_df.reset_index()
gaf_df = gaf_df.iloc[:, [1,4]]
gaf_df.columns = ['GeneID', 'go_term']

go_desc_df = load_go_descriptions()
gaf_df = gaf_df.merge(go_desc_df, how='left')

# load differential expression data and perform hypergeometric test
for comp in comparisons['contrast']:
    de_data = pd.read_csv(f"results/genediff/{comp}.csv")
    sig_genes = de_data.query("padj < 0.05 and FC > 2")['GeneID']
    
    gsea_df = go_hypergeometric(sig_genes, gaf_df)
    gsea_df.to_csv(f"results/gsea/genediff/{comp}.de.tsv", sep="\t")
    
    if selection:
        fst_data = pd.read_csv(f"results/variantAnalysis/selection/FstPerGene.tsv", sep="\t")
        fst_comp_df = fst_data.loc[:, ['GeneID', f'{comp}_zFst']].sort_values(by=f'{comp}_zFst', ascending=False).dropna()
        n_genes = fst_comp_df.shape[0]
        percentile_5 = int(n_genes* 0.05) #5th percentile
        fst_genes = fst_comp_df.iloc[:percentile_5].loc[:,'GeneID'].to_numpy()
        gsea_df = go_hypergeometric(fst_genes, gaf_df)
        gsea_df.to_csv(f"results/gsea/fst/{comp}.fst.tsv", sep="\t")
        