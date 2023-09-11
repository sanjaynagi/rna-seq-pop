#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


#' Finds genes that are differentially expressed in the same direction
#' Across two DE comparisons. For example, a gene may be significantly 
#' upregulated in the control v intermediate phenotype comparison. 
#' And also in the intermediate phenotype vs full phenotype comparison.

library(data.table)
library(tidyverse)
library(glue)

# read metadata and get contrasts
comps = snakemake@params['comps'][[1]]
padj_threshold = snakemake@params['pval']
isoform_level = snakemake@params['isoform_level']
gene_level = snakemake@params['gene_level']
upper_fc = snakemake@params['fc']
lower_fc = 1/as.numeric(upper_fc) # if someone wants a FC threshold of 2, need to have lower threshold of 0.5.

for (cont in comps){
  sus = str_split(cont, "_")[[1]][1]           # get control/susceptible name (such as Kisumu)
  intermediate = str_split(cont, "_")[[1]][2]     # get intermediate
  res = str_split(cont, "_")[[1]][3]              # get last of string, which is resistant/case 

  if (gene_level == TRUE){

    #### Gene diff ####
    one = fread(glue("results/genediff/{intermediate}_{res}.csv"))
    up1  = one %>% filter(FC > upper_fc, padj < padj_threshold)
    down1 = one %>% filter(FC < lower_fc, padj < padj_threshold)
    
    two = fread(glue("results/genediff/{sus}_{intermediate}.csv"))
    up2  = two %>% filter(FC > upper_fc, padj < padj_threshold)
    down2 = two %>% filter(FC < lower_fc, padj < padj_threshold)
    
    intersectdown = inner_join(down1, down2, by="GeneID", suffix=c("_field", "_lab"))
    intersectup = inner_join(up1, up2, by="GeneID", suffix=c("_field", "_lab"))
    
    fwrite(intersectup, glue("results/genediff/{sus}_{intermediate}_{res}.up.progressive.tsv"), sep="\t")
    fwrite(intersectdown, glue("results/genediff/{sus}_{intermediate}_{res}.down.progressive.tsv"), sep="\t")
  }

  if (isoform_level == TRUE){
    #### Isoforms #### 
    one = fread(glue("results/isoformdiff/{intermediate}_{res}.csv"))
    up1  = one %>% filter(FC > upper_fc, qval < padj_threshold)
    down1 = one %>% filter(FC < lower_fc, qval < padj_threshold)
    
    two = fread(glue("results/isoformdiff/{sus}_{intermediate}.csv"))
    up2  = two %>% filter(FC > upper_fc, qval < padj_threshold)
    down2 = two %>% filter(FC < lower_fc, qval < padj_threshold)
    
    intersectdown = inner_join(down1, down2, by="TranscriptID", suffix=c("_field", "_lab"))
    intersectup = inner_join(up1, up2, by="TranscriptID", suffix=c("_field", "_lab"))
    
    fwrite(intersectup, glue("results/isoformdiff/{sus}_{intermediate}_{res}.up.progressive.tsv"), sep="\t")
    fwrite(intersectdown, glue("results/isoformdiff/{sus}_{intermediate}_{res}.down.progressive.tsv"), sep="\t")
  }

}


sessionInfo()