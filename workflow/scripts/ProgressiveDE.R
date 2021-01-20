#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(data.table)
library(tidyverse)
library(glue)

#read metadata and get contrasts
comps = snakemake@params['comps']

#comps = c(paste0(comps[[1]][1:3], collapse = "_"), paste0(comps[[1]][4:6], collapse = "_"))
####
for (cont in comps){
  res = str_split(cont, "_")[[1]][3] #get first of string, which is control 
  intermediate = str_split(cont, "_")[[1]][2] #get case 
  sus = str_split(cont, "_")[[1]][1]

  #### Gene diff ####
  one = fread(glue("results/genediff/{intermediate}_{res}.csv"))
  up1  = one %>% filter(FC > 1, padj < 0.05)
  down1 = one %>% filter(FC < 1, padj < 0.05)
  
  two = fread(glue("results/genediff/{sus}_{intermediate}.csv"))
  up2  = two %>% filter(FC > 1, padj < 0.05)
  down2 = two %>% filter(FC < 1, padj < 0.05)
  
  intersectdown = inner_join(down1, down2, by="GeneID", suffix=c("_field", "_lab"))
  intersectup = inner_join(up1, up2, by="GeneID", suffix=c("_field", "_lab"))
  
  fwrite(intersectup, glue("results/genediff/{sus}_{intermediate}_{res}.up.progressive.tsv"), sep="\t")
  fwrite(intersectdown, glue("results/genediff/{sus}_{intermediate}_{res}.down.progressive.tsv"), sep="\t")

  #### Isoforms #### 
  
  one = fread(glue("results/isoformdiff/{intermediate}_{res}.csv"))
  up1  = one %>% filter(FC > 1, qval < 0.05)
  down1 = one %>% filter(FC < 1, qval < 0.05)
  
  two = fread(glue("results/isoformdiff/{sus}_{intermediate}.csv"))
  up2  = two %>% filter(FC > 1, qval < 0.05)
  down2 = two %>% filter(FC < 1, qval < 0.05)
  
  intersectdown = inner_join(down1, down2, by="GeneID", suffix=c("_field", "_lab"))
  intersectup = inner_join(up1, up2, by="GeneID", suffix=c("_field", "_lab"))
  
  fwrite(intersectup, glue("results/isoformdiff/{sus}_{intermediate}_{res}.up.progressive.tsv"), sep="\t")
  fwrite(intersectdown, glue("results/isoformdiff/{sus}_{intermediate}_{res}.down.progressive.tsv"), sep="\t")
}


sessionInfo()
