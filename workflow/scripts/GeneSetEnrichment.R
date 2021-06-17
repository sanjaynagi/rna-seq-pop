#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


#' This script runs GSEA analysis on each contrast in the workflow
#' on gene lists, ranked using fold-change from DE data, Fst, PBS 
#' (optional) and differential SNPs analysis

# GSEA RNA-Seq Ag
library(GO.db)
library(gage)
library(fgsea)
library(data.table)
library(glue)
library(tidyverse)

## functions ##
runEnrich = function(rankedList, GeneSetList, outName){
  
  for (set in c("kegg", "GO")){
    # pathway enrichment with fgsea for both kegg and  GO terms
    fgseaRes = fgsea(pathways = GeneSetList[[set]], 
                     stats = na.omit(rankedList),
                     minSize  = 15,
                     maxSize  = 500,
                     nperm=10000)
    
    fgseaRes %>% arrange(padj) %>% fwrite(., file = glue("results/gsea/{outName}.{set}.tsv"), sep = "\t")
    
    topPathwaysUp = fgseaRes[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown = fgseaRes[ES < 0][head(order(pval), n=10), pathway]
    topPathways = c(topPathwaysUp, rev(topPathwaysDown))
    pdf(glue("results/gsea/{outName}.{set}.fgsea.pdf"))
    plotGseaTable(GeneSetList[[set]][topPathways], na.omit(rankedList), fgseaRes)
    null = dev.off()
    print(glue("{set} complete"))
  }
}


###### configuration - metadata and parameters ######
metadata = fread(snakemake@input[['metadata']]) %>% as.data.frame()
comparisons = data.frame("contrast" = snakemake@params[['DEcontrasts']])
gaffile = snakemake@input[['gaf']]
variantCalling = snakemake@params[['VariantCalling']]
pbs = snakemake@params[['pbs']]
diffsnps = snakemake@params[['diffsnps']]
pbscomps = snakemake@params[['pbscomps']]
replaceString = snakemake@params[['replaceStringKegg']]
speciesID = snakemake@params[['KeggSpeciesID']]

 
GeneSetList = list()
######### GO Terms #########
# read in gaf file and select to appropriate columns and rename
go = fread(gaffile, sep = "\t", skip = 1, header = FALSE) %>% 
  as_tibble() %>% dplyr::select(c(2,5)) %>% distinct() %>% dplyr::rename("GeneID" = V2, "GO.id" = V5)
# download GO terms and descriptions as we need descriptions
goterms = Term(GOTERM) %>% enframe() %>% dplyr::rename("GO.id" = "name", "description" = "value")
# join with our go terms
go = inner_join(go, goterms) %>% dplyr::select("description", "GeneID")
# turn into a large list of terms and all the genes with that term
GOlist = go %>% pivot_wider(names_from="description", values_from="GeneID") %>% transpose()
GeneSetList[["GO"]] = GOlist[[1]]

#### KEGG Pathways ####
kg.geneset = kegg.gsets(speciesID)
GeneSetList[['kegg']] = kg.geneset$kg.sets[kg.geneset$sigmet.idx]
# optional str_replace due to strange KEGG naming
if (!is.null(replaceString)){
  GeneSetList[['kegg']] = lapply(GeneSetList[['kegg']], str_remove, replaceString)
}


##### gene DE ####
for (comp in comparisons$contrast){
  print(glue("Running KEGG and GO enrichment analyses for {comp}"))
  # make ranked list using DE results, rank based on log2foldchange
  rank = fread(glue("results/genediff/{comp}.csv")) %>% 
    arrange("log2FoldChange") %>% 
    dplyr::select(c('GeneID', all_of("log2FoldChange"))) %>% 
    deframe()
  runEnrich(rankedList = rank, GeneSetList = GeneSetList, outName = glue("genediff/{comp}.DE"))
}


### Fst ####
if (variantCalling == TRUE){
  for (comp in comparisons$contrast){
  print(glue("Running KEGG and GO enrichment analyses for Fst {comp}"))
  # make ranked list using DE results, rank based on log2foldchange
  rank = fread("results/variants/FstPerGene.tsv") %>% 
    distinct() %>% 
    dplyr::select(GeneID, glue("{comp}_zFst")) %>% 
    deframe()

  rank = rank %>% sort(decreasing = TRUE)
  runEnrich(rankedList = rank, GeneSetList = GeneSetList, outName = glue("/fst/{comp}.FST"))
  }
}



### PBS ####
if (pbs == TRUE){
  for (comp in pbscomps){
  print(glue("Running KEGG and GO enrichment analyses for PBS {comp}"))
  # make ranked list using DE results, rank based on log2foldchange
  rank = fread("results/variants/PbsPerGene.tsv") %>% 
    distinct() %>% 
    drop_na() %>% 
    dplyr::select(GeneID, glue("{comp}PBS")) %>% 
    deframe()
  
  rank = rank %>% abs() %>% sort(decreasing = TRUE)
  
  runEnrich(rankedList = rank, GeneSetList = GeneSetList, outName = glue("/pbs/{comp}.PBS"))
  }
}


#### diff SNPs ####
if (diffsnps == TRUE){
  for (comp in comparisons$contrast){
  print(glue("Running KEGG and GO enrichment analyses for diffSNPs {comp}"))
  # make ranked list using DE results, rank based on log2foldchange
  rank = fread(glue("results/variants/diffsnps/{comp}.kissDE.tsv"), sep="\t") %>% 
    arrange("Deltaf/DeltaPSI") %>% 
    dplyr::select(c('GeneID', all_of("Deltaf/DeltaPSI"))) %>% 
    deframe()
  runEnrich(rankedList = rank, GeneSetList = GeneSetList, outName = glue("/diffsnps/{comp}.diffsnps"))
  }
}

sessionInfo()
