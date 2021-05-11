#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


#' This script uses the sleuth R package to perform differential isoform analysis
#' This will be performed for each 
#' Sleuth uses the bootstraps of kallisto to estimate uncertainty 

## Differential expression
library(tidyverse)
library(sleuth)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)
library(EnhancedVolcano)

print("------------- Kallisto - Sleuth - RNASeq isoform Differential expression ---------")
#### read data ####

samples = fread(snakemake@input[['samples']], sep="\t") %>% 
  as.data.frame() %>% 
  dplyr::rename('sample' = "samples")

#add path column for sleuth object
samples$path = paste0("results/quant/",samples$sample)

#read metadata and get contrasts
gene_names = fread(snakemake@input[['gene_names']], sep="\t") %>% 
  dplyr::rename("GeneID" = "Gene_stable_ID")

#contrasts
contrastsdf = data.frame("contrast" = snakemake@params[['DEcontrasts']])
contrasts = contrastsdf$contrast

######### subset data and run DESeq for each combo we need, store in xlsx ########
results_list = list()
names_list = list()
#loop through each chosen contrast and perform DE 
for (cont in contrasts){
  control = str_split(cont, "_")[[1]][1] #get first of string, which is control (or susceptible)
  case = str_split(cont, "_")[[1]][2] #get second of string, which is case (or resistant)
  controls = which(samples$treatment %in% control) #get indices of case/control
  cases = which(samples$treatment %in% case)
  
  ##subset to subcounts of our contrast
  subsamples = samples[c(controls, cases),]
  
  #make treatment a factor with the 'susceptible' as reference
  subsamples$treatment = as.factor(as.character(subsamples$treatment))
  subsamples$treatment = relevel(subsamples$treatment, control)

  print(glue("Running Sleuth differential isoform expression analysis on {cont}"))
  so = sleuth_prep(subsamples, extra_bootstrap_summary = TRUE, 
                    transformation_function = function(x) log2(x + 0.5))
  so = sleuth_fit(so, ~treatment, 'full')
  #fit reduced model necessary for isoform analysis
  so = sleuth_fit(so, ~1, 'reduced')
  #fit likelihood ratio test
  so = sleuth_wt(so, which_beta = paste0("treatment",case))
  
  results = sleuth_results(so, test =paste0("treatment",case), 'reduced:full', show_all = FALSE) %>% 
    dplyr::rename("GeneID" = "target_id") %>% 
    dplyr::mutate("FC" = (2^b))
  
  # Join DE results with normal gene names
  results = unique(left_join(results, gene_names))
  fwrite(results, glue("results/isoformdiff/{cont}.csv")) # write to csv 

  # Store names of contrasts and results tables for xlsx report sheets
  results_list[[cont]] = results
  names_list[[cont]] = cont
  
    # volcano plot for each comparison, using EnhancedVolcano. First make vector of labels which is AGAPs unless a gene name exists
  labels = results %>% mutate("Gene_name" = case_when(Gene_name == "" ~ GeneID,
                                     Gene_name == NA ~ GeneID,
                                     TRUE ~ Gene_name)) %>% select(Gene_name) %>% deframe()
  pdf(glue("results/isoformdiff/Volcano_plot_{cont}.pdf"))
  print(EnhancedVolcano(results_list[[cont]],
                  lab=labels,
                  x='b',
                  y='pval',
                  title = cont))
  garbage = dev.off()
  print(glue("{cont} complete!"))
}

#### write to excel file on diff sheets #### 
sheets = unlist(names_list)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = FALSE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file=snakemake@output[['xlsx']], overwrite = TRUE)
  
sessionInfo()
