#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(data.table)
library(glue)
library(openxlsx)

#### allele imbalance ####
metadata = fread(snakemake@input[['metadata']], sep="\t")
samples = metadata$sampleID
# Read IR mutation data 
mutation_data = fread(snakemake@input[['mutations']], sep="\t")

all_list = list()
mean_list = list()

# Loop through each mutation, finding allele coverage at each position
for (m in mutation_data$Name){
  
  base = mutation_data[mutation_data$Name == m]$ALT
  propstring = glue("proportion{base}")
  base2 = mutation_data[mutation_data$Name == m]$ALT2 #### add in if second alt
  propstring2 = glue("proportion{base2}")

  #### load allele balance data ####
  allele_list = list()
  # read data for each sample and subset to what we want
  for (sample in samples){
    allele_list[[sample]] = fread(glue("results/variantAnalysis/variantsOfInterest/counts/{sample}_{m}_allele_counts.tsv"))[,c(1:8)] #for each sample read data, and store first 8 columns 
    allele_list[[sample]]$sample = sample                                            #add sample column
    allele_list[[sample]]$treatment = metadata$treatment[samples == sample]         #add treatment column
    allele_list[[sample]]$mutation = m
    allele_list[[sample]]$gene = mutation_data[mutation_data$Name == m]$Gene
    
    cover = allele_list[[sample]] %>% select(A,C,G,T) %>% rowSums()
    allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring := (!!sym(base))/cover) #new column, proportion of Alts to total.  
    
    if (!base2 %in% c("", NA)){
    allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring2 := (!!sym(base2))/cover) #new column, proportion of Alts to total.  
    }
  }
  
  # We have 24 separate dataframes in a list, bind them together into one big dataframe
  alleles = rbindlist(allele_list, fill = TRUE)
  # average across replicates 
  mean_alleles = alleles %>% group_by(chrom, pos,ref, mutation, treatment) %>% summarise_at(.vars = c("cov","A","C","G","T", propstring)
                                                                                      , .funs = c(mean="mean"))
  #write to file, reorder mean_kdr_alleles
  fwrite(alleles, glue("results/variantAnalysis/variantsOfInterest/csvs/{m}_alleleBalance.csv"))
  fwrite(mean_alleles, glue("results/variantAnalysis/variantsOfInterest/csvs/mean_{m}_alleleBalance.csv"))
  
  all_list[[m]] = alleles
  mean_list[[m]] = mean_alleles
}

#### write to excel file on diff sheets #### 
results_list = all_list
sheets = unique(mutation_data$Name)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = TRUE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file=snakemake@output[['alleleBalance']], overwrite = TRUE)


### mean balance ####
#### write to excel file on diff sheets #### 
results_list = mean_list
sheets = unique(mutation_data$Name)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = TRUE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file=snakemake@output[['mean_alleleBalance']], overwrite = TRUE)

sessionInfo()
