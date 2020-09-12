#!/usr/bin/env Rscript
.libPaths(c("/home/snagi/miniconda3/lib/R/library", "~/R/x86_64-conda_cos6-linux-gnu-library/3.6", "/home/sanj/R/x86_64-pc-linux-gnu-library/3.6", 
            "/usr/local/lib/R/site-library", "/usr/lib/R/site-library",  "/usr/lib/R/library"))
library(tidyverse)
library(data.table)
library(glue)
library(openxlsx)

#### allele imbalance ####
metadata = fread("config/samples.tsv", sep="\t")
samples= metadata$samples
#read IR mutation data 
mutation_data = fread("resources/IRmutations.tsv", sep="\t")

#human readable treatment column

all_list = list()
mean_list = list()

for (m in mutation_data$Name){
  
  base = mutation_data[mutation_data$Name == m]$ALT
  propstring = glue("proportion{base}")
  base2 = mutation_data[mutation_data$Name == m]$ALT2 #### add in if second alt
  propstring2 = glue("proportion{base2}")

  #### load allele balance data ####
  allele_list = list()
  # read data for each sample and subset to what we want
  for (sample in samples){
    allele_list[[sample]] = fread(glue("results/allele_balance/counts/{sample}_{m}_allele_counts.tsv"))[,c(1:8)] #for each sample read data, and store first 8 columns 
    allele_list[[sample]]$sample = sample                                            #add sample column
    allele_list[[sample]]$treatment = metadata$treatment[samples == sample]         #add treatment column
    allele_list[[sample]]$mutation = m
    allele_list[[sample]]$gene = mutation_data[mutation_data$Name == m]$Gene
    
    allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring := (!!sym(base))/`cov`) #new column, proportion of As to total.  
    
    if (!base2 %in% c("", NA)){
    allele_list[[sample]] = allele_list[[sample]] %>% mutate(!!propstring2 := (!!sym(base2))/`cov`) #new column, proportion of As to total.  
    }
  }
  
  # We have 24 separate dataframes in a list, bind them together into one big dataframe
  alleles = rbindlist(allele_list, fill = TRUE)
  # average across replicates 
  mean_alleles = alleles %>% group_by(chrom, pos,ref, mutation, treatment) %>% summarise_at(.vars = c("cov","A","C","G","T", propstring)
                                                                                      , .funs = c(mean="mean"))
  #write to file, reorder mean_kdr_alleles
  fwrite(alleles, glue("results/allele_balance/csvs/{m}_allele_balance.csv"))
  fwrite(mean_alleles, glue("results/allele_balance/csvs/mean_{m}_allele_balance.csv"))
  
  all_list[[m]] = alleles
  mean_list[[m]] = mean_alleles
}

# might want to group in spreadhseets rather than one sheet each
#genelist = list()

#for (gene in unique(mutation_data$Gene)){
#  for (m in mutation_data$Name){
#  
#    if (mutation_data[mutation_data$Name == m]$Gene == gene){
#      genelist[[gene]][[m]] = all_list[[m]]
#    }
#  }
#}
#genelist


#### write to excel file on diff sheets #### 
results_list = all_list
sheets = unique(mutation_data$Name)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = TRUE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file="results/allele_balance/allele_balance.xlsx", overwrite = TRUE)


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
saveWorkbook(wb,file="results/allele_balance/mean_allele_balance.xlsx", overwrite = TRUE)

