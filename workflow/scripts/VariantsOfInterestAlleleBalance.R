#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(data.table)
library(glue)
library(openxlsx)

load_metadata <- function(metadata_path) {
  # Check the file extension and load metadata accordingly
  if (tools::file_ext(metadata_path) == "xlsx") {
    metadata <- readxl::read_excel(metadata_path)
  } else if (tools::file_ext(metadata_path) == "tsv") {
    metadata <- data.table::fread(metadata_path, sep = "\t")
  } else if (tools::file_ext(metadata_path) == "csv") {
    metadata <- data.table::fread(metadata_path, sep = ",")
  } else {
    stop("Metadata file must be .xlsx, .tsv, or .csv")
  }
  return(metadata)
}

#### allele imbalance ####
metadata = load_metadata(snakemake@input[['metadata']])
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
  
  # now lets calculate population means and lower and upper CIs
  alleles_per_pop_list = list()
  for (pop in unique(metadata$treatment)){
    
    alleles_per_pop = alleles %>% filter(treatment == pop) 
    
    pop_prop = sum(alleles_per_pop[, ..base])/ sum(alleles_per_pop[, 'cov'])
    error = sqrt((pop_prop*(1-pop_prop))/nrow(alleles_per_pop))*1.96
    lower = pmax(pop_prop - error, 0)
    upper = pmin(pop_prop + error, 1)
    
    alleles_per_pop_list[[pop]] = alleles_per_pop %>% mutate(!!propstring := pop_prop, lowerCI = lower, upperCI = upper)
    
    # average across replicates
    if (!base2 %in% c("", NA)){
      pop_prop2 = sum(alleles_per_pop[, ..base2])/ sum(alleles_per_pop[, 'cov'])
      error2 = sqrt((pop_prop*(1-pop_prop))/nrow(alleles_per_pop))*1.96
      lower2 = pmax(pop_prop - error, 0)
      upper2 = pmin(pop_prop + error, 1)
      
      alleles_per_pop_list[[pop]] = alleles_per_pop_list[[pop]] %>% mutate(!!propstring2 := pop_prop2, lowerCI_2 = lower2, upperCI_2 = upper2)
    }
  }
  
  mean_alleles = rbindlist(alleles_per_pop_list, fill=TRUE)
  
  # average across replicates
  if (!base2 %in% c("", NA)){
    mean_alleles = mean_alleles %>% 
      group_by(chrom, pos, ref, mutation, treatment, !!sym(propstring), lowerCI, upperCI, !!sym(propstring2), lowerCI_2, upperCI_2) %>% 
      summarise_at(.vars = c("cov","A","C","G","T"), .funs = c(mean="mean"), na.rm = TRUE) %>% 
      select(chrom, pos, ref, mutation, treatment, cov_mean, A_mean, C_mean, G_mean, T_mean, lowerCI, upperCI, !!propstring, !!propstring2, lowerCI_2, upperCI_2)
  } else {
    mean_alleles = mean_alleles %>% 
      group_by(chrom, pos, ref, mutation, treatment, !!sym(propstring), lowerCI, upperCI) %>% 
      summarise_at(.vars = c("cov","A","C","G","T"), .funs = c(mean="mean"), na.rm = TRUE) %>% 
      select(chrom,pos, ref, mutation, treatment, cov_mean, A_mean, C_mean, G_mean, T_mean, lowerCI, upperCI, !!propstring)
  }
  
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
