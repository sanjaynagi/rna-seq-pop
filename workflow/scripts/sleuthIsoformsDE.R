#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

args = commandArgs(trailingOnly=TRUE)
.libPaths(c("/home/snagi/miniconda3/lib/R/library", "~/R/x86_64-conda_cos6-linux-gnu-library/3.6", "/home/sanj/R/x86_64-pc-linux-gnu-library/3.6", 
          "/usr/local/lib/R/site-library", "/usr/lib/R/site-library",  "/usr/lib/R/library"))

## Differential expression
library(sleuth)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)

print("------------- Kallisto - Sleuth - RNASeq isoform Differential expression ---------")
#### read data ####
print("Reading metadata file...")

samples = fread(snakemake@input[[1]], sep="\t") %>% 
  as.data.frame() %>% 
  rename('sample' = "samples")

#add path column for sleuth object
samples$path = paste0("results/quant/",samples$sample)

#read metadata and get contrasts
gene_names = fread(snakemake@input[[2]], sep="\t") %>% 
  rename("GeneID" = "Gene_stable_ID")

#contrasts
contrastsdf = fread(snakemake@input[[3]])
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
    rename("GeneID" = "target_id") %>% 
    mutate("FC" = (2^b))
  
  #join DE results with normal gene names
  results_list[[cont]] = unique(left_join(results, gene_names))
  fwrite(results_list[[cont]], glue("results/isoformdiff/{cont}.csv")) #write to csv 
  
  #store names of conts for xlsx report sheets
  names_list[[cont]] = cont
  
  #volcano plot for each comparison, first filter to remove very lowly expressed genes 
  pdf(glue("results/isoformdiff/Volcano_plot_{cont}.pdf"))
  print(ggplot(results_list[[cont]], aes(x=b, y=-log10(qval))) + geom_point() +
          geom_hline(yintercept = (-log10(0.05)), linetype='dashed') +
          geom_vline(xintercept = (log2(2)), linetype='dashed') + 
          geom_vline(xintercept = log2(0.5), linetype='dashed') + 
          geom_text_repel(data=subset(results_list[[cont]], abs(FC) > 2 & qval < 0.0000001), aes(label=Gene_name)) +
          theme_light())
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
saveWorkbook(wb,file=snakemake@output[[1]], overwrite = TRUE)
  




  










































