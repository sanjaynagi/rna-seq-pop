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
samples = fread(snakemake@input[[2]]) %>% as.data.frame() %>% rename('sample' = "samples")
#add path column for sleuth object
#samples$path = paste0("results/quant/",samples$sample) # do i need this?

######### subset data and run DESeq for each combo we need, store in xlsx ########
results_list = list()
names_list = list()
# series of loops to get combinations of comparisons - eg for each species in the dataset, 
# for each strain in the dataset, for each insecticide in the dataset, do a differential expression 
# analysis of case v control with DESEq2 and make PCAs and volcano plots.
for (sp in unique(samples$species)){
  df_samples = samples %>% filter(species == sp)
  lab_samples = df_samples %>% filter(lab == 'TRUE')   #keep lab samples separate
  df_samples = df_samples %>% filter(lab == 'FALSE')  #remove lab samples
  #print(lab_samples)
  for (st in unique(df_samples$strain)){
    ai_list = list()
    ctrl_list = list()
    ai_list[[st]] = df_samples[df_samples$strain == st,]
    
    ### subset dataframes to control and case 
    ctrl_list[[st]] = ai_list[[st]] %>% filter(insecticide == 'none')
    ai_list[[st]] = ai_list[[st]] %>% filter(insecticide != 'none')
    
    #in case there is no exposed samples...
    if (nrow(ai_list[[st]]) > 0){
      ais = unique(ai_list[[st]]$insecticide)
    } else {
      ais = unique(ctrl_list[[st]]$insecticide)
    }
    
    #for each insecticide or (none) do a comparison
    for (i in ais){
      df2 = list()
      df_s = ai_list[[st]] %>% filter(insecticide == i)
      # sometimes we will have lab controls, sometimes not, so the following if statements should allow us to do both comparisons automatically and separately
      # the results will be stored in the names list results_list. We can then use 
      # if for this insecticide/strain combo, there are both no lab samples, and no insecticide exposed samples, skip this iteration of the loop
      # important in case we have an unbalanced dataset - ie, not all strains have all insecticides
      if ((nrow(lab_samples) < 1) & (nrow(df_s) < 1)){
        next
      }
    
      #if we have exposed samples, store df of unexposed v exposed 
      if (nrow(df_s) > 1){
        print("Using exposed strain...")
        name = paste0(c(unique(ctrl_list[[st]]$treatment),unique(df_s$treatment)), collapse = '_') #name of comparison
        df2[[name]] = bind_rows(ctrl_list[[st]], df_s)
      }
      
      #if we have lab strain then store separately, lab v unexposed field 
      if (nrow(lab_samples) > 1){
        for (s in unique(lab_samples$strain)){ #in case of multiple lab strains
          print(glue("Using lab strain...{s}"))
          lab = lab_samples[lab_samples$strain == s,]
          name2 = paste0(c(unique(lab$treatment), unique(ctrl_list[[st]]$treatment)), collapse = '_') #name of comparison
          df2[[name2]] = bind_rows(lab, ctrl_list[[st]])
        }
      }
      
    
      #for each set in the list df2, do DE. for  example either ngousso v tiefora or tiefora_control v tiefora_exposed 
      for (comparison in names(df2)){


	if ('TRUE' %in% df2[[comparison]]$lab){
                sus = unique(df2[[comparison]][df2[[comparison]]$lab == 'TRUE',]$treatment)
                res = unique(df2[[comparison]][df2[[comparison]]$lab == 'FALSE',]$treatment)
        } else {
                sus = unique(df2[[comparison]][df2[[comparison]]$cohort == 'control',]$treatment)
                res = unique(df2[[comparison]][df2[[comparison]]$cohort == 'case',]$treatment)
        }

        #set susceptible to be the reference 
        df2[[comparison]]$treatment = relevel(factor(df2[[comparison]]$treatment), ref = sus)

        ## Run sleuth for each comparison
        print(glue("Running Sleuth differential isoform expression analysis on {comparison}"))
        so = sleuth_prep(df2[[comparison]], extra_bootstrap_summary = TRUE, 
                         transformation_function = function(x) log2(x + 0.5))
        so = sleuth_fit(so, ~treatment, 'full')
        #fit reduced model necessary for isoform analysis
        so = sleuth_fit(so, ~1, 'reduced')
        #fit likelihood ratio test
        so = sleuth_wt(so, which_beta = paste0("treatment",res))
        
        results = sleuth_results(so, test =paste0("treatment",res), 'reduced:full', show_all = FALSE) %>% 
          rename("Gene_stable_ID" = "target_id") %>% 
          mutate("FC" = (2^b))
        
        #add in gene names
        gene_names = fread("data/gene_names.tsv", sep="\t")
        #join DE results with normal gene names
        results_list[[comparison]] = unique(left_join(results, gene_names))
        fwrite(results_list[[comparison]], glue("results/isoformdiff/{comparison}.csv")) #write to csv 
        
        #store names of comparisons for xlsx report sheets
        names_list[[comparison]] = comparison
        
        #volcano plot for each comparison, first filter to remove very lowly expressed genes 
        pdf(glue("results/isoformdiff/Volcano_plot_{comparison}.pdf"))
        print(ggplot(results_list[[comparison]], aes(x=b, y=-log10(qval))) + geom_point() +
                geom_hline(yintercept = (-log10(0.05)), linetype='dashed') +
                geom_vline(xintercept = (log2(2)), linetype='dashed') + 
                geom_vline(xintercept = log2(0.5), linetype='dashed') + 
                geom_text_repel(data=subset(results_list[[comparison]], abs(FC) > 2 & qval < 0.0000001), aes(label=Gene_name)) +
                theme_light())
        garbage = dev.off()
        print(glue("{comparison} complete!"))
      }
    }
  }
}


#### write to excel file on diff sheets #### 
sheets = unlist(names_list)
wb <- createWorkbook("Workbook")

for (i in 1:length(sheets)){
  addWorksheet(wb, glue("{sheets[[i]]}"))
  writeData(wb, sheets[i], results_list[[i]], rowNames = FALSE, colNames = TRUE)
}
#### save workbook to disk once all worksheets and data have been added ####
saveWorkbook(wb,file="results/isoformdiff/RNA-Seq_isoformdiff.xlsx", overwrite = TRUE)
  




  










































