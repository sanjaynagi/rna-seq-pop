#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

args = commandArgs(trailingOnly=TRUE)
.libPaths(c("/home/snagi/miniconda3/lib/R/library", "~/R/x86_64-conda_cos6-linux-gnu-library/3.6", "/home/sanj/R/x86_64-pc-linux-gnu-library/3.6", 
          "/usr/local/lib/R/site-library", "/usr/lib/R/site-library",  "/usr/lib/R/library"))

## Differential expression
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)

##### functions ######
vst_pca = function(counts, samples, colourvar, name="PCA_", st="", comparison=""){
 
  #make DESeq dataset
    dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = samples, 
                               design = ~ treatment)

  ###### estimate paramters and normalise 
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  vsd = varianceStabilizingTransformation(dds)
  vstcounts = assay(vsd)
  vstcounts = vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]
  
  #### write pca of samples to pdf
  conds=as.factor(as.character(samples[,colourvar]))
  cond_colours = brewer.pal(length(levels(conds)),"Set2")[conds]
  names(cond_colours)=samples[,colourvar]
  pdf(glue("results/{name}{st}{comparison}.pdf"))
  pca2=prcomp(t(vstcounts),center=TRUE)
  plot(pca2$x, col=cond_colours,  pch=19, cex=2, main=glue("{name}{st}{comparison} (VST)"))
  text(pca2$x, as.vector(colnames(counts)), pos=3, cex=0.4)
  dev.off()
  
  return(list(vstcounts, dds))
}

print("------------- Kallisto - DESeq2 - RNASeq Differential expression ---------")
#### read data ####
print("Reading metadata file...")
samples = fread("data/samples.tsv") %>% as.data.frame()

#Read data
df = list()
for (sample in samples$samples){
  df[[sample]]= fread(glue("results/quant/{sample}/abundance.tsv"), sep = "\t")
}

counts = data.frame('Geneid' = df[[1]]$target_id)

#get read counts for each gene and fill table
for (sample in samples$samples){
  reads = df[[sample]]$est_counts
  counts = cbind(counts, reads)
}

#rename columns
colnames(counts) = c("Geneid", samples$samples)

## aggregate to gene level
counts$Geneid = str_remove(counts$Geneid, "-RA|-RB|-RC|-RD|-RE|-RF|-RG|-RH|-RI|-RJ|-RK|-RL|-RM|-RN")
counts = counts %>% group_by(Geneid) %>% summarise_all(sum)

### count total reads counted
counts = counts %>% column_to_rownames('Geneid')
samples$total_reads_ag = colSums(counts)

#make cohort a factor with 'control' as reference
samples$cohort = as.factor(as.character(samples$cohort))
samples$cohort = relevel(samples$cohort, "control")

#### get totals and plot ####
tot_counts = apply(counts, 2, sum) %>% 
  as.data.frame() %>% 
  rownames_to_column('sample')

colnames(tot_counts) = c("sample", "total_reads")

print("Counting and plotting total reads per sample...")
pdf("results/quant/total_reads_counted.pdf")
ggplot(tot_counts, aes(x=sample, y=total_reads, fill=samples$treatment)) + 
  geom_bar(stat='identity') + 
  theme_light() +
  ggtitle("Total reads counted (mapped to Ag transcriptome (PEST))") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(plot.title = element_text(hjust = 0.5))
null = dev.off()
 
counts = counts %>% rownames_to_column('Geneid') %>% 
  mutate_if(is.numeric, round) %>% column_to_rownames('Geneid')


############ Plots with all data ########################
vstcounts = vst_pca(counts, samples, 'strain', "PCA")[[1]]
#calculate correlations between samples based on the count data 
correlations = cor(vstcounts)
#### 
pdf("results/heatmap_correlations.pdf")
pheatmap(correlations)
garbage = dev.off()


######### subset data and run DESeq for each combo we need, store in xlsx ########
results_list = list()
names_list = list()
# series of loops to get combinations of comparisons - eg for each species in the dataset, 
# for each strain in the dataset, for each insecticide in the dataset, do a differential expression 
# analysis of case v control with DESEq2 and make PCAs and volcano plots.
for (sp in unique(samples$species)){
  df_samples = samples %>% filter(species == sp)
  lab_samples = df_samples %>% filter(lab == 'TRUE') #keep lab samples separate
  df_samples = df_samples %>% filter(lab == 'FALSE') #remove lab samples
  #print(lab_samples)
  for (st in unique(df_samples$strain)){
    ai_list = list()
    ctrl_list = list()
    ai_list[[st]] = df_samples[df_samples$strain == st,]
    
    #if the strain has both case and control (i.e. exposed v unexposed)
    if (length(unique(ai_list[[st]]$cohort)) > 1){
      print(glue("Running PCA for all {st} samples"))
      #do a pca for each strain 
      subcounts = counts[colnames(counts) %in% ai_list[[st]]$samples]
      #perform PCA on the data at strain level
      vstcounts = vst_pca(subcounts, ai_list[[st]], 'treatment', "PCA_", st=st)[[1]]
    }
      
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
      
      # if for this insecticide/strain combo, there are both no lab samples, and no insecticide exposed samples, skip this iteration of the loop
      # important in case we have an unbalanced dataset - ie, not all strains have all insecticides
      if ((nrow(lab_samples) < 1) & (nrow(df_s) < 1)){
        next
      }
      
      #for each set in the list df2, do DE. for  example either ngousso v tiefora or tiefora_control v tiefora_exposed 
      for (comparison in names(df2)){
        
        write(glue("{comparison}\n"), file="data/DE.comparison.list", append=TRUE)

        #subset count data to be columns we need for this comparison
        subcounts = counts[colnames(counts) %in% df2[[comparison]]$samples]
        #perform pca for each comparison
        print(glue("Running DESeq2 differential expression analysis on {comparison}"))
        
        if ('TRUE' %in% df2[[comparison]]$lab){

          sus = unique(df2[[comparison]][df2[[comparison]]$lab == 'TRUE',]$treatment)
          res = unique(df2[[comparison]][df2[[comparison]]$lab == 'FALSE',]$treatment)
	
          #make treatment a factor with the 'susceptible' as reference
          df2[[comparison]]$treatment = as.factor(as.character(df2[[comparison]]$treatment))
          df2[[comparison]]$treatment = relevel(df2[[comparison]]$treatment, sus)
          dds = vst_pca(subcounts, df2[[comparison]], colourvar='treatment', "PCA_", comparison=comparison)[[2]]
  
        } else {

            sus = unique(df2[[comparison]][df2[[comparison]]$cohort == 'control',]$treatment)
            res = unique(df2[[comparison]][df2[[comparison]]$cohort == 'case',]$treatment)

            #make treatment a factor with the 'susceptible' as reference
            df2[[comparison]]$treatment = as.factor(as.character(df2[[comparison]]$treatment))
            df2[[comparison]]$treatment = relevel(df2[[comparison]]$treatment, sus)
            dds = vst_pca(subcounts, df2[[comparison]], colourvar='treatment', "PCA_", comparison=comparison)[[2]]
        }

        #set pvalue threshold and perform  DE with each comparison
        p_threshold=0.05
        cds = nbinomWaldTest(dds)
        results = results(cds) %>% as.data.frame() 
        results = results[order(results$padj),] #order by pvalue 
        results = results %>% rownames_to_column("Gene_stable_ID") %>% mutate("FC" = (2^log2FoldChange))
        
        ###absolute difference
        #### remove _ and digits from columns in order to sum them. then get difference of counts and join with DE results
        readdiff = subcounts 
        colnames(readdiff) = str_replace(colnames(readdiff), "_", "") # remove digits from column names
        colnames(readdiff) = str_replace(colnames(readdiff), "[:digit:]+", "") # remove digits from column names        
        readdiff = data.frame(t(rowsum(t(readdiff), group = colnames(readdiff), na.rm = T))) #transpose and get rowsums for each group
        readdiff$absolute_diff = readdiff[,res] - readdiff[,sus] #get difference
        readdiff = data.frame(readdiff) %>% rownames_to_column('Gene_stable_ID')
        results = unique(left_join(results, readdiff[,c('Gene_stable_ID','absolute_diff')]))
          
        #add in gene names
        gene_names = fread("data/gene_names.tsv", sep="\t")
        #join DE results with normal gene names
        results_list[[comparison]] = unique(left_join(results, gene_names))
        fwrite(results_list[[comparison]], glue("results/diff/{comparison}.csv")) #write to csv 
    
        #store names of comparisons for xlsx report sheets
        names_list[[comparison]] = comparison
        
        #volcano plot for each comparison, first filter to remove very lowly expressed genes 
        a=results_list[[comparison]] %>% filter(`baseMean` > 20)
        pdf(glue("results/diff/Volcano_plot_{comparison}.pdf"))
        print(ggplot(a, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() +
                geom_hline(yintercept = (-log10(0.05)), linetype='dashed') +
                geom_vline(xintercept = (log2(2)), linetype='dashed') + 
                geom_vline(xintercept = log2(0.5), linetype='dashed') + 
                geom_text_repel(data=subset(results_list[[comparison]], abs(log2FoldChange) > 2 & padj < 0.0000001), aes(label=Gene_name)) +
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
saveWorkbook(wb,file="results/diff/RNA-Seq_diff.xlsx", overwrite = TRUE)

