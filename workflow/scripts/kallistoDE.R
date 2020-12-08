#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Differential expression
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(data.table)
library(ggrepel)
library(openxlsx)
library(glue)
library(RColorBrewer)
library(EnhancedVolcano)

#read metadata and get contrasts
samples = fread(snakemake@input[[1]], sep="\t") %>% as.data.frame()
gene_names = fread(snakemake@input[[2]], sep="\t") %>% 
  rename("GeneID" = "Gene_stable_ID")

contrastsdf = fread(snakemake@input[[3]])
contrasts = contrastsdf$contrast

##### functions ######
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}

vst_pca = function(counts, samples, colourvar, name="PCA_", st="", comparison=""){
 
  #make DESeq dataset
  dds = DESeqDataSetFromMatrix(countData = counts, 
                               colData = samples, 
                               design = ~ treatment)

  ###### estimate paramters and normalise 
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  vsd = varianceStabilizingTransformation(dds)
  normcounts = counts(dds, normalized=TRUE)
  vstcounts = assay(vsd)
  vstcounts = vstcounts[order(apply(vstcounts,1,sum),decreasing =TRUE),]
  
  #### write pca of samples to pdf
  conds=as.factor(as.character(samples[,colourvar]))
  cond_colours = brewer.pal(length(levels(conds)),"Set2")[conds]
  names(cond_colours)=samples[,colourvar]
  pdf(glue("results/plots/{name}{st}{comparison}.pdf"))
  pca2=prcomp(t(vstcounts),center=TRUE)
  plot(pca2$x, col=cond_colours,  pch=19, cex=2, main=glue("{name}{st}{comparison} (VST)"))
  text(pca2$x, as.vector(colnames(counts)), pos=3, cex=0.4)
  dev.off()
  
  return(list(vstcounts, dds, normcounts))
}

#### main ####
cat("\n", "------------- Kallisto - DESeq2 - RNASeq Differential expression ---------", "\n")
#### read data ####
#### Read counts for each sample
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
colnames(counts) = c("GeneID", samples$samples)

## aggregate to gene level
counts$GeneID = substr(counts$GeneID, 1, 10) #get first 10 letters, (remove -RA,-RB etc of transcripts)
counts = counts %>% group_by(GeneID) %>% summarise_all(sum)

### count total reads counted
counts = counts %>% column_to_rownames('GeneID')
#### get count statistics for each sample and plot ####
count_stats = apply(counts, 2, sum) %>% enframe(name="Sample", value="total_counts") # total counts
ngenes = nrow(counts)
count_stats$genes_zerocounts = apply(counts, 2, function(x){sum(x==0)}) # genes with zero counts 
count_stats$genes_lessthan10counts = apply(counts, 2, function(x){sum(x<10)}) # genes with less than 10 counts
count_stats = count_stats %>% mutate("proportion_zero" = genes_zerocounts/ngenes,
                                     "proportion_low" = genes_lessthan10counts/ngenes)
count_stats %>% fwrite(., "results/quant/count_statistics.tsv",sep="\t")

print("Counting and plotting total reads per sample...")
pdf("results/quant/total_reads_counted.pdf")
ggplot(count_stats, aes(x=Sample, y=total_counts, fill=samples$treatment)) + 
  geom_bar(stat='identity') + 
  theme_light() +
  ggtitle("Total reads counted (mapped to Ag transcriptome (PEST))") +
  theme(axis.text.x = element_text(angle=45)) +
  theme(plot.title = element_text(hjust = 0.5))
null = dev.off() 
 
# round numbers to be whole
counts = counts %>% rownames_to_column('GeneID') %>% 
  mutate_if(is.numeric, round) %>% column_to_rownames('GeneID')

############ Plots with all data ########################
res = vst_pca(counts, samples, colourvar = 'strain', name="PCA")
vstcounts = res[[1]]
normcounts = res[[3]]

#calculate correlations between samples based on the count data 
correlations = cor(vstcounts)

### write out raw and normalised counts 
counts %>% rownames_to_column("GeneID") %>% fwrite(., "results/quant/rawcounts.tsv", sep="\t", row.names = FALSE)
normcounts %>% as.data.frame() %>% rownames_to_column("GeneID") %>% round_df(., 1) %>% fwrite(., "results/quant/normcounts.tsv", sep="\t", row.names = FALSE)

#### 
pdf("results/plots/heatmap_correlations.pdf")
pheatmap(correlations)
garbage = dev.off()

if(!("lab" %in% colnames(samples)))
{
  samples$lab = 'FALSE'
}

######### subset data and run DESeq for each combo we need, store in xlsx ########
results_list = list()
names_list = list()
# for each strain in the dataset, make PCA plot 
# analysis of case v control with DESEq2 and make PCAs and volcano plots.
for (sp in unique(samples$species)){
  df_samples = samples %>% filter(species == sp)
  df_samples = df_samples %>% filter(lab == 'FALSE') #remove lab samples
  #print(lab_samples)
  for (st in unique(df_samples$strain)){
    ai_list = list()
    ai_list[[st]] = df_samples[df_samples$strain == st,]
    
    #if the strain has both case and control (i.e. exposed v unexposed)
    if (length(unique(ai_list[[st]]$cohort)) > 1){
      cat(glue("\n Running PCA for all {st} samples \n"))
      #do a pca for each strain 
      subcounts = counts[colnames(counts) %in% ai_list[[st]]$samples]
      #perform PCA on the data at strain level
      vstcounts = vst_pca(subcounts, ai_list[[st]], 'treatment', "PCA_", st=st)[[1]]
    }
  }
}


for (cont in contrasts){
  control = str_split(cont, "_")[[1]][1] #get first of string, which is control 
  case = str_split(cont, "_")[[1]][2] #get case 
  controls = which(samples$treatment %in% control)
  cases = which(samples$treatment %in% case)
  
  ##subset to subcounts
  subcounts = counts[,c(controls, cases)]
  subsamples = samples[c(controls, cases),]

  #make treatment a factor with the 'susceptible' as reference
  subsamples$treatment = as.factor(as.character(subsamples$treatment))
  subsamples$treatment = relevel(subsamples$treatment, control)
  cat("\n", glue("--- Running DESeq2 differential expression analysis on {cont} ---"), "\n")
  
  dds = vst_pca(subcounts, subsamples, colourvar='treatment', "PCA_", comparison=cont)[[2]]
  #set pvalue threshold and perform  DE with each comparison
  p_threshold=0.05
  cds = nbinomWaldTest(dds)
  results = results(cds) %>% as.data.frame() 
  results = results[order(results$padj),] #order by pvalue 
  results = results %>% rownames_to_column("GeneID") %>% mutate("FC" = (2^log2FoldChange))
  
  ###absolute difference
  #### Get rowsums of counts, grouping by case/control. Then get difference of counts and join with DE results
  readdiff = data.frame(t(rowsum(t(subcounts), group = subsamples$treatment, na.rm = T))) #transpose and get rowsums for each group
  readdiff$absolute_diff = readdiff[,case] - readdiff[,control] #get difference
  readdiff = data.frame(readdiff) %>% rownames_to_column('GeneID')
  results = unique(left_join(results, readdiff[,c('GeneID','absolute_diff')]))
  
  #join DE results with normal gene names
  results_list[[cont]] = unique(left_join(results, gene_names))
  fwrite(results_list[[cont]], glue("results/genediff/{cont}.csv")) #write to csv 
  
  #store names of comparisons for xlsx report sheets
  names_list[[cont]] = cont
  
  #volcano plot for each comparison, first filter to remove very lowly expressed genes 
  a=results_list[[cont]] %>% filter(`baseMean` > 20)
  pdf(glue("results/genediff/Volcano_plot_{cont}.pdf"))
  EnhancedVolcano(results,
                  lab=results$Gene_name,
                  x='log2FoldChange',
                  y='pvalue',
                  title = cont)
  garbage = dev.off()
  cat("\n", glue("{cont} complete!"), "\n")
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

sessionInfo()
