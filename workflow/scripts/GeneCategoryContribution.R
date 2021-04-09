#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

### Contribution of IR gene categories to overall gene expression
library(tidyverse)
library(data.table)
library(glue)

##### define functions ######
round_df = function(df, digits) {
  
  #' This function rounds all the numeric columns of a data.frame
  nums = vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] = round(df[,nums], digits = digits)
  (df)
}

samples = fread(snakemake@input[['samples']], sep="\t") %>% as.data.frame()
#samples = fread("config/samples.tsv", sep="\t") %>% as.data.frame()
gene_names = fread("resources/gene_names.tsv", sep="\t")[1:14459,] %>% distinct()

# Read in normalised counts and sum biological replicates
counts = fread("results/quant/normcounts.tsv", sep="\t") %>% column_to_rownames("GeneID")
counts = data.frame(t(rowsum(t(counts), group = samples$treatment, na.rm = T)))

# Get total counts for each treatment
totalCounts = colSums(counts)

# Get lists of P450s, COEs and GSTs based on gene description
p450genes = gene_names[str_detect(gene_names$Gene_description, "P450"),] %>% select(Gene_stable_ID, Gene_name)
p450genes = p450genes[!p450genes$Gene_name %in% c("CYP4G16", "CYP4G17")] %>% distinct() # remove these genes from P450 column as I have added to cuticular
coegenes = gene_names[str_detect(gene_names$Gene_description, "carboxylesterase"),] %>% select(Gene_stable_ID)
coegenes = bind_rows(coegenes, data.frame("Gene_stable_ID" = "AGAP006227"))
gstgenes = gene_names[str_detect(gene_names$Gene_description, c("glutathione S"))] %>% select(Gene_stable_ID, Gene_name, Gene_description)


genegrps = list()
for (file in list.files("resources/annotation/")){
  group = str_remove(file, ".txt") 

  genegrps[[group]] = fread(glue("resources/annotation/{file}"), sep="\t") %>% 
    rename("Gene_stable_ID" = "Gene stable ID") %>% 
    select(`Gene_stable_ID`) %>% distinct()

}

# Add manually searched groups
genegrps[["P450s"]] = p450genes
genegrps[["COEs"]] = coegenes
genegrps[["GSTs"]] = gstgenes


# Loop through groups and calculate % contribution to total counts
perc_contrib = list()
for (group in names(genegrps)){
  subcounts = counts %>% 
    rownames_to_column("GeneID") %>% 
    filter(GeneID %in% genegrps[[group]]$`Gene_stable_ID`) %>% 
    column_to_rownames("GeneID")

  perc_contrib[[group]] = data.frame((colSums(subcounts) / totalCounts) *100) %>% 
    rownames_to_column("sample")
  colnames(perc_contrib[[group]]) = c("sample", group)
  
}

# Join dfs
perc_contrib_df = purrr::reduce(perc_contrib, merge) %>% round_df(3)

# Loop through each category and plot % contribution to total counts
for (group in names(genegrps)){
  plt =  ggplot(perc_contrib_df, aes_string(x="sample", y=glue("{group}"), fill="sample")) + 
    geom_bar(stat="identity") + 
    ggtitle(glue("% of read counts {group}")) +
    theme_light()
  
  pdf(glue("results/plots/percentageContribution_{group}.pdf"))
  print(plt)
  dev.off()
}

# Write to file 
fwrite(perc_contrib_df, 
       "results/quant/percentageContributionGeneCategories.tsv", sep="\t")



sessionInfo()
