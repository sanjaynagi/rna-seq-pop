#! /usr/bin/env RScript
#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

library(data.table)
library(kissDE)
library(glue)
library(tidyverse)
library(rtracklayer)

######## parse files #############
chroms = snakemake@params[[1]]
metadata = fread(snakemake@input[[1]])
contrasts = fread(snakemake@input[[3]], header = TRUE) 
comparisons = contrasts %>% separate(contrast, into = c("control", "case"), sep = "_")
mincounts = snakemake@params[[3]]
pval = snakemake@params[[4]]
gffchromprefix = snakemake@params[[5]]

gffpath = snakemake@input[[2]]
#load gff with rtracklayer and filter to genes only 
gff = rtracklayer::import(gffpath) %>% 
  as.data.frame() %>% 
  filter(type == 'gene') %>% 
  select(seqnames, start, end, ID) %>% 
  dplyr::rename("chrom" = 'seqnames') %>% 
  arrange(chrom, start) %>% 
  as.data.table()
#remove prefix for chrom column, so it matches the bed file 
gff$chrom = str_remove(gff$chrom, gffchromprefix) #could be snakemake param if needed

#for each contrast/comparison, perform differential SNP analysis 
for (i in 1:nrow(comparisons)){
  name = contrasts[i,]
  control = comparisons$control[i]
  case = comparisons$case[i]
  
  samples = metadata[metadata$treatment %in% c(case, control)]$samples
  print(glue("Extracting allele tables for {case}, {control}"))
  
  sample_list = list()
  for (sample in samples){

    chrom_list = list()
    for (chrom in chroms){
      #read in data 
      alleles  = fread(glue("results/variants/alleleTables/{sample}.chr{chrom}.allele.table"))
      #sum lowercase and uppercase alleles, make new column (refcount)
      chrom_list[[chrom]] = alleles %>% 
        mutate("A" = A+a,"T" = T+t,"G" = G+g,"C" = C+c) %>% 
        select(-c(a,t,c,g,V15)) %>% 
        mutate(refcount = case_when(ref == 'T' ~ T,
                                    ref == 'G' ~ G, 
                                    ref == 'A' ~ A,
                                    ref == 'C' ~ C))
    }
    #bind chroms together
    alleles = rbindlist(chrom_list)
    #this step find the most abundant ALT base, whether A,T,C or G, assigns its count as 'altcount'
    list_ = list()
    for (base in c("A", "C", "G", "T")){
      a = alleles %>% filter(ref == base) %>% select(A,T,C,G) %>% select(-base)
      b = alleles %>% filter(ref == base)
      col = colnames(a)[apply(a, 1 , which.max)]
      b$alt = col
      b$altcount = apply(a, 1, max)
      list_[[base]] = b
    }
    
    alleles = rbindlist(list_) %>% arrange(chr, loc)
    #pivot to get separate rows for refcount and altcount as preferred by kissDE package
    sample_list[[sample]] = alleles %>% select(chr, loc, ref,alt, refcount, altcount) %>% 
      pivot_longer(c(refcount, altcount),names_to="type", values_to=sample) %>% 
      mutate(eventsName = paste0("chr",chr,"_",loc,"_",ref,">",alt)) %>% 
      select(eventsName,type, sample)
    
  }
  #merge all counts across samples 
  counts = sample_list %>% 
    purrr::reduce(full_join, by=c("eventsName", "type")) %>% 
    mutate("eventsLength" = 1) %>% 
    replace(is.na(.), 0) %>% 
    select(-type) %>% 
    relocate(eventsLength, .after=eventsName) %>% 
    as.data.frame()
  
  #sum refs and alts together for each snp to use as filter
  print(glue("Filtering SNP count data at a minimum of {mincounts} summed across samples"))
  summed = counts %>% 
    group_by(eventsName, eventsLength) %>% 
    summarise_all(sum)
  #get rowsums of counts, to use as filter
  rowcounts = summed %>% ungroup %>% select(-c(eventsName, eventsLength)) %>% rowSums()
  summed = summed[rowcounts >= mincounts,]
  print(glue("Retaining {nrow(summed)} SNPs, out of {length(rowcounts)}"))
  #remove snps that have a count across samples less than mincounts(100)
  counts = counts[counts$eventsName %in% summed$eventsName,]
  
  conditionsde = c(rep(control, 3), rep(case,3))
  #rename counts columns for kissde
  colnames(counts) = colnames(counts) %>% 
    str_replace("1", "rep1") %>% 
    str_replace("2", "rep2") %>% 
    str_replace("3", "rep3")
  #run kissde quality control
  print(glue("Running kissDE QC for {name}"))
  qualityControl(counts, conditionsde, storeFigs = glue("results/variants/diffsnps/kissDEfigs_{name}"))
  
  #Run kissde algorithm to find differentially expressed variants
  de_Vars = diffExpressedVariants(counts, conditions = conditionsde)                    

  #parse results to more readable, filterable form 
  results = de_Vars[[1]] %>% separate(ID, into = c("chrom", "pos", "REF>ALT"), sep="_") %>% 
    mutate("pos" = as.numeric(pos), "chrom" = str_remove(chrom, "chr")) %>% 
    arrange(chrom, pos)
  
  ##### intersect (data.table::foverlap) gff and results to get gene names and descriptions of snps #######
  bed = results %>% 
    select(chrom, pos, Adjusted_pvalue, `Deltaf/DeltaPSI`, `REF>ALT`) %>% 
    mutate("chrom" = as.character(chrom), "start" = as.numeric(pos) -1, "end" = as.numeric(pos)) %>% 
    select(chrom, start, end, Adjusted_pvalue, `Deltaf/DeltaPSI`,`REF>ALT`, -pos) %>% 
    as.data.table()
  
  #data table fast overlaps, useful to find intersections and join  
  setkey(gff, chrom, start, end)
  de_variants = foverlaps(bed, gff, 
                          by.x=c("chrom", "start", "end"), 
                          by.y=c("chrom", "start", "end"),
                          type = "within",
                          nomatch = 0L)

  #write to file
  print(glue("Writing {name} results to results/variants/diffsnps/"))
  results %>% fwrite(., glue("results/variants/diffsnps/{name}.normcounts.tsv"), sep="\t", row.names=FALSE)
  de_variants %>% dplyr::rename("GeneID" = "ID") %>% fwrite(., glue("results/variants/diffsnps/{name}.kissDE.tsv"), sep="\t", row.names=FALSE)
      
  de_variants %>% 
    filter(Adjusted_pvalue <= pval) %>% 
    fwrite(., glue("results/variants/diffsnps/{name}.sig.kissDE.tsv"), sep="\t", row.names=FALSE)
}
