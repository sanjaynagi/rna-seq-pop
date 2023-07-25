#! /usr/bin/env RScript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


#' This script reads in tables of allele counts, produced with samtools mpileup.
#' Using the R package kissDE, it performs a differential SNP testing 
#' analysis, testing if certain alleles are biased to be found in one 
#' treatment group 

library(tidyverse)
library(data.table)
library(kissDE)
library(glue)
library(rtracklayer)


######## parse inputs #############
contigs = snakemake@params[['contigs']]
metadata = fread(snakemake@input[['metadata']])
contrasts = data.frame("contrast" = snakemake@params[['DEcontrasts']])
comparisons = contrasts %>% separate(contrast, into = c("control", "case"), sep = "_")
mincounts = snakemake@params[['mincounts']]
pval = snakemake@params[['pval_flt']]
gene_names = fread(snakemake@input[['geneNames']], sep="\t") %>% select(-TranscriptID) %>% distinct()

gffpath = snakemake@input[['gff']]
#load gff with rtracklayer and filter to genes only 
gff = rtracklayer::import(gffpath) %>% 
  as.data.frame() %>% 
  filter(type == 'gene') %>% 
  select(seqnames, start, end, ID) %>% 
  dplyr::rename("contig" = 'seqnames') %>% 
  arrange(contig, start) %>% 
  as.data.table()
#remove prefix for contig column, so it matches the bed file 

#for each contrast/comparison, perform differential SNP analysis 
for (i in 1:nrow(comparisons)){
  name = contrasts[i,]
  control = comparisons$control[i]
  case = comparisons$case[i]
  
  nrepscontrol =  sum(metadata$treatment %in% control)
  nrepscase = sum(metadata$treatment %in% case)
  
  samples = metadata[metadata$treatment %in% c(case, control),]$sampleID
  print(glue("Extracting allele tables for {case}, {control}"))
  
  sample_list = list()
  for (sample in samples){
    
    chrom_list = list()
    for (contig in contigs){
      #read in data 
      alleles  = fread(glue("results/variantAnalysis/alleleTables/{sample}.chr{contig}.allele.table"))
      #sum lowercase and uppercase alleles, make new column (refcount)
      chrom_list[[contig]] = alleles %>% 
        mutate("A" = A+a,"T" = T+t,"G" = G+g,"C" = C+c) %>% 
        select(-c(a,t,c,g,V15)) %>% 
        mutate(refcount = case_when(ref == 'T' ~ T,
                                    ref == 'G' ~ G, 
                                    ref == 'A' ~ A,
                                    ref == 'C' ~ C))
    }
    #bind contigs together
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
  
  conditionsde = c(rep(control, nrepscontrol), rep(case, nrepscase))
  #rename counts columns for kissde
  #colnames(counts) = colnames(counts) %>% 
  #  str_replace("1", "rep1") %>% 
  #  str_replace("2", "rep2") %>% 
  #  str_replace("3", "rep3")
  #run kissde quality control
  print(glue("Running kissDE QC for {name}"))
  qualityControl(counts, conditionsde, storeFigs = glue("results/variantAnalysis/diffsnps/kissDEfigs_{name}"))
  
  #Run kissde algorithm to find differentially expressed variants
  de_Vars = diffExpressedVariants(counts, conditions = conditionsde)                    
  
  #parse results to more readable, filterable form 
  results = de_Vars[[1]] %>% separate(ID, into = c("contig", "pos", "REF>ALT"), sep="_") %>% 
    mutate("pos" = as.numeric(pos), "contig" = str_remove(contig, "chr")) %>% 
    arrange(contig, pos)
  
  ##### intersect (data.table::foverlap) gff and results to get gene names and descriptions of snps #######
  bed = results %>% 
    select(contig, pos, Adjusted_pvalue, `Deltaf/DeltaPSI`, `REF>ALT`) %>% 
    mutate("contig" = as.character(contig), "start" = as.numeric(pos) -1, "end" = as.numeric(pos)) %>% 
    select(contig, start, end, Adjusted_pvalue, `Deltaf/DeltaPSI`,`REF>ALT`, -pos) %>% 
    as.data.table()
  
  # Data table fast overlaps, useful to find intersections and join  
  setkey(gff, contig, start, end)
  de_variants = foverlaps(bed, gff, 
                          by.x=c("contig", "start", "end"), 
                          by.y=c("contig", "start", "end"),
                          type = "within",
                          nomatch = 0L) %>% 
    dplyr::rename("GeneID" = "ID")

    # Write to file
  print(glue("Writing {name} results to results/variantAnalysis/diffsnps/"))
  results %>% fwrite(., glue("results/variantAnalysis/diffsnps/{name}.normcounts.tsv"), sep="\t", row.names=FALSE)

  de_variants = left_join(de_variants, gene_names) %>% distinct()
  
  fwrite(de_variants, glue("results/variantAnalysis/diffsnps/{name}.kissDE.tsv"), sep="\t", row.names=FALSE)
  
  de_variants %>% 
    filter(Adjusted_pvalue <= pval) %>% 
    fwrite(., glue("results/variantAnalysis/diffsnps/{name}.sig.kissDE.tsv"), sep="\t", row.names=FALSE)
}

sessionInfo()