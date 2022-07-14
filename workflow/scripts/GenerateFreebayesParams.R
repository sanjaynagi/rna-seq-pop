#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Generate freebayes params ##
# 1) make bed for parallelising freebayes
library(dplyr)
library(data.table)
library(glue)

# read inputs
chroms = snakemake@params[['chroms']]
chunks = snakemake@params[['chunks']]
fai = fread(snakemake@input[['index']])

# select chroms we want, and start, end columns
fai = fai[fai$V1 %in% chroms, c(1,2)]

# for each chromsome
for (chrom in chroms){
   #subset index to desired chrom
   f = fai[fai$V1 == chrom]
   #get sequence of n chunks from 0 to length of chrom
   bedseq = round(seq(0, f$V2, length.out = chunks))
   
   #for each chunk
   for (i in 1:(chunks-1)){
      #write bed file, one for each chunk/interval, which will be passed as input to freebayes
      row = c(chrom, bedseq[i], bedseq[i+1])
      data.frame(row) %>% t() %>% fwrite(., glue("resources/regions/genome.{chrom}.region.{i}.bed"), sep="\t", col.names = FALSE)
   }
}



# 2) Make bamlist and populations.tsv file

metadata = fread(snakemake@params[['metadata']], sep="\t")

metadata$bams = paste0("results/alignments/", metadata$sampleID,".bam")

metadata %>% 
   select(bams, strain) %>% 
   fwrite(., snakemake@output[['pops']], sep="\t", row.names = FALSE, col.names = FALSE)

metadata %>% 
   select(bams) %>% 
   fwrite(.,snakemake@output[['bamlist']], sep="\t", row.names = FALSE, col.names = FALSE)

sessionInfo()