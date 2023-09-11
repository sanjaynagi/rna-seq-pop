#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Generate freebayes params ##
# 1) make bed for parallelising freebayes
library(dplyr)
library(data.table)
library(glue)


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

# read inputs
contigs = snakemake@params[['contigs']]
chunks = snakemake@params[['chunks']]
fai = fread(snakemake@input[['index']])

# select contigs we want, and start, end columns
fai = fai[fai$V1 %in% contigs, c(1,2)]

# for each chromsome
for (contig in contigs){
   #subset index to desired contig
   f = fai[fai$V1 == contig]
   #get sequence of n chunks from 0 to length of contig
   bedseq = round(seq(0, f$V2, length.out = chunks))
   
   #for each chunk
   for (i in 1:(chunks-1)){
      #write bed file, one for each chunk/interval, which will be passed as input to freebayes
      row = c(contig, bedseq[i], bedseq[i+1])
      data.frame(row) %>% t() %>% fwrite(., glue("resources/regions/genome.{contig}.region.{i}.bed"), sep="\t", col.names = FALSE)
   }
}



# 2) Make bamlist and populations.tsv file

metadata = load_metadata(snakemake@params[['metadata']])
metadata$bams = paste0("results/alignments/", metadata$sampleID,".bam")

metadata %>% 
   select(bams, strain) %>% 
   fwrite(., snakemake@output[['pops']], sep="\t", row.names = FALSE, col.names = FALSE)

metadata %>% 
   select(bams) %>% 
   fwrite(.,snakemake@output[['bamlist']], sep="\t", row.names = FALSE, col.names = FALSE)

sessionInfo()