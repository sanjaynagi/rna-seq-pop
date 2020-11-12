#make bed for parallelising freebayes

library(tidyverse)
library(data.table)
library(glue)

chroms = snakemake@params[['chroms']]
fai = fread(snakemake@input[['index']])

fai = fai[fai$V1 %in% chroms, c(1,2)]

range = snakemake@params['chunks']

for (chrom in chroms){
 f = fai[fai$V1 == chrom]
 
 for (i in seq(range)){
   
   bedseq = round(seq(0, f$V2, length.out = 6))
   
   row = c(chrom, bedseq[i], bedseq[i+1])
   
   data.frame(row) %>% t() %>% fwrite(., glue("resources/regions/genome.{chrom}.region.{i}.bed"), sep="\t", col.names = FALSE)
   }
}
  
