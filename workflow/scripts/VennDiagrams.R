
library(VennDiagram)
library(data.table)
library(tidyverse)
library(glue)


#### Read data ####

#compare sig DE genes and top 5% fst genes?

de = fread("results/diff/ContAngola_FenAngola.csv")

fst = fread("results/variants/Fst_PBS.tsv")

sigde = de %>% filter(padj < 0.05)

perc = 5
highfst = fst %>% 
  arrange(desc(ContAngola_FenAngola_Fst)) %>% 
  filter(row_number() < nrow(fst) * 0.05)


