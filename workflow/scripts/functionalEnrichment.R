#GSEA RNA-Seq Ag
BiocManager::install("clusterProfiler")
BiocManager::install("gage")
BiocManager::install("KEGGREST")

library(gage)
library(KEGGREST)
library(clusterProfiler)
library(fgsea)
library(data.table)
library(glue)
library(tidyverse)
  



    ########## GO Terms #########

go = fread("resources/reference/VectorBase-48_AaegyptiLVP_AGWG_GO.gaf", sep = "\t", skip = 1, header = FALSE) %>% 
  as_tibble()

go = go[,c(2,5)]



examplePathways

############## GSEA ################

samples = fread("config/samples.tsv") %>% as.data.frame()
comparisons = fread("resources/DE.comparison.list", header = FALSE)

de = fread("../RNA_Seq_Ag/analysis/diff/ContAbo_MalaAbo.csv")
de = de[order(log2FoldChange)]
rank = de[,c('Gene_stable_ID', 'log2FoldChange')] %>% column_to_rownames('Gene_stable_ID')
######
rownames(rank) = rownames(rank) %>% str_replace("AGAP", "AgaP_AGAP")
head(rank)

#------------------- get a.gambiae kegg pathways -------------------------
kg.aga=kegg.gsets("aga")
kg.aga.gs=kg.aga$kg.sets[kg.aga$sigmet.idx]
kg.aga

#kegg pathway enrichment
kegg_enrich <- gage(rank, gsets=kg.aga.gs)
kegg_enrich

keggdf = as.data.frame(kegg_enrich$greater)

kegg_sig<-sigGeneSet(kegg_enrich, outname="aga.gage")
fwrite(rbind(kegg_enrich$greater, kegg_enrich$less), file = "aga.kegg.tsv", sep = "\t")


### fgsea ##########

rank2 = setNames(rank$log2FoldChange, rownames(rank))

str(rank2)
str(exampleRanks)
fgseaRes <- fgsea(pathways = kg.aga.gs, 
                  stats    = rank,
                  minSize  = 15,
                  maxSize  = 50,
                  nperm = 1000)

kg.aga.gs
head(str(examplePathways))
head(str(kg.aga.gs))
head(rank2)
head(exampleRanks)

  plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam=0.5)



