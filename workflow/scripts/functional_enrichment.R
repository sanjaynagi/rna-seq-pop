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
  

############## GSEA ################

samples = fread("config/samples.tsv") %>% as.data.frame()
comparisons = fread("resources/DE.comparison.list", header = FALSE)

de = fread("results/diff/ContAbo_PiriAbo.csv")

sort(de$absolute_diff, decreasing = TRUE)

de = de[order(-absolute_diff)]
de = de %>% column_to_rownames('Gene_stable_ID')
diff = de %>% select('absolute_diff')



######
data(examplePathways)
data(exampleRanks)



fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500,
                  nperm = 10000)



plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam=0.5)





#------------------- get a.gambiae kegg pathways -------------------------
kg.aga=kegg.gsets("aga")
kg.aga.gs=kg.aga$kg.sets[kg.aga$sigmet.idx]


#kegg pathway enrichment
Ag_gage <- gage(diff, gsets=kg.aga.gs)

Ag_gage$greater

?gage
conv <- keggConv("aga", "ncbi-geneid")
head(conv)
?keggConv

for (com in comparisons){
}

write.table(aga.gage$greater, file = "aga.kegg.txt", sep = "\t")
write.table(rbind(aga.gage$greater, aga.gage$less), file = "aga.kegg.txt", sep = "\t")

aga.gage.sig<-sigGeneSet(aga.gage, outname="aga.gage")
