biocManager::install("kissDE")

library(data.table)
library(kissDE)
library(tidyverse)

?kissDE

?diffExpressedVariants


g280s = fread("analysis/allele_balance/csvs/G280S_allele_balance.csv")

?fisher.test()


counts = data.frame("eventsName" = c("SNP1", "SNP1", "SNP2", "SNP2"),
                    "eventsLength" = as.numeric(c(100, 100, 100,100)),
                    "cond1rep1" = c(14, 50, 22, 53), 
                    "cond1rep2" = c(10, 37, 25,47),
                    "cond1rep3" = c(9, 40, 15, 45),
                    "cond2rep1" = c(44, 11, 70, 10),
                    "cond2rep2" = c(49, 6, 80, 10),
                    "cond2rep3" = c(50, 2, 64, 7))
str(counts)

conditionsde = c(rep("cond1", 3), rep("cond2",3))

dev.off()
qualityControl(counts, conditionsde)

counts[3:8]  = counts[3:8]^2


de_Vars = diffExpressedVariants(counts, conditions = conditionsde)                    


results  = de_Vars[[1]]
  




