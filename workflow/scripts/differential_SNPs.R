library(data.table)
library(kissDE)
library(glue)
library(tidyverse)

?kissDE
?diffExpressedVariants

######## parse files #############
chroms = c(1,2,3) #args.chrom
metadata = fread("config/samples.tsv")
names = fread("resources/DE.comparison.list", header = FALSE) 
comparisons = names %>% separate(V1, into = c("control", "case"), sep = "_")


for (i in 1:nrow(comparisons)){
  name = names[i,]
  control = comparisons$control[i]
  case = comparisons$case[i]
  
  samples = metadata[metadata$treatment %in% c(case, control)]$samples

  for (sample in samples){
    print(sample)
    chrom_list = list()
    for (chrom in chroms){
      alleles  = fread(glue("results/variants/alleleTable/{sample}.chr{chrom}.allele.table"))
      chrom_list[[chrom]] = alleles %>% 
        mutate("A" = A+a,"T" = T+t,"G" = G+g,"C" = C+c) %>% 
        select(-c(a,t,c,g,V15)) %>% 
        mutate(refcount = case_when(ref == 'T' ~ T,
                                    ref == 'G' ~ G, 
                                    ref == 'A' ~ A,
                                    ref == 'C' ~ C))
    }
    
    alleles = rbindlist(chrom_list)
    
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
    
  }
  sample_list = list()
  sample_list[[sample]] = alleles %>% select(chr, loc, ref,alt, refcount, altcount) %>% 
      pivot_longer(c(refcount, altcount),names_to="type", values_to=sample) %>% 
      mutate(eventsName = paste0("chr",chr,"_",loc,"_",ref,">",alt)) %>% select(eventsName, sample)
    
}


str(sample_list)








chrom_list = list()

for (chrom in chroms){
  alleles  = fread(glue("results/variants/alleleTable/{sample}.chr{chrom}.allele.table"))
  chrom_list[[chrom]] = alleles %>% 
    mutate("A" = A+a,"T" = T+t,"G" = G+g,"C" = C+c) %>% 
    select(-c(a,t,c,g,V15)) %>% 
    mutate(refcount = case_when(ref == 'T' ~ T,
                                ref == 'G' ~ G, 
                                ref == 'A' ~ A,
                                ref == 'C' ~ C))
  }
  
alleles = rbindlist(chrom_list)


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
  
sample_list[[sample]] = alleles %>% select(chr, loc, ref,alt, refcount, altcount) %>% 
    pivot_longer(c(refcount, altcount),names_to="type", values_to=sample) %>% 
    mutate(eventsName = paste0("chr",chr,"_",loc,"_",ref,">",alt)) %>% select(eventsName, type, sample)
  


  

?pivot_longer
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
  




