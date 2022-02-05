## ---------------------------
##
## Script name:  This is script for merging of neg and pos 
##
## Purpose of script:
##
## Author: Dr. JijunYu
##
## Date Created: 2022-01-14
##
## Copyright (c) Timothy Farewell, 2022
## Email: jijunyuedu@outlook.com  

######  library and input #### 
library(tidyverse)
neg <- readr::read_delim("./output/neg/neg.neoantigen.tsv")
pos <- readr::read_delim("./output/pos/final.pos.neoantigens.filtered.tsv")

neg_filter <- neg %>% 
  select(wild_peptide,mutate_peptide,Allele,wild_ic50,wild_rank,mutate_ic50,mutate_rank) %>% 
  rename_with(tolower) %>% 
  filter(mutate_ic50 <= 500) %>% 
  mutate(type = "neg")
pos_filter <- pos %>% 
  filter(mutate_ic50<= 500) %>% 
  mutate(type = "pos")

Neoantigen.all <- rbind(pos_filter,neg_filter)
write.table(Neoantigen.all,file = "output/StimNeoantigen/neoantigen_all_0115.tsv",
            quote = F,col.names = T,row.names = F)
