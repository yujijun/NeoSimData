library(tidyverse)
neoantigen_all <- readr::read_delim("./output/StimNeoantigen/neoantigen_all_0115.tsv")
neoantigen_all_neg <- neoantigen_all %>% 
  filter(type == "neg")
neoantigen_all_pos <- neoantigen_all %>% 
  filter(type == "pos")
NeoStim_list <- list()
for(i in 1:1000){
  neg_tmp <- neoantigen_all_neg %>% 
    dplyr::sample_n(900)
  pos_tmp <- neoantigen_all_pos %>% 
    dplyr::sample_n(100)
  NeoStim_tmp <- rbind(neg_tmp,pos_tmp)
  NeoStim_list[[i]] <- NeoStim_tmp
}
save(NeoStim_list,file = "./output/StimNeoantigen/NeoStim_list.RData")
