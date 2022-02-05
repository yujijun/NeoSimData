#### Description #### 
# Put mutation and wild peptides together 

#### Hyperparameter and library #### 
input_path <- "./output/neg/"
output_path <- "./output/neg/"
library(tidyverse)

#### input #####
mutate.input <- read.delim(file = paste0(output_path,"fraction.mutate.fsa.out"),sep = "\t")
wild.input <- read.delim(file = paste0(output_path,"fraction.wild.fsa.out"),sep = "\t")
colnames(mutate.input) <- paste0("mutate_",colnames(mutate.input))
colnames(wild.input) <- paste0("wild_",colnames(wild.input))
mutate.input <- mutate.input %>% 
  distinct() %>% 
  mutate(id = paste0(mutate.input$mutate_seq_num,"_",mutate.input$mutate_start,"_",mutate.input$mutate_end))
wild.input <- wild.input %>% 
  distinct() %>% 
  mutate(id = paste0(wild.input$wild_seq_num,"_",wild.input$wild_start,"_",wild.input$wild_end))
all <- mutate.input %>% 
  left_join(wild.input,by = c("id" = "id")) %>% 
  select(mutate_allele,id,mutate_seq_num,mutate_start,mutate_end,mutate_length,
         mutate_peptide,wild_peptide,mutate_ic50,mutate_rank,wild_peptide,wild_ic50,wild_rank)
colnames(all) <- c("Allele","id","seq_num","start","end","length","mutate_peptide","wild_peptide","mutate_ic50",
                   "mutate_rank","wild_ic50","wild_rank")
write.table(all,file = "./output/neg/neg.neoantigen.tsv",sep = "\t",quote = T,col.names = T,row.names = T)
################ reformat netMHCpan result##### 
## reformat function 
# fraction_reformat <- function(fraction.mutate){
#   fraction.mutate.1 <- as.data.frame(fraction.mutate[!str_detect(string = fraction.mutate[,1],pattern = "-----"),])
#   fraction.mutate.2 <- str_split_fixed(string = fraction.mutate.1[,1],pattern = "\\s+",n = Inf)
#   # fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
#   #   mutate(V1 = rep(protein.mutate %>%  pull(V1),times = protein.times))
#   fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
#     select(!V1)
#   colnames(fraction.mutate.3)[1:14] <- c("Pos","MHC","Peptide",
#                                          "Core","Of","Gp","G1",
#                                          "Ip","I1","Icore","Identity",
#                                          "Score_EL","%Rank_EL","BindLevel")  
#   return(fraction.mutate.3)
# }
# ##Running 
# mutate.result <- fraction_reformat(mutate.input)
# 
# wild.result <- fraction_reformat(wild.input)
# write.table(wild.result,file = paste0(output_path,"fraction.wild.fsa.out.1.Rformat.tsv"),
#             sep = "\t",row.names = F,col.names = T,quote = F)
# write.table(mutate.result,file = paste0(output_path,"fraction.mutate.fsa.out.1.Rformat.tsv"),
#             sep = "\t",row.names = F,col.names = T,quote = F)
# 
# #### merge mutation and wild result ####
# ## mutation peptides
# colnames(reformat.fraction.mutation) <- reformat.fraction.mutation[1,]
# reformat.fraction.mutation <- reformat.fraction.mutation[-1,]
# reformat.fraction.mutation <- reformat.fraction.mutation %>% 
#   select("MHC","Peptide","Identity","Score_EL","%Rank_EL","BindLevel","V16") %>% 
#   unite(col = "BindLevel", c(BindLevel,V16),sep = " ",remove = TRUE)
# colnames(reformat.fraction.mutation) <- paste0("mut_",colnames(reformat.fraction.mutation))
# ## wild peptide 
# colnames(reformat.fraction.wild) <- reformat.fraction.wild[1,]
# reformat.fraction.wild <- reformat.fraction.wild[-1,]
# reformat.fraction.wild <- reformat.fraction.wild %>% 
#   select("MHC","Peptide","Identity","Score_EL","%Rank_EL","BindLevel","V16") %>% 
#   unite(col = "BindLevel", c(BindLevel,V16),sep = " ",remove = TRUE)
# colnames(reformat.fraction.wild) <- paste0("wild_",colnames(reformat.fraction.wild))
# 
# neg.neoantigens <- cbind(reformat.fraction.wild,
#                          reformat.fraction.mutation)
# 
# neg.neoantigens.relocate <- neg.neoantigens %>% 
#   mutate(mut_MHC = NULL) %>% 
#   relocate(mut_Peptide,.after = wild_Peptide) %>% 
#   relocate(mut_Identity,.after = wild_Identity) %>% 
#   relocate(mut_Score_EL,.after = wild_Score_EL) %>% 
#   relocate(`mut_%Rank_EL`,.after = `wild_%Rank_EL`) %>% 
#   relocate(mut_BindLevel,.after = wild_BindLevel) %>% 
#   relocate(wild_MHC,.after = mut_Peptide)
# 
# neg.neoantigens.filtered <- neg.neoantigens.relocate %>% 
#   group_by(wild_Identity) %>% 
#   arrange(desc(`mut_%Rank_EL`),desc(`wild_%Rank_EL`)) %>% 
#   slice(1) %>% 
#   ungroup()

# write.table(neg.neoantigens.filtered,file = "./output/final.neg.neoantigens.filtered.tsv",
#             sep = "\t",row.names = F,quote = F)
  