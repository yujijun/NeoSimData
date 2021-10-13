#### Description #### 
# Put mutation and wild peptides together 

#### Hyperparameter and library #### 
input_path <- "./output/neg/"
output_path <- "./output/neg/"
library(tidyverse)

#### Main ####
reformat.fraction.mutation <- read.table("./output/ncbi/fraction.mutate.fsa.out.1.Rformat.tsv",
                                sep = "\t")
colnames(reformat.fraction.mutation) <- reformat.fraction.mutation[1,]
reformat.fraction.mutation <- reformat.fraction.mutation[-1,]
reformat.fraction.mutation <- reformat.fraction.mutation %>% 
  select("MHC","Peptide","Identity","Score_EL","%Rank_EL","BindLevel","V16") %>% 
  unite(col = "BindLevel", c(BindLevel,V16),sep = " ",remove = TRUE)
colnames(reformat.fraction.mutation) <- paste0("mut_",colnames(reformat.fraction.mutation))

reformat.fraction.wild <- read.table("./output/ncbi/fraction.wild.fsa.out.1.Rformat.tsv",
                                sep = "\t")
colnames(reformat.fraction.wild) <- reformat.fraction.wild[1,]
reformat.fraction.wild <- reformat.fraction.wild[-1,]
reformat.fraction.wild <- reformat.fraction.wild %>% 
  select("MHC","Peptide","Identity","Score_EL","%Rank_EL","BindLevel","V16") %>% 
  unite(col = "BindLevel", c(BindLevel,V16),sep = " ",remove = TRUE)
colnames(reformat.fraction.wild) <- paste0("wild_",colnames(reformat.fraction.wild))

neg.neoantigens <- cbind(reformat.fraction.wild,
                         reformat.fraction.mutation)

neg.neoantigens.relocate <- neg.neoantigens %>% 
  mutate(mut_MHC = NULL) %>% 
  relocate(mut_Peptide,.after = wild_Peptide) %>% 
  relocate(mut_Identity,.after = wild_Identity) %>% 
  relocate(mut_Score_EL,.after = wild_Score_EL) %>% 
  relocate(`mut_%Rank_EL`,.after = `wild_%Rank_EL`) %>% 
  relocate(mut_BindLevel,.after = wild_BindLevel) %>% 
  relocate(wild_MHC,.after = mut_Peptide)

neg.neoantigens.filtered <- neg.neoantigens.relocate %>% 
  group_by(wild_Identity) %>% 
  arrange(desc(`mut_%Rank_EL`),desc(`wild_%Rank_EL`)) %>% 
  slice(1) %>% 
  ungroup()

write.table(neg.neoantigens.filtered,file = "./output/final.neg.neoantigens.filtered.tsv",
            sep = "\t",row.names = F,quote = F)
  