#### Description #### 
# Put mutation and wild peptides together 

#### hyperparameter and library #### 
input_path <- ""
output_path <- ""
library(tidyverse)

#### Main #### 
reformat.fraction.mutation <- read.table("./output_20210824/mutate.netMHCoutput.tsv",
                                         sep = "\t")
reformat.fraction.mutation <- reformat.fraction.mutation %>% 
  unite(col = "BindLevel",c(V17,V18),sep = "",remove = TRUE)
colnames(reformat.fraction.mutation) <- reformat.fraction.mutation[1,]
reformat.fraction.mutation <- reformat.fraction.mutation[-1,]
reformat.fraction.mutation <- reformat.fraction.mutation %>% 
  select("MHC","Peptide","Aff(nM)","BindLevel")
colnames(reformat.fraction.mutation) <- paste0("mut_",colnames(reformat.fraction.mutation))

reformat.fraction.wild <- read.table("./output_20210824/wild.netMHCoutput.tsv",
                                     sep = "\t")
reformat.fraction.wild <- reformat.fraction.wild %>% 
  unite(col = "BindLevel",c(V17,V18),sep = "",remove = TRUE)
colnames(reformat.fraction.wild) <- reformat.fraction.wild[1,]
reformat.fraction.wild <- reformat.fraction.wild[-1,]
reformat.fraction.wild <- reformat.fraction.wild %>% 
  select("MHC","Peptide","Aff(nM)","BindLevel")
colnames(reformat.fraction.wild) <- paste0("wild_",colnames(reformat.fraction.wild))

pos.neoantigens <- cbind(reformat.fraction.wild,
                         reformat.fraction.mutation)

pos.neoantigens.relocate <- pos.neoantigens %>% 
  mutate(mut_MHC = NULL) %>% 
  relocate(mut_Peptide,.after = wild_Peptide) %>% 
  relocate(mut_BindLevel,.after = wild_BindLevel) %>% 
  relocate(wild_MHC,.after = mut_Peptide) %>% 
  relocate(`mut_Aff(nM)`,.before = mut_BindLevel)

pos.neoantigens.filtered <- pos.neoantigens.relocate %>% 
  filter(wild_BindLevel == "" & mut_BindLevel != "")

write.table(pos.neoantigens.filtered,file = "./output_20210824/final.pos.neoantigens.filtered.tsv",
            sep = "\t",row.names = F,quote = F)
