#### Description #### 
# Reformat NetMHCpan result and Put mutation and wild result together. 

#### Hyperparameter and library #### 
input_path <- "./output/pos/"
output_path <- "./output/pos/"
library(tidyverse)

#### Reformat netMHCpan result ##### 
mutate.input <- read.delim(file = paste0(output_path,"pos.mutate.peptide.out.1"),sep = "\n")
wild.input <- read.delim(file = paste0(output_path,"pos.wild.peptide.out.1"),sep = "\n")
fraction_reformat <- function(fraction.mutate){
  fraction.mutate.1 <- as.data.frame(fraction.mutate[!str_detect(string = fraction.mutate[,1],pattern = "-----"),])
  fraction.mutate.2 <- str_split_fixed(string = fraction.mutate.1[,1],pattern = "\\s+",n = Inf)
  fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
    select(!V1) 
  colnames(fraction.mutate.3) <- NULL
  rownames(fraction.mutate.3) <- NULL
  return(fraction.mutate.3)
}
mutate.output <- fraction_reformat(mutate.input)
wild.output <- fraction_reformat(wild.input)
write.table(mutate.output,
            file = paste0(output_path,"mutate.netMHCoutput.tsv"),
            sep = "\t",row.names = F,col.names = F,quote = F)
write.table(wild.output,
            file = paste0(output_path,"wild.netMHCoutput.tsv"),
            sep = "\t",row.names = F,col.names = F,quote = F)

#### Put wild and mutation result together #### 
reformat.fraction.mutation <- reformat.fraction.mutation %>% 
  unite(col = "BindLevel",c(V17,V18),sep = "",remove = TRUE)
colnames(reformat.fraction.mutation) <- reformat.fraction.mutation[1,]
reformat.fraction.mutation <- reformat.fraction.mutation[-1,]
reformat.fraction.mutation <- reformat.fraction.mutation %>% 
  select("MHC","Peptide","Aff(nM)","BindLevel")
colnames(reformat.fraction.mutation) <- paste0("mut_",colnames(reformat.fraction.mutation))

reformat.fraction.wild <- read.table(paste0(input_path,"wild.netMHCoutput.tsv"),
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

write.table(pos.neoantigens.filtered,file = paste0(output_path,"final.pos.neoantigens.filtered.tsv"),
            sep = "\t",row.names = F,quote = F)
