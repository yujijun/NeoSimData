#### Description ######## 
# formatting result form pos_main.R and running netMHCpan
#### library and hyperparameter ##### 
library(tidyverse)
library(Biostrings)

input_path <- "./output/pos/"
output_path <- "./output/pos/"
############### 1 input and preprocess input dataset #####
#### 1.1 load final result ####
load(paste0(input_path,"final.result.RData"))
#### 1.2 convert final result list into a dataframe ####
for(i in length(final_result):1){
  if(!is.null(final_result[[i]])){
    final_result[[i]] <- final_result[[i]] %>% 
      mutate(index = i)
  }
}
finalResult.df <- do.call(rbind.data.frame,final_result) 

#### 1.3 proteinomics dataset preprocess####
library(Biostrings)
proteinomics <- readAAStringSet(filepath = "./input/uniprot-proteome_UP000005640.fasta")
omics_seq <- as.character(proteinomics,use.names =F)
omics_name <- names(proteinomics)
omics_seq <- as.data.frame(omics_seq) 
omics_matrix <- omics_seq %>% 
  mutate(gene_name = omics_name,.before = 1) %>% 
  mutate(seq_length = str_length(omics_seq)) %>% 
  mutate(ncut = ntile(row_number(gene_name),10)) %>% 
  group_split(ncut) 
omics_matrix <- purrr::map(omics_matrix,~add_column(.x,
                                                    cumsum_length = cumsum(.x$seq_length)))
omicsMatrix.df <- do.call(rbind.data.frame,omics_matrix)
omicsMatrix.df <- omicsMatrix.df %>% 
  mutate(gene_name_v1 = stringr::str_match(gene_name,pattern = "\\|\\w*_HUMAN")) %>% 
  mutate(gene_name_v1 = stringr::str_remove(gene_name_v1,pattern = "\\|")) %>% 
  mutate(gene_name_v1 = stringr::str_remove(gene_name_v1,pattern = "_HUMAN")) %>% 
  relocate(gene_name_v1,.after = gene_name) %>% 
  mutate(startPos = c(1,na.omit(lag(cumsum_length)) + 1)) %>% 
  dplyr::rename(endPos = cumsum_length) %>% 
  select(gene_name_v1,startPos,endPos,ncut,seq_length,omics_seq,gene_name)
save(omicsMatrix.df,file = paste0(output_path,"/omics_matrix.RData"))
#### 1.4 Tcell dataset preprocess #####
Tcell_v3 <- read.delim("./input/tcell_full_v3.csv",sep = ",",header = T)
Tcell_v3_clean <- Tcell_v3 %>%  
  dplyr::select(matches("Epitope")|matches("Host")|matches("Assay")|matches("MHC"))
colnames(Tcell_v3_clean) <- Tcell_v3_clean[1,]
Tcell_v3_clean <- Tcell_v3_clean[-1,]
Tcell_v3_clean <- janitor::clean_names(Tcell_v3_clean)
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(!if_any(.cols =  c(description,organism_name,name,allele_name),.fns = `==` , ""))
Tcell_v3_clean <- Tcell_v3_clean %>% 
  mutate(pep_length = stringr::str_length(description))
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(pep_length <= 50)
### add count
Tcell_v3_clean <- Tcell_v3_clean %>% group_by(description) %>% 
  add_count()
Tcell_v3_human <-Tcell_v3_clean %>% 
  filter(str_detect(name,pattern = "Homo") & 
           stringr::str_detect(qualitative_measure,pattern = "Positive") &
           pep_length == 9) 

Tcell_v3_human_peptides <- Tcell_v3_human %>% pull(description) %>% unique() #所以这里需要去除重复
Tcell_v3_human_peptides <- Tcell_v3_human_peptides[!(str_detect(Tcell_v3_human,pattern = "\\+") |
                                     str_detect(Tcell_v3_human,pattern = "\\(") |
                                     str_detect(Tcell_v3_human,pattern = "-"))]

peptide_list <- Tcell_v3_human_peptides %>% as.list()
save(Tcell_v3_human,file = paste0(output_path,"/Tcell_v3_human.RData"))





