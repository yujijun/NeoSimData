########### Description ####
#create time: Fri Jul 23 14:56:19 2021
# Description: The script is used to generate simulated neopeptides.

########### hyperparameter and library #####
input_path <- "./input/"
output_path <- "./output/pos"
library(tidyverse)
library(Biostrings)
library(seqinr)
addname <- function(li){
  # addname for each list 
  for(i in 1:length(li)){
    names(li[[i]]) <- names(li)[i]
  }
  return(li)
}
time1 <- Sys.time()
timeRecord <- function(start_time){
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}
########## 1 dataset preprocess #####
#### 1.1 generate protein omics list ####
# reformat uniprot-proteome_UP000005640.fasta file
proteinomics <- readAAStringSet(filepath = paste0(input_path,"uniprot-proteome_UP000005640.fasta"))
omics_seq <- as.character(proteinomics,use.names =F)
omics_name <- names(proteinomics)
omics_seq <- as.data.frame(omics_seq) 
omics_matrix <- omics_seq %>% 
  mutate(gene_name = omics_name,.before = 1) %>% 
  mutate(seq_length = str_length(omics_seq)) %>% 
  mutate(ncut = ntile(row_number(gene_name),10)) %>% #将打的string 切成10等分
  group_split(ncut) 
omics_matrix <- purrr::map(omics_matrix,~add_column(.x,
                                                    sum_length = cumsum(.x$seq_length)))
omics_seq_list <- purrr::map(omics_matrix,~stringr::str_c(.x$omics_seq,collapse = ""))
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
#### 1.2 IEDB varified peptides list ####
# generate verified peptide list 
Tcell_v3 <- read.delim(stringr::str_c(input_path,"tcell_full_v3.csv"),sep = ",",header = T)
Tcell_v3_clean <- Tcell_v3 %>%  
  dplyr::select(matches("Epitope")|matches("Host")|matches("Assay")|matches("MHC"))
colnames(Tcell_v3_clean) <- Tcell_v3_clean[1,]
Tcell_v3_clean <- Tcell_v3_clean[-1,]
Tcell_v3_clean <- janitor::clean_names(Tcell_v3_clean)
#Delete peptides if there is a missing in description,organism_name,name,allele_name.
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(!if_any(.cols =  c(description,organism_name,name,allele_name),.fns = `==` , ""))

### Add length of peptides
Tcell_v3_clean <- Tcell_v3_clean %>% 
  mutate(pep_length = stringr::str_length(description))

### filter bigger than 50 
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(pep_length <= 50)

### add count
Tcell_v3_clean <- Tcell_v3_clean %>% group_by(description) %>% 
  add_count()

### peptides from homo sapiens and positive qualitative measure value.
Tcell_v3_human <-Tcell_v3_clean %>% 
  filter(str_detect(name,pattern = "Homo") & 
           stringr::str_detect(qualitative_measure,pattern = "Positive") &
           pep_length == 9) %>% 
  pull(description) %>% unique() #delete unique
Tcell_v3_human <- Tcell_v3_human[!(str_detect(Tcell_v3_human,pattern = "\\+") |
                                     str_detect(Tcell_v3_human,pattern = "\\(") |
                                     str_detect(Tcell_v3_human,pattern = "-"))]
save(Tcell_v3_human,file = paste0(output_path,"/Tcell_v3_human.RData"))
peptide_list <- Tcell_v3_human %>% as.list()

########### 2 All Needed function ####
#### 2.1 Function to create a string pattern ####
Create.pattern <- function(string){
  #this is a function to create all pattern of a string
  require(stringr)
  string.pattern.1 <- rep(string,str_length(string))
  str_sub(string.pattern.1,start = 1:str_length(string),end = 1:str_length(string)) <- "-"
  string.pattern.2 <- rep(string.pattern.1,str_length(string.pattern.1))
  str_sub(string.pattern.2,start = 1:str_length(string),end = 1:str_length(string)) <- "-"
  string.pattern.2 <- unique(string.pattern.2)
  string.pattern.3 <- setdiff(string.pattern.2,string.pattern.1)
  string.pattern.3 <- str_replace_all(string.pattern.3,pattern = "-",replacement = "[:alpha:]") 
  names(string.pattern.3) <- rep(string,length(string.pattern.3))
  return(string.pattern.3)
}

#### 2.2 match all string pattern list with a longer string ####
All.match.strings <- function(dataset){
  require(stringr)
  require(dplyr)
  pattern.list <- All.pattern
  match.result <- purrr::map(pattern.list,~str_match_all(string = dataset,
                                                         .x))
  for(i in length(match.result):1){
    if (identical(unlist(match.result[[i]]),character(0))){
      match.result[[i]] <- NULL
    }else{
      match.result[[i]] <- as.data.frame(match.result[[i]][[1]]) 
      match.result[[i]] <- match.result[[i]] %>% 
        mutate(mutate_peptide = names(match.result)[i]) 
    }
  }
  #match.result <- unlist(match.result)
  locate.list <- purrr::map(pattern.list,~str_locate_all(string = dataset,
                                                         .x))
  for(i in length(locate.list):1){
    if (identical(unlist(locate.list[[i]]),integer(0))){
      locate.list[[i]] <- NULL
    }else{
      locate.list[[i]] <- as.data.frame(locate.list[[i]][[1]])
    }
  }
  
  
  if(length(match.result) >= 1){
    all.list <- list()
    for(i in 1:length(match.result)){
      all.list[[i]] <- (cbind(match.result[[i]],locate.list[[i]]))
    }
    all.result <- do.call(rbind.data.frame,all.list)
    all.result <- all.result %>% 
      distinct()
  }else{
    all.result <- NULL
  }

  

  #locate.result <- str_locate_all(string = dataset,pattern = pattern)
  # locate.result <- matrix(unlist(locate.list),ncol = 2,byrow = F)
  # if(!is_empty(match.result)){
  #   all.result <- cbind(match.result,locate.result)
  #   colnames(all.result) <- c("Wild_peptide","Start","End")
  #   all.result <- as.data.frame(all.result) %>% 
  #     rownames_to_column() %>% 
  #     mutate(Mutate_peptide = stringr::str_replace(string = rowname,
  #                                                  pattern = "\\..*",
  #                                                  replacement = ""),
  #            .before = 1) %>%
  #     select(-rowname) %>% 
  #     distinct() %>% 
  #     mutate(protein = rep(names(dataset),nrow(.)))
  #   return(all.result)
  # }
  return(all.result)
}

########### 3 Running ####
#### 3.1 create all pattern ####
All.pattern <- purrr::map(peptide_list,Create.pattern) %>% 
  unlist() %>% as.list() %>% addname()
#### 3.2 matching ####
print("It's ready to running")
final_result <- purrr::map(omics_seq_list,~All.match.strings(dataset = .x))
save(final_result,file = paste0(output_path,"final.result.RData"))

########### 4 Result proprocess ########
#### 4.1 final peptides from which proteins #### 
findGeneFun <- function(start,end,index){
  gene_name <- c()
  pos <- omicsMatrix.df[which(start >= omicsMatrix.df$startPos 
                              & end < omicsMatrix.df$endPos
                              & omicsMatrix.df$ncut %in% index),]
  if(nrow(pos) == 1){
    gene_name <- c(gene_name,pos$gene_name_v1)
  }else{
    gene_name <- c(gene_name, NULL)
  }
  return(gene_name)
}
findGeneFun.v2 <- function(starts,ends,indexs){
  starts <- as.list(starts)
  ends <- as.list(ends)
  indexs <- as.list(indexs)
  l <- list(start = starts,
            end = ends,
            index = indexs)
  gene_names_list <- purrr::pmap(l,findGeneFun)
  return(gene_names_list)
}

gene_names <- findGeneFun.v2(starts = finalResult.df$start,
                             ends = finalResult.df$end,
                             indexs = finalResult.df$index)
finalResult.df$geneName <- gene_names
#### 4.2 peptides from which HLA allele #####
Tcell_v3_human_match <- Tcell_v3_human %>% 
  select(description,allele_name) %>% 
  distinct()
finalResult.df.1 <- finalResult.df %>% 
  left_join(Tcell_v3_human_match,by = c("mutate_peptide" = "description")) %>% 
  distinct()
finalResult.df.filter <- finalResult.df.1 %>% 
  select(c(V1,mutate_peptide,allele_name))
colnames(finalResult.df.1)[1] <- "wildPeptide"
colnames(finalResult.df.1)[2] <- "mutatePeptide"
save(finalResult.df.1,file = paste0(output_path,"/final_result_all.RData"))
save(finalResult.df.filter,file = paste0(output_path,"/final_result_filter.RData"))

########### 5 filter out HLA alleles ##################
finalResult.df.2 <- finalResult.df.filter[str_detect(string = finalResult.df.filter$allele_name,pattern = "\\w\\w\\w-\\w\\*\\d\\d:\\d\\d"),]
finalResult.df.3 <- finalResult.df.2 %>% 
  distinct() %>% 
  rename(wild_peptide = V1) %>% 
  mutate(allele_name = str_match(allele_name,pattern = "\\w\\w\\w-\\w\\*\\d\\d:\\d\\d"))
# take HLA-A*02:01 as an example, users could choose their own desired HLA alleles
finalResult.df.3.0201.wild <- finalResult.df.3 %>% 
  filter(allele_name == "HLA-A*02:01") %>% 
  pull(wild_peptide)
finalResult.df.3.0201.mutate <- finalResult.df.3 %>% 
  filter(allele_name == "HLA-A*02:01") %>% 
  pull(mutate_peptide)
write.table(finalResult.df.3.0201.wild,
            file = paste0(output_path,"pos.wild.peptide"),
            sep = "\t",row.names = F,col.names = F,quote = F)
write.table(finalResult.df.3.0201.mutate,
            file = paste0(output_path,"pos.mutate.peptide"),
            sep = "\t",row.names = F,col.names = F,quote = F)

########### 6 Test Running Time #### 
time.all <- timeRecord(time1)
time.list <- list(time.all = time.all)
save(time.list,file = paste0(output_path,"time.record.RData"))


########### Next: Running netMHCpan #### 
# Running `Pos_netMHCpanRuning.sh`

