########### Description ####
#create time: Fri Jul 23 14:56:19 2021
########### hyperparameter and library #####
output_path <- "./output/pos"
input_path <- "./input/"
library(tidyverse)
library(Biostrings)
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
########## 1. dataset preprocess #####
#### 1.1 generate protein omics list ####
# reformat uniprot-proteome_UP000005640.fasta file
library(Biostrings)
proteinomics <- readAAStringSet(filepath = paste0(input_path,"uniprot-proteome_UP000005640.fasta"))
omics_seq <- as.character(proteinomics,use.names =F)
omics_name <- names(proteinomics)
#omics_all <- stringr::str_c(omics_seq,collapse = "")
omics_seq <- as.data.frame(omics_seq) 
omics_matrix <- omics_seq %>% 
  mutate(gene_name = omics_name,.before = 1) %>% 
  mutate(seq_length = str_length(omics_seq)) %>% 
  mutate(ncut = ntile(row_number(gene_name),10)) %>% #将打的string 切成10等分
  group_split(ncut) 
omics_matrix <- purrr::map(omics_matrix,~add_column(.x,
                                                    sum_length = cumsum(.x$seq_length)))
omics_seq_list <- purrr::map(omics_matrix,~stringr::str_c(.x$omics_seq,collapse = ""))

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

peptide_list <- Tcell_v3_human %>% as.list()

########### 2. All Needed function ####
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

########### 3. Running ####
#### 3.1 create all pattern ####
All.pattern <- purrr::map(peptide_list,Create.pattern) %>% 
  unlist() %>% as.list() %>% addname()
#### 3.2 matching ####
print("It's ready to running")
final_result <- purrr::map(omics_seq_list,~All.match.strings(dataset = .x))
save(final_result,file = paste0(output_path,"final.result.RData"))
time.all <- timeRecord(time1)
time.list <- list(time.all = time.all)
save(time.list,file = paste0(output_path,"time.record.RData"))
