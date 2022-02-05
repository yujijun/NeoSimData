########### Description ####
# create time: Fri Jul 23 14:56:19 2021
# Description: The script is used to generate simulated neopeptides.
Main <- function(ProteinSeq = "./input/uniprot-proteome_UP000005640.fasta",
                 pepMatrix = "./input/deredundancy.allhuman.tsv",
                 peptide_colname = "description",
                 output_name = "./output/pos/"){
  library(tidyverse)
  library(Biostrings)
  library(seqinr)
  library(data.table)
  library(foreach)
  library(parallel)
  library(doSNOW)
  # create pattern with one character as pattern
  Create.pattern.1 <- function(string){
    require(stringr)
    string.pattern.1 <- rep(string,str_length(string))
    str_sub(string.pattern.1,start = 1:str_length(string),end = 1:str_length(string)) <- "-"
    string.pattern.2 <- str_replace_all(string.pattern.1,pattern = "-",replacement = "[:alpha:]")
    names(string.pattern.2) <- rep(string,length(string.pattern.2))
    return(string.pattern.2)
  }
  # create pattern with two character as pattern
  Create.pattern.2 <- function(string){
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
  #match all string pattern list with a longer string 
  matchstrings.single <- function(patternExm, stringExm){
    require(stringr)
    require(tidyverse)
    matchResult <- stringr::str_match_all(patternExm,string = stringExm) %>% 
      unlist()
    return(matchResult)
  }
  # A function for formating of final result
  fmt_match <- function(matchrst = final_match){
    matchrst_dt <- as.data.frame(matchrst) %>% 
      rownames_to_column() %>% 
      magrittr::set_colnames(c("name","wild_peptides")) %>% 
      as.data.table()
    matchrst_dt[,repnum := stringr::str_extract(string = name,pattern = "\\d")]
    matchrst_dt[,pattern_NoNum := stringr::str_remove(string = name,pattern = "\\d")]
    matchrst_dt <- matchrst_dt %>% 
      separate(col = pattern_NoNum,
               into = c("mutate_peptides","mutate_peptides_pattern"),
               sep = "_",remove = F) %>% 
      filter(wild_peptides != mutate_peptides) %>% 
      unique(by = c("wild_peptides","mutate_peptides")) %>% 
      select(wild_peptides,mutate_peptides)
    return(matchrst_dt)
  }
  ## input
  proteinomics <- readAAStringSet(filepath = ProteinSeq)
  peptideMatrix <- fread(input = pepMatrix)
  #### 1  preprocess of protein sequence 
  omics_seq <- as.character(proteinomics,use.names =F)
  omics_name <- names(proteinomics)
  omics_seq <- as.data.frame(omics_seq) 
  omics_matrix <- omics_seq %>% 
    mutate(gene_name = omics_name,.before = 1) %>% 
    mutate(seq_length = str_length(omics_seq))
  omics_matrix <- omics_matrix %>% 
    mutate(sum_length = cumsum(seq_length))
  Allstring <- stringr::str_c(omics_matrix$omics_seq,collapse = "")
  #### 2 preprocess of peptide list 
  peptide_list <- peptideMatrix[[peptide_colname]] %>% 
    as.list()
  #### 3 match by pattern 
  All.pattern <- purrr::map(peptide_list,Create.pattern.1) %>% 
    unlist() %>% as.list()
  parameterDT <- as.data.table(unlist(All.pattern)) %>% 
    magrittr::set_colnames(c("Pattern"))
  coreN <- parallel::detectCores(logical=F)
  cl <- parallel::makeCluster(coreN, type="PSOCK")
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
  sink()
  opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
  dt_peptdesc <- foreach::foreach(i=1:nrow(parameterDT), .inorder=T, .options.snow=opts)%dopar%{
    matchstrings.single(parameterDT$"Pattern"[[i]], Allstring)
  } 
  listname <- paste0(names(All.pattern),"_",unlist(All.pattern))
  names(dt_peptdesc) <- listname
  final_match <- unlist(dt_peptdesc)
  close(pb)
  parallel::stopCluster(cl)
  gc();gc()
  #### 4 format matching result
  final_result <- fmt_match(final_match)
  final_result <- final_result %>% 
    left_join(peptideMatrix,by = c("mutate_peptides"= peptide_colname))
  #### 5 output 
  save(final_result,file = paste0(output_name))
  #return(final_result)
}

## 

Main(pepMatrix = "./input/deredundancy.allclean.tsv",
             output_name = "./output/pos/allcleanNeoStim.Rdata")

Main(pepMatrix = "./input/deredundancy.allhuman.tsv",
     output_name = "./output/pos/allhumanNeoStim.Rdata")

Main(pepMatrix = "./input/deredundancy.allmouse.tsv",
     output_name = "./output/pos/allmouseNeoStim.Rdata")


## load and format for Pos_netMHcpanRunning.sh running ####

load("./output/pos/allcleanNeoStim.Rdata")
wild_seq <- final_result$wild_peptides
wild_name <- paste0("Seq_",1:length(wild_seq))
write.fasta(sequences = as.list(wild_seq),
            names = wild_name,
            file.out =  paste0(output_path,"pos.fraction.wild.fsa"))

mutate_seq <- final_result$mutate_peptides
mutate_name <- paste0("Seq_",1:length(mutate_seq))
write.fasta(sequences = as.list(mutate_seq),
            names = mutate_name,
            file.out =  paste0(output_path,"pos.fraction.mutate.fsa"))
