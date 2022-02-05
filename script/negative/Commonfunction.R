#### 0 code convert ####
codon <- c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA",
           "TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
           "CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT",
           "ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC",
           "GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG",
           "UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG","UAU","UAC","UAA",
           "UAG","UGU","UGC","UGA","UGG","CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG",
           "CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG","AUU","AUC","AUA","AUG","ACU",
           "ACC","ACA","ACG","AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG","GUU","GUC",
           "GUA","GUG","GCU","GCC","GCA","GCG","GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG")

amino <- c("F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
           "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
           "T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
           "D","D","E","E","G","G","G","G",
           "F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
           "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
           "T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
           "D","D","E","E","G","G","G","G")

trans_from <- c("A", "T", "G", "C")
trans_to   <- c("T", "A", "C", "G")
names(trans_from) <- trans_to



#### 1. 从GFF文件中匹配出gene_name所对应转录组信息 ####
gene_match_fun <- function(gene_name,ncbi_GFF_clean){
  # info 从GFF文件中匹配出gene_name所对应转录组信息
  # input string gene_name 
  # output dataframe matched transcript_id of gene name 
  match_id <- ncbi_GFF_clean %>% 
    filter(gene %in% gene_name) %>% 
    select(gene,transcript_id,type,start,end,strand) %>% 
    na.omit() %>% 
    distinct() 
  return(match_id)
}

#### 2.1 从ncbi_mrna_name_matrix 中为df中的transcript_id 匹配相应seq 和 seqlength信息 ####
# 并将匹配结果添加会df
match_seq_fun <- function(df,ncbi_mrna_name_matrix){
  # info 用transcript_id 从ncbi_mrna_name_matrix 匹配对应的序列信息
  # input  dataframe include transcript_id from output of gene_match_fun function 
  # output dataframe of match_id_merge added seq info
  match_matrix <- ncbi_mrna_name_matrix[df$transcript_id,]
  #match_matrix <- na.omit(match_matrix)
  if(nrow(na.omit(match_matrix)) == nrow(df)){
    match_seq <- match_matrix %>% pull(seq)
    match_seqLength <- match_matrix %>% pull(seqLength)
    df$seq <- match_seq
    df$mRNALength <- match_seqLength
    return(df)
  }else if(nrow(na.omit(match_matrix)) == 0){
    print("None of the sequence was matched !!!!")
    df$seq <- "None of the sequence was match"
    df$mRNALength <- 0
    return(df)
  }else{
    print("Part of the transcript_id was matched")
    match_seq <- match_matrix %>% pull(seq)
    match_seqLength <- match_matrix %>% pull(seqLength)
    df$seq <- match_seq
    df$mRNALength <- match_seqLength
    return(df)
  }
}


#### 3.1 判断外显子组的和是否等于对应转录本的长度 ####
Length_match_fun <- function(df) {
  # info 查看一个transctript_id 对应的所有的exon的区间和 是否 等于mRNA的string长度
  # input a dataframe from match_id_merge_mrna_new list 
  # output a same dataframe but add one more column T or F
  df <- df %>% 
    mutate(fragmentLength = end - start + 1,.after = end)
  df_update <- df %>%   
    filter(type == "exon") %>% 
    group_by(ID,transcript_id) %>%
    mutate(exonFragmentSum = sum(fragmentLength)) %>% 
    mutate(judge1 = dplyr::if_else(mRNALength == exonFragmentSum,
                                   true = T,
                                   false = F))
  return(df_update)
}

#### 3.2 判断突变位置是否在转录本区间内 ####
pos_judge_fun <- function(df){
  #info start <pos < end  in df 
  #input a vector of pos, a dataframe from match_id_merge_mrna_new list 
  #output df add one more column with judge of pos
  df_new <- df %>% 
    mutate(judge2 = dplyr::if_else(start <= POS & POS <= end,T,F))
  return(df_new)
}
#### 3.3 判断vcf突变位置的参考氨基酸 和 转录本对应位置的氨基酸是否相同 ####
stringsame_judge_fun <- function(df){
  ## info 
  # 1. calculate position of vcf pos in exon 
  # 2. string_find <- find the string str_sub(match_id_merge_mrna_new[[1]]$seq[1],start = 2312,end = 2312)
  # 3. string_find  == ref or not
  # input: ref a string; a dataframe from match_id_merge_mrna_new list 
  # output: df add one more column with judge of ref string
  df_new <- df %>% 
    relocate(strand,.after = ALT) %>% 
    group_by(ID,transcript_id) %>% 
    arrange(ID,transcript_id,start) %>% 
    mutate(startStr = cumsum(c(1,as.vector(na.omit(lag(fragmentLength))))),.after = end) %>% 
    mutate(endStr = cumsum(fragmentLength),.after = startStr) 
  
  df_new2 <- df_new %>% 
    mutate(vcfRefStringPos = if_else(judge2 == TRUE,POS - start + startStr,NULL),.after = endStr) %>% 
    mutate(vcfReffind = str_sub(string = seq,
                                start = vcfRefStringPos,
                                end = vcfRefStringPos),.after = REF) %>% 
    mutate(judge3 = if_else((REF == vcfReffind & strand == "+"),T,F),.after = judge2)
    #mutate(judge3 = if_else(judge3 == F & REF == names(grep(vcfReffind,x = trans_from,value =T)) & strand == "-"),
    #       T,F)
    vcfReffind_transfer <- function(vcfReffind){
      if(!is.na(vcfReffind)){
        vcfReffind_transfer_tmp <- names(grep(vcfReffind,x = trans_from,value =T))
      }else{
        vcfReffind_transfer_tmp <- NA
      }
      return(vcfReffind_transfer_tmp)
    }

  vcfReffind_transfer_result <- purrr::map(as.list(df_new2$vcfReffind),vcfReffind_transfer)
  df_new2$vcfReffind_transfer = vcfReffind_transfer_result
  df_new2 <- df_new2 %>%
    relocate(vcfReffind_transfer,.after = vcfReffind)
  df_new3 <- df_new2 %>%
    mutate(judge3 = if_else(judge3 == T | (judge3 == F & strand == "-" & REF == vcfReffind_transfer),T,F),.after = judge2)
  return(df_new3)
}
#(REF == names(grep(vcfReffind,x = trans_from,value =T)) & strand == "-")
#### 3.4 strand "-" 处理
#stringsame_judge_fun_strand <- function(df){
#   df$judge3[which(df$strand == "-" & is.na(vcfReffind))]
# }
#### 3.4 添加蛋白质突变位点 #####
add_mutate_pos <- function(df){
  df_new <- df %>% 
    mutate(peptideMutatePos = ceiling(vcfRefStringPos/3),.after = vcfRefStringPos)
  return(df_new)
}
#### 4.1 根据突变位点，将orginal ref 转换成 mutate ref ####
make_mutated_dna <- function(df){
  mutate_fun <- function(ALT,orginal_seq,RefPos){
    alt_split <- stringr::str_split(ALT,pattern = ",")[[1]]
    alt_length <- length(alt_split)
    alt_seq <- NULL
    for(i in 1:alt_length){
      orginal_seq_tmp <- orginal_seq
      str_sub(orginal_seq_tmp,start = RefPos,end = RefPos) <- alt_split[i]
      alt_seq <- c(alt_seq,orginal_seq_tmp)
    }
    alt_seq_connect <- stringr::str_c(alt_seq,collapse = ",")
    return(alt_seq_connect)
  }
  df_new <- df %>% 
    rowwise() %>% 
    mutate(alt_seq = mutate_fun(ALT,seq,vcfRefStringPos))
  df_new2 <- df_new %>% 
    tidyr::separate(col = alt_seq,into = c("mutateseq1","mutateseq2","mutateseq3"),sep = ",")
  return(df_new2)
}

#### 4.2 将 ref DNA sequence 转换成 proteins for wild peptide ####
convertDNAintoProtein<- function(df){
  require(tidyverse)
  # generate normal peptide 
  string_transfer <- function(string){
    if(!is.na(string)){
      protein <- NULL
      while(nchar(string) >= 3){
        protein <- c(protein, amino[match(substr(string, 1, 3), codon)])
        string <- substr(string, 4, nchar(string))
      }
      protein_c <- stringr::str_c(protein,collapse = "")
      return(protein_c)
    }else{
      protein_c <- NA
      return(protein_c)
    }
  } 
  strings_transfer <- function(strings){
    peptides <- purrr::map(as.list(strings),string_transfer)
    return(peptides)
  }

  df_new <- df %>% 
    mutate(across(.cols = matches("seq"),.fns = ~strings_transfer(.x),
                  .names = "{.col}_proteins")) #就是这个across函数写的有问题
  return(df_new)
}

#### bioparallele ####
convertDNAintoProtein.para<- function(df){
  require(tidyverse)
  require(foreach)
  require(doSNOW)
  library(parallel)
  string_transfer.s <- function(strings){
    # code 
    codon <- c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA",
               "TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
               "CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT",
               "ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC",
               "GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG",
               "UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG","UAU","UAC","UAA",
               "UAG","UGU","UGC","UGA","UGG","CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG",
               "CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG","AUU","AUC","AUA","AUG","ACU",
               "ACC","ACA","ACG","AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG","GUU","GUC",
               "GUA","GUG","GCU","GCC","GCA","GCG","GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG")
    
    amino <- c("F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
               "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
               "T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
               "D","D","E","E","G","G","G","G",
               "F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
               "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
               "T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
               "D","D","E","E","G","G","G","G")
    
    trans_from <- c("A", "T", "G", "C")
    trans_to   <- c("T", "A", "C", "G")
    names(trans_from) <- trans_to
    
    # generate normal peptide #
    string_transfer <- function(string){
      if(!is.na(string)){
        protein <- NULL
        while(nchar(string) >= 3){
          protein <- c(protein, amino[match(substr(string, 1, 3), codon)])
          string <- substr(string, 4, nchar(string))
        }
        protein_c <- stringr::str_c(protein,collapse = "")
        return(protein_c)
      }else{
        protein_c <- NA
        return(protein_c)
      }
    } 

    coreN = parallel::detectCores(logical=F)
    parameterDT <- data.table::data.table(strings) %>%
      magrittr::set_colnames(c("Strings"))
    cl <- parallel::makeCluster(coreN,type = "PSOCK")
    doSNOW::registerDoSNOW(cl)
    sink(tempfile())
    pb <- pbapply::timerProgressBar(max=nrow(parameterDT), style=1)
    sink()
    opts <- list(progress=function(n){pbapply::setTimerProgressBar(pb, n)})
    Strings_proteins <- foreach::foreach(i=1:nrow(parameterDT), .inorder=T, .options.snow=opts)%dopar%{
      string_transfer(parameterDT$"Strings"[[i]])
    } %>% unlist()  #Attention, this could not be unique
    close(pb)
    parallel::stopCluster(cl)
    gc();gc()
    return(Strings_proteins)
  }
  df_new <- df %>% 
    mutate(across(.cols = matches("seq"),.fns = ~string_transfer.s(.x),
                  .names = "{.col}_proteins")) #就是这个across函数写的有问题
  return(df_new)
} 

#### 4.3 将突变的DNA序列 转变为突变的protein 序列 for mutated peptide ####
#make_mutated_peptide <- function(df){
  # 学会将一个复杂的问题拆分成几个小函数
#   make_mutated_peptide1 <- function(dna_trans_mut, amino = amino, codon = codon){
#     peptide <- NULL
#     while(nchar(dna_trans_mut) >= 3){
#       a <- amino[match(substr(dna_trans_mut, 1, 3), codon)]
#       peptide <- c(peptide, a)
#       dna_trans_mut <- substr(dna_trans_mut, 4, nchar(dna_trans_mut))
#     }
#     peptide_single <- stringr::str_c(peptide,collapse = "")
#     return(peptide_single)
#   }
#   
#   make_mutated_peptide2 <- function(seq_test){
#     if(!is.na(seq_test)){
#       seq_split <- stringr::str_split_fixed(seq_test,pattern = ",",n = Inf)
#       seq_split_length <- length(seq_split)
#       peptides <- NULL
#       for(i in seq_split_length){
#         peptide_tmp <- make_mutated_peptide1(seq_split[i],amino = amino,codon = codon)
#         peptides <- c(peptides,peptide_tmp)
#       }
#       peptides <- stringr::str_c(peptides,collapse = ",")
#       return(peptides)
#     }else{
#       peptides = NA
#       return(peptides)
#     }
#   }
#   
#   df_new <- df %>% 
#     rowwise() %>% 
#     mutate(alt_protein = make_mutated_peptide2(alt_seq))
#   return(df_new)
#   
# }
#### 5.1 将蛋白质序列剪切成需要最大肽段的长度 ####
geneProteinFraction <- function(df,max_peptide_length){
  generate_fraction <- function(proteinSeq,MutationPos){
    if(!is.na(MutationPos)){
      proteinSeq <- stringr::str_split(proteinSeq,pattern = "",n = Inf)[[1]]
      proteinSeq_start <- MutationPos - max_peptide_length + 1 
      if(proteinSeq_start < 1) proteinSeq_start <- 1
      proteinSeq_end <- MutationPos + max_peptide_length - 1
      if(proteinSeq_end > length(proteinSeq)) proteinSeq_end <- length(proteinSeq)
      peptide_choose <- proteinSeq[proteinSeq_start:proteinSeq_end]
      peptide_choose_c <- stringr::str_c(peptide_choose,collapse = "")
      return(peptide_choose_c)
    }else{
      peptide_choose_c <- NA
      return(peptide_choose_c)
    }
  }
  generate_fractions <- function(proteinSeqs,MutationPos){
    fractionProteinSeqs <- purrr::map2(proteinSeqs,as.list(MutationPos),generate_fraction)
    #fractionProteinSeqs <- unlist(fractionProteinSeqs)
    return(fractionProteinSeqs)
  }
  # df$mutateseq1_proteins_fraction <- generate_fractions(proteinSeqs = df$mutateseq1_proteins,
  #                                                          MutationPos = df$peptideMutatePos)
  # df$mutateseq2_proteins_fraction <- generate_fractions(proteinSeqs = df$mutateseq2_proteins,
  #                                                       MutationPos = df$peptideMutatePos)
  # df$mutateseq3_proteins_fraction <- generate_fractions(proteinSeqs = df$mutateseq3_proteins,
  #                                                       MutationPos = df$peptideMutatePos)
  # df$seq_proteins_fraction <- generate_fractions(proteinSeqs = df$seq_proteins,
  #                                                       MutationPos = df$peptideMutatePos)
  # 
  parameters <- df$peptideMutatePos
  df_new <- df %>%
    mutate(across(.cols = matches("proteins"),
                  .fns = ~generate_fractions(.x,MutationPos = parameters),
                  .names = "{.col}_fraction"))
  return(df_new)
}

#### 5.2 将蛋白质切片切成k-mer ####
generateKmer <- function(df,max_peptide_length){
  peptide_length = max_peptide_length
  k_mer_fun <- function(peptide_fraction,max_peptide_length = peptide_length){
    if(!is.na(unlist(peptide_fraction))){
      # peptide_fraction which is the result generate from generate_fraction function
      #peptide_fraction_split <- stringr::str_split(peptide_fraction,pattern = "",n = Inf)[[1]]
      start_para <- 1:(nchar(peptide_fraction) - max_peptide_length + 1)
      end_para <- max_peptide_length : nchar(peptide_fraction)
      peptide_all <- stringr::str_sub(peptide_fraction,start = start_para,end = end_para)
      #peptide_all <- unlist(peptide_all)
      return(peptide_all)
    }else{
      peptide_all <- NA
      return(peptide_all)
    }
  }
  
  k_mer_funs <- function(peptide_fractions){
    peptides <- purrr::map(peptide_fractions,k_mer_fun)
    #is.list(peptides)
    return(peptides)
  }
  df_new <- df %>% 
    mutate(across(.cols = matches("fraction"),
                  .fns = k_mer_funs,
                  .names = "{.col}_allpeptides"))
  return(df_new)
}

#### 6 Time record function #### 
timeRecord <- function(start_time){
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}