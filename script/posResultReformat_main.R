#### Description ######## 
# 本项目主要用来规整pos 数据集输出的结果。
# Next：
#1. pos_data_main.R 函数运行结束后
#2. input: 从`1.1 load final result `位置导入final result
#2. change output_path if needed
#3. Running all code from Start into End 
#4. get output file: omics_matrix.RData\Tcell_v3_human.RData
#\final_result_all.RData\final_result_filter.RData

####-----------Start---------------####
#### library and hyperparameter ##### 
library(tidyverse)
library(Biostrings)
library(seqinr)
input_path <- "/mnt/data/yjj/pos_data/input"
output_path <- "./output_20210824/"
############### 1 input and preprocess input dataset #####
#### 1.1 load final result ####
load("./output_20210824/final.result.RData")
#### 1.2 convert final result list into a dataframe ####
for(i in length(final_result):1){
  if(!is.null(final_result[[i]])){
    final_result[[i]] <- final_result[[i]] %>% 
      mutate(index = i)
  }
}
finalResult.df <- do.call(rbind.data.frame,final_result) 

#### 1.3 proteinomics dataset preprocess [这个不需要重复运行] ####
library(Biostrings)
proteinomics <- readAAStringSet(filepath = "./uniprot-proteome_UP000005640.fasta")
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
#### 1.4 Tcell dataset preprocess [这个不要重复运行] #####
Tcell_v3 <- read.delim(stringr::str_c(input_path,"tcell_full_v3.csv"),sep = ",",header = T)
### 列选择
Tcell_v3_clean <- Tcell_v3 %>%  
  dplyr::select(matches("Epitope")|matches("Host")|matches("Assay")|matches("MHC"))
### 更改列名
colnames(Tcell_v3_clean) <- Tcell_v3_clean[1,]
Tcell_v3_clean <- Tcell_v3_clean[-1,]
Tcell_v3_clean <- janitor::clean_names(Tcell_v3_clean)
### 行选择
#如果数据中无论是肽段、或者来源微生物的名称或者是宿主的名称，或者是allele 中存在缺失，则直接将结果删除。
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(!if_any(.cols =  c(description,organism_name,name,allele_name),.fns = `==` , ""))
### 将peptide长度添加进去
Tcell_v3_clean <- Tcell_v3_clean %>% 
  mutate(pep_length = stringr::str_length(description))
### filter bigger than 50 
Tcell_v3_clean <- Tcell_v3_clean %>% 
  filter(pep_length <= 50)
### add count
Tcell_v3_clean <- Tcell_v3_clean %>% group_by(description) %>% 
  add_count()
### 抽取出人来源的term
# peptide 具有重复性
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


############### 2 final matrix proprocess ########
#### 2.1 还需要清楚是从那个蛋白质中 call 出来的 #### 
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
#### 2.2 还需要清楚是那个亚型 #####
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
####-------------End----------------####

############ Next ##################
#### filter out finalResult.df file ##### 
finalResult.df.2 <- finalResult.df.filter[str_detect(string = finalResult.df.filter$allele_name,pattern = "\\w\\w\\w-\\w\\*\\d\\d:\\d\\d"),]
finalResult.df.3 <- finalResult.df.2 %>% 
  distinct() %>% 
  rename(wild_peptide = V1) %>% 
  mutate(allele_name = str_match(allele_name,pattern = "\\w\\w\\w-\\w\\*\\d\\d:\\d\\d"))

#### Calculate and load NetMHCpan result ####
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

###########################################################
#### reformat netMHCpan result ##### 
mutate.input <- read.delim(file = paste0(output_path,"pos.mutate.peptide.out.1"),sep = "\n")
wild.input <- read.delim(file = paste0(output_path,"pos.wild.peptide.out.1"),sep = "\n")
fraction_reformat <- function(fraction.mutate){
  fraction.mutate.1 <- as.data.frame(fraction.mutate[!str_detect(string = fraction.mutate[,1],pattern = "-----"),])
  fraction.mutate.2 <- str_split_fixed(string = fraction.mutate.1[,1],pattern = "\\s+",n = Inf)
  # fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
  #   mutate(V1 = rep(protein.mutate %>%  pull(V1),times = protein.times))
  fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
    select(!V1) 
  # colnames(fraction.mutate.3)[1:14] <- c("Pos","MHC","Peptide",
  #                                        "Core","Of","Gp","G1",
  #                                        "Ip","I1","Icore","Identity",
  #                                        "Score_EL","%Rank_EL","BindLevel")  
  colnames(fraction.mutate.3) <- NULL
  rownames(fraction.mutate.3) <- NULL
  #fraction.mutate.3 <- fraction.mutate.3[,-1]
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

# /mnt/data/meng/software/netMHCpan-4.1/netMHCpan -BA  -p pos.mutate.peptide > pos.mutate.peptide.out
#### join NetMHCpan result with final_result ####
#### filter peptides by nm50  ####