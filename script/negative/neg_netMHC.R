################ Description ##### 
# this function is used for df5.5 -> fraction -> netMHCpan running -> reformat netMHC result
################ library and hyperparameter ####
library(seqinr)
library(tidyverse)
output_path <- "output/neg"
input_path <- "output/neg"
################ Convert df5.5 into fraction.mutate / wild fsa ##### 
load(paste0(input_path,"ncbi_final_df.5.5.RData"))
df.5.6 <- df.5.5 %>% 
  select(seq_proteins_fraction,mutateseq1_proteins_fraction) %>% 
  distinct() 
df.5.6$seq_proteins_fraction <- unlist(df.5.6$seq_proteins_fraction)
df.5.6$mutateseq1_proteins_fraction <- unlist(df.5.6$mutateseq1_proteins_fraction)
df.5.6 <- df.5.6 %>% 
  filter(!(seq_proteins_fraction == mutateseq1_proteins_fraction))
fraction.choose.wild <- df.5.6$seq_proteins_fraction
write.fasta(sequences = as.list(fraction.choose.wild),
            names = fraction.choose.wild,
            file.out =  paste0(output_path,"fraction.wild.fsa"))
fraction.choose.mutate <- df.5.6$mutateseq1_proteins_fraction
write.fasta(sequences = as.list(fraction.choose.mutate),
            names = fraction.choose.mutate,
            file.out =  paste0(output_path,"fraction.mutate.fsa"))

################ Running netMHCpan #####
# running Neg_netMHCpanRunning.sh

################ fraction reformat function  ##### 
# fraction.mutate <- read.delim(file = "./output/test/fraction.fsa.out.test.1",sep = "\n")
# protein.mutate <- read.delim(file = "output/test/fraction.tsv",header = F)
fraction_reformat <- function(fraction.mutate){
  fraction.mutate.1 <- as.data.frame(fraction.mutate[!str_detect(string = fraction.mutate[,1],pattern = "-----"),])
  fraction.mutate.2 <- str_split_fixed(string = fraction.mutate.1[,1],pattern = "\\s+",n = Inf)
  # fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
  #   mutate(V1 = rep(protein.mutate %>%  pull(V1),times = protein.times))
  fraction.mutate.3 <- as.data.frame(fraction.mutate.2) %>% 
    select(!V1)
  colnames(fraction.mutate.3)[1:14] <- c("Pos","MHC","Peptide",
                              "Core","Of","Gp","G1",
                              "Ip","I1","Icore","Identity",
                              "Score_EL","%Rank_EL","BindLevel")  
  return(fraction.mutate.3)
}
#result <- fraction_reformat(fraction.mutate,protein.mutate,26)
#save(fraction.mutate.final,file = paste0("output/test/","fraction.mutate.final.RData"))


################ Running ##### 
mutate.input <- read.delim(file = paste0(output_path,"fraction.mutate.fsa.out.1"),sep = "\n")
#mutate.fraction <- read.delim(file = paste0(output_path,"fraction.mutate.list"),header = F)
#mutate.number <- read.delim(file = paste0(output_path,"fraction.mutate.number"),header = F)
#mutate.number.1 <- as.numeric(str_split_fixed(mutate.number$V1,pattern = "\\s+",n = Inf)[,18])
mutate.result <- fraction_reformat(mutate.input)
write.table(mutate.result,file = paste0(output_path,"fraction.mutate.fsa.out.1.Rformat.tsv"),
            sep = "\t",row.names = F,col.names = T,quote = F)

wild.input <- read.delim(file = paste0(output_path,"fraction.wild.fsa.out.1"),sep = "\n")
# wild.fraction <- read.delim(file = paste0(output_path,"fraction.wild.list"),header = F)
# wild.number <- read.delim(file = paste0(output_path,"fraction.wild.number"),header = F)
# wild.number.1 <- as.numeric(str_split_fixed(wild.number$V1,pattern = "\\s+",n = Inf)[,18])
wild.result <- fraction_reformat(wild.input)
write.table(wild.result,file = paste0(output_path,"fraction.wild.fsa.out.1.Rformat.tsv"),
            sep = "\t",row.names = F,col.names = T,quote = F)


