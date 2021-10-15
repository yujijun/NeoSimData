################ Description ##### 
################ library and hyperparameter ####
output_path <- "./output/neg"
input_path <- "./output/neg"
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

#### Running netMHCpan #####
# running Neg_netMHCpanRunning.sh


