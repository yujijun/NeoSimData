#### Description ######
# Created negative neoantigen datasets by SNP information from NCBI
# Create Time ：2021-08-22 10:19
# Update Time: Thu Jan 13 11:46:33 2022

#### hyperparameter and library #### 
input_path <- "./input/"
output_path <- "./output/neg/"
library(tidyverse)
library(Biostrings)
library(seqinr)
source("./script/negative/Commonfunction.R")
time1 <- Sys.time()
testGeneNumber <- 3 ### All genes need so much time,we can just test some genes.

####-------------------------------Running----------------------------------####
#### 1 input vcf\GFF\mRNArefseq and data preprocessing ####
ncbi_vcf <- read.table(file = stringr::str_c(input_path,"common_all_20180418_snp.recode.maf5.nsm.vcf"))
colnames(ncbi_vcf) <- c("CHROM",	"POS",	"ID",	"REF",
                        "ALT",	"QUAL",	"FILTER",	"INFO")
ncbi_GFF <- read.table(file = stringr::str_c(input_path,"GRCh38_latest_genomic.gff"),
                       nrows = -1,sep = "\t") 
colnames(ncbi_GFF) <- c("seqid","source","type",
                        "start","end","score",
                        "strand","phase","attributes")
ncbi_refmrna <- readDNAStringSet(filepath = stringr::str_c(input_path,"GRCh38_latest_rna.fna"))
ncbi_mrna_name <- names(ncbi_refmrna)
ncbi_mrna_seq <- as.character(ncbi_refmrna)
#### 1.1 preprocess of vcf  
ncbi_vcf_clean <- ncbi_vcf %>% 
  mutate(GENEINFO = stringr::str_match(pattern = "GENEINFO=.*?;",INFO)) %>% 
  tidyr::separate(col = GENEINFO,into = c("Gene1","Gene2"),sep = "\\|",remove = F) %>% 
  mutate(Gene1_v1 = stringr::str_remove(Gene1,pattern = "GENEINFO=")) %>% 
  tidyr::separate(col = Gene1_v1,into = c("gene1Name","gene1ID"),sep = ":",remove=T) %>% 
  tidyr::separate(col = Gene2, into = c("gene2Name","gene2ID"),remove = T,sep = ":") %>%   
  mutate(gene1ID = stringr::str_remove(gene1ID,pattern = ";")) %>% 
  mutate(gene2ID = stringr::str_remove(gene2ID,pattern = ";"))
#### 1.2 preprocess of GFF 
ncbi_GFF_clean <- ncbi_GFF %>% 
  mutate(ID = stringr::str_match(pattern = "ID=.*?;",attributes)) %>% 
  mutate(Parent = stringr::str_match(pattern = "Parent=.*?;",attributes)) %>% 
  mutate(Dbxref = stringr::str_match(pattern = "Dbxref=.*?;",attributes)) %>%
  mutate(Name = stringr::str_match(pattern = "Name=.*?;",attributes)) %>%
  mutate(gdkey = stringr::str_match(pattern = "gbkey=.*?;",attributes)) %>%
  mutate(gene = stringr::str_match(pattern = "gene=.*?;",attributes)) %>%
  mutate(product = stringr::str_match(pattern = "product=.*?;",attributes)) %>%
  mutate(pseudo = stringr::str_match(pattern = "pseudo=.*?;",attributes)) %>%
  mutate(transcript_id = stringr::str_match(pattern = "transcript\\_id=.*",attributes)) 
ncbi_GFF_clean <- ncbi_GFF_clean %>%   
  mutate(across(.cols = c("ID","Parent"),
                .fns = ~stringr::str_remove(string = .x,pattern = ".+?-")),
         .keep = "all") %>% 
  mutate(transcript_id =stringr::str_remove(string = transcript_id,pattern = "transcript_id=")) %>% 
  mutate(gene =stringr::str_remove(string = gene,pattern = "gene=")) %>% 
  mutate(gene = stringr::str_remove(string = gene,pattern = ";"))
#### 1.3 Preprocess of mRNArefseq 
ncbi_mrna_clean <- 
  as.data.frame(ncbi_mrna_name) %>% 
  mutate(id = str_split_fixed(ncbi_mrna_name,pattern = " ", n =2)[,1]) %>% 
  mutate(seq = ncbi_mrna_seq) %>% 
  mutate(seqLength = stringr::str_length(seq)) %>% 
  column_to_rownames(var = "id")

#### 2 Find All gene info ####
gene_name_list <- ncbi_vcf_clean %>% pull(gene1Name) %>% unique()
names(gene_name_list) <- ncbi_vcf_clean %>% pull(gene1Name) %>% unique()
df <- gene_match_fun(gene_name_list,ncbi_GFF_clean) # gene name match transcript_id
df.2 <- match_seq_fun(df, ncbi_mrna_clean) # transcript_id match seq 
#### 3 merge gene info matrix into vcf matrix #### 
df.3 <- dplyr::full_join(ncbi_vcf_clean,df.2,by = c("gene1Name" ="gene"))
#### 4 Determine whether the VCF matches the reference sequence ####
df.4.1 = Length_match_fun(df.3)
df.4.2 = pos_judge_fun(df.4.1)
df.4.3 = stringsame_judge_fun(df.4.2)
print("Running in df.4.3")
df.4.4 = df.4.3 %>% 
  filter(judge3 == T)


#### 5 Convert DNAseq into proteins and kmer peptides ####
df.5.1 = add_mutate_pos(df.4.4)
df.5.2 = make_mutated_dna(df.5.1)

time.df.5.2 = timeRecord(time1)

df.5.3 = convertDNAintoProtein.para(df.5.2) # Need more time

time.df.5.3 = timeRecord(time1)

print("Running in df.5.3")
df.5.4 = geneProteinFraction(df.5.3,max_peptide_length = 9)
df.5.5 = generateKmer(df.5.4,max_peptide_length = 9)
#save(df.4.3,file = paste0(output_path,"ncbi_final_df.4.3.RData")) # 这个文件保存需要花费很长时
save(df.5.5,file = paste0(output_path,"ncbi_final_df.5.5.RData")) 

####-------------------------------Outing-----------------------------------####
#### 6 output sequence ####
df.5.6 <- df.5.5 %>% 
  select(seq_proteins_fraction,mutateseq1_proteins_fraction) %>% 
  distinct() 
df.5.6$seq_proteins_fraction <- unlist(df.5.6$seq_proteins_fraction)
df.5.6$mutateseq1_proteins_fraction <- unlist(df.5.6$mutateseq1_proteins_fraction)
df.5.6 <- df.5.6 %>% 
  filter(!(seq_proteins_fraction == mutateseq1_proteins_fraction)) %>% 
  filter(!(str_detect(string = seq_proteins_fraction,pattern = "X"))) %>% 
  filter(!(str_detect(string = mutateseq1_proteins_fraction,pattern = "X")))

fraction.choose.wild <- df.5.6$seq_proteins_fraction
write.fasta(sequences = as.list(fraction.choose.wild),
            names = fraction.choose.wild,
            file.out =  paste0(output_path,"fraction.wild.fsa"))
fraction.choose.mutate <- df.5.6$mutateseq1_proteins_fraction
write.fasta(sequences = as.list(fraction.choose.mutate),
            names = fraction.choose.mutate,
            file.out =  paste0(output_path,"fraction.mutate.fsa"))
# Running time
time.all <- timeRecord(time1)
time.list <- list(time.df.5.2 = time.df.5.2,
                  time.df.5.3 = time.df.5.3,
                  time.all = time.all)
save(time.list,file = paste0(output_path,"time.record.RData"))

#### Reference ####
#1 format about GFF：
#https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#Refseq categories:https://en.wikipedia.org/wiki/RefSeq


