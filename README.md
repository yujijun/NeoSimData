# NeoSimData

## Description

This is a repository to generate simulated neopeptides.

To construct the positive dataset, we initially collected tcell_full_v3.zip dataset from IEDB dataset[1], then filtered peptides which were reported to exhibit positive T-cell response in humans. To find epitopes that can be generated by hypothetical somatic mutations in human genome, we compared the
sequences of the epitopes with those of the 20198 reviewed Swiss-
Prot human proteins [2] to epitopes whose sequences
are differed from the best-match by up to 2 amino acids; the best-match
protein was further used as a corresponding wild-type.

For negative dataset, we initially collected 22 245 variants from common
(minor allele frequency [MAF] > 0.05) non-synonymous single nucleotide
polymorphisms (SNPs) from dbSNP v.141 with the assumption
that widespread peptide variants would not lead to an immunogenic response.
The HLA allele for negative dataset was randomly selected from
HLA alleles of the positive set. 

For all the obtained peptides datasets, it is recommended to use the netMHCpan algorithm[3] to calculate the affinity, and use the affinity value for further screening.

This method has been proven effective and can be used in machine learning algorithms for Determination of immunogenic neoantigens[4], while our method was an update version, user can simulate the desired HLA alleles according to their own needs and use netMHCpan for the second round of screening. The results of netMHCpan show that majority of the modeled positive peptides have lower binding affinity(nM), and the negative peptides have higher binding affinity(nM), which further proves the effectiveness of this method.

In this tutorial,we use `HLA-A*02:01` as an example.

## Requirement 

### R Environment

R version 4.0.3 (2020-10-10)

attached base packages:

[1] stats      
[2] graphics     
[3] grDevices     
[4] utils        
[5] datasets    
[6] methods  
[7] base     

other attached packages:

 [1] seqinr_4.2-8       
 [2] Biostrings_2.58.0  
 [3] XVector_0.30.0     
 [4] IRanges_2.24.1     
 [5] S4Vectors_0.28.1   
 [6] BiocGenerics_0.36.1     
 [7] forcats_0.5.1        
 [8] stringr_1.4.0      
 [9] dplyr_1.0.7        
[10] purrr_0.3.4        
[11] readr_2.0.2        
[12] tidyr_1.1.4        
[13] tibble_3.1.5       
[14] ggplot2_3.3.5      
[15] tidyverse_1.3.1  

### MHC peptides binding prediction software

[1] NetMHCpan-4.1   
[2] NetMHCIIpan-4.0

# Usage:
```
git clone git@github.com:yujijun/NeoSimData.git
cd /NeoSimData 
cd /input 
./download.sh #download all needed input files.

#### negative neopeptides dataset:
./script/neg_main.R  ## Generate protein fractions of wild and mutation 
./script/Neg_netMHCpanRunning.sh ## Running NetMHCpan for MHC peptide affinity information. attention: Change into your own NetMHCpan path. 
./script/neg_format.R ## Format and filter NetMHCpan Running result, put all peptides togethers. 

#### postive neopeptides dataset:
./script/pos_main.R  ## Generate protein fractions of wild and mutation 
./script/Pos_netMHCpanRunning.sh ## Running NetMHCpan for MHC peptide affinity information. attention: Change into your own NetMHCpan path. 
./script/pos_format.R ## Format and filter NetMHCpan Running result, put all peptides togethers. 
```

## Output example

All of the result was after filtered from NetMHCpan result.

![Postive result example](https://github.com/yujijun/NeoSimData/blob/main/output/Positive_example.png)

![Negative result example](https://github.com/yujijun/NeoSimData/blob/main/output/Negative_example.png)

## Reference:

[1] Vita R, Mahajan S, Overton JA, Dhanda SK, Martini S, Cantrell JR, Wheeler DK, Sette A, Peters B. The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res. 2018 Oct 24. doi: 10.1093/nar/gky1006. [Epub ahead of print] PubMed PMID: 30357391.

[2] The Universal Protein Resource (UniProt) Nucleic Acides Res.2008 Jan; 36(Dataset issue):D190-D195.

[3] NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379

[4] Neopepsee: accurate genome-level prediction of
neoantigens by harnessing sequence and amino acid
immunogenicity information. Annals of Oncology 29: 1030–1036, 2018
doi:10.1093/annonc/mdy022


