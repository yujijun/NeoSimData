# NeoSimData
This is a repository to generate simulated neopeptides datasets.

To construct the positive dataset, we initially collected tcell_full_v3.zip dataset from IEDB dataset, then filtered peptides which were reported to exhibit positive T-cell response in humans. To find epitopes that can be generated by hypothetical somatic mutations in human genome, we compared the
sequences of the epitopes with those of the 20198 reviewed Swiss-
Prot human proteins [19] to epitopes whose sequences
are differed from the best-match by up to 2 amino acids; the best-match
protein was further used as a corresponding wild-type.

For negative dataset, we initially collected 22 245 variants from common
(minor allele frequency [MAF]0.05) non-synonymous single nucleotide
polymorphisms (SNPs) from dbSNP v.141 with the assumption
that widespread peptide variants would not lead to an immunogenic response.
The HLA allele for negative dataset was randomly selected from
HLA alleles of the positive set. 

For all the obtained peptides datasets, it is recommended to use the netMHCpan algorithm to calculate the affinity, and use the affinity value for further screening.

This method has been proven effective and can be used in machine learning algorithms for Determination of immunogenic neo-antigens[], while our method was an update version, user can simulate the desired HLA alleles according to their own needs and use netMHCpan for the second round of screening.In this tutorial,we use HLA-A*02:01 as an example.

# Requirement 
Need to install netMHCpan and change 
# Usage:
1. Download needed datasets
2. running all dataset.

Reference:


