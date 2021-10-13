cd ./output/pos
### Running netMHCpan of postive peptides to get affinity prediction result.
/mnt/data/meng/software/netMHCpan-4.1/netMHCpan -BA  -p pos.mutate.peptide > pos.mutate.peptide.out
/mnt/data/meng/software/netMHCpan-4.1/netMHCpan -BA  -p pos.wild.peptide > pos.wild.peptide.out

###
### preprocess about output file####
grep -v "^#" pos.mutate.peptide.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > pos.mutate.peptide.out.1
grep -v "^#" fraction.wild.fsa.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > pos.wild.peptide.out.1

cd .../..