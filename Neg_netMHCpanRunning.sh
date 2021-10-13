cd ./output/neg
### Running netMHCpan of negative fraction to get affinity prediction result.
# Attention: Users need to change `/mnt/data/meng/software/netMHCpan-4.1/netMHCpan` with their own netMHCpan software pathway
/mnt/data/meng/software/netMHCpan-4.1/netMHCpan  fraction.wild.fsa > fraction.wild.fsa.out
/mnt/data/meng/software/netMHCpan-4.1/netMHCpan  fraction.mutate.fsa > fraction.mutate.fsa.out
### preprocess about output file####
grep -v "^#" fraction.mutate.fsa.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > fraction.mutate.fsa.out.1
grep -v "^#" fraction.wild.fsa.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > fraction.wild.fsa.out.1

cd ../..