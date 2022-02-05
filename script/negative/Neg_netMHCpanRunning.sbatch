#!/bin/bash 
#SBATCH -J Neg_netMHCpanRunning
#SBATCH -p xhacnormalb
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -o out.%j
#SBATCH -e err.%j

source activate base 
python /public/home/wumeng01/NeoantigenML/QBRC-Neoantigen-Pipeline/code/mhc_i/src/predict_binding.py netmhcpan_ba HLA-A*02:01 9 /public/home/yujijun01/NeoStim/fraction.wild.fsa > /public/home/yujijun01/NeoStim/fraction.wild.fsa.out

python /public/home/wumeng01/NeoantigenML/QBRC-Neoantigen-Pipeline/code/mhc_i/src/predict_binding.py netmhcpan_ba HLA-A*02:01 9 /public/home/yujijun01/NeoStim/fraction.mutate.fsa > /public/home/yujijun01/NeoStim/fraction.mutate.fsa.out

### preprocess about output file####
#grep -v "^#" fraction.mutate.fsa.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > fraction.mutate.fsa.out.1
#grep -v "^#" fraction.wild.fsa.out |  grep -v "^ Pos" | grep -v "^HLA" | grep -v "^Protein" > fraction.wild.fsa.out.1
