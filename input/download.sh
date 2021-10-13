# This is a script for downloading of input datasets.
### Download postive datasets 
# uniprot-proteome_UP000005640.fasta
wget https://www.uniprot.org/uniprot/?include=false&format=fasta&compress=yes&force=true&query=proteome:UP000005640
# tcell_full_v3.csv 
wget https://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip


### Download negative datasets
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
  
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz

wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz

### decompress gzip 
gzip -d *.gz
unar *.zip

#### filtered common vcf file #####
grep "NSM" common_all_20180418.vcf | grep "G5A" > common_all_20180418_snp.recode.maf5.nsm.vcf


