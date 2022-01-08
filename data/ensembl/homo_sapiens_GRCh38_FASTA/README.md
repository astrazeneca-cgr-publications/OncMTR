```
# download fasta files for each chromosome from Ensembl (human build GRCh38)
./download.sh

## [IMPORTANT]: replace all header with >chr1, >chr2, etc. (or >1, >2, etc.)

# concatenate chromosome fasta files into a single one
./concat.sh

# unzip concatenated fa.gz
gunzip Homo_sapiens.GRCh38.dna.all_chromosomes.fa.gz

# create sequence dictionary (.dict) file
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.all_chromosomes.fa

# create .fai index file
samtools faidx Homo_sapiens.GRCh38.dna.all_chromosomes.fa
```

# Retrieve sequences from indexed human genome
``` 
faToTwoBit Homo_sapiens.GRCh38.dna.all_chromosomes.fa hsa38.2bit   

twoBitInfo hsa38.2bit stdout | sort -k2rn > hsa38.chrom.sizes   

# Test: twoBitToFa hsa38.2bit:1:1-10 /dev/stdout 
```
