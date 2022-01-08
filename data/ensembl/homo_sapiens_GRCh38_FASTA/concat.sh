#!/bin/bash 

concat_file=Homo_sapiens.GRCh38.dna.all_chromosomes.fa
rm -f $concat_file


for i in `seq 1 22`;
do
	cat Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz >> $concat_file
done

i=X
cat Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz >> $concat_file

i=Y
cat Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz >> $concat_file


mv $concat_file ${concat_file}.gz
