#!/bin/bash 
#SBATCH --time=24:00:00 
#SBATCH -o clinvar_pathogenic.snpeff.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16  

ml Java

input_dir="../../clinvar/vcf_GRCh38"
snpeff_dir=$( grep "^snpeff_dir" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )
out_dir=$( grep "^snpeff_res" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )

mkdir -p $out_dir

#java -Xmx4g -jar $snpeff_dir/snpEff.jar -noStats -v -t GRCh38.86 $input_dir/clinvar.pathogenic_subset.vcf > $out_dir/clinvar.pathogenic_subset.snpeff_annotation.vcf

java -Xmx4g -jar $snpeff_dir/snpEff.jar -noStats -v -t GRCh38.92 -c $snpeff_dir/snpEff.config $input_dir/clinvar.pathogenic_subset.vcf > $out_dir/clinvar.pathogenic_subset.snpeff_annotation.vcf
