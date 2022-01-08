#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH -o all_clinvar.snpeff.out
##SBATCH --mem-per-cpu=4G
##SBATCH -n 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

ml Java

snpeff_dir="../../utils/snpEff"
clinvar_dir="../../data/clinvar/vcf_GRCh38/filtered-clinvar-files"


#java -Xmx4g -jar ../../utils/snpEff/snpEff.jar -noStats -v -t GRCh38.86 $clinvar_dir/full_clinvar_df.snpeff_input.vcf > $clinvar_dir/full_clinvar_df.snpeff_output.vcf

java -Xmx4g -jar ../../utils/snpEff/snpEff.jar -noStats -v -t GRCh38.92 -c $snpeff_dir/snpEff.config $clinvar_dir/full_clinvar_df.snpeff_input.vcf > $clinvar_dir/full_clinvar_df.snpeff_output.vcf
