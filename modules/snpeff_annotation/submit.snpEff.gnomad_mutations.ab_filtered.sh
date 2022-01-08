#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -o gnomad.ab_filtered.snpeff.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16

#ml Java
module load Java/1.8.0_221

input_dir=$( grep "^gnomad_dir" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )
snpeff_dir=$( grep "^snpeff_dir" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )
out_dir=$( grep "^snpeff_res" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )

mkdir -p $out_dir


# PASS variants only
#java -Xmx4g -jar $snpeff_dir/snpEff.jar -noStats -v -t GRCh38.86 $input_dir/snpeff_input.gnomad_exons_mutations.pass.ab_filtered.vcf > $out_dir/gnomad_exons_mutations.snpeff_annotation.pass.ab_filtered.vcf

# Ensembl v. GRCh38.92
java -Xmx4g -jar $snpeff_dir/snpEff.jar -noStats -v -t GRCh38.92 -c $snpeff_dir/snpEff.config $input_dir/snpeff_input.gnomad_exons_mutations.pass.ab_filtered.vcf > $out_dir/gnomad_exons_mutations.snpeff_annotation.pass.ab_filtered.vcf
