#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH -o parse_snpeff_output.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=3


snpeff_results=$( grep "^snpeff_res" ../../.conf | awk 'BEGIN {FS="\t"} { print $2 }' )

# -- GnomAD
python -u parse_snpeff_output.py $snpeff_results/gnomad_exons_mutations.snpeff_annotation.pass.vcf &
# (AB-filtered)
python -u parse_snpeff_output.py $snpeff_results/gnomad_exons_mutations.snpeff_annotation.pass.ab_filtered.vcf &



# -- All possible mutations
python -u parse_snpeff_output.py $snpeff_results/all_exons_mutations.snpeff_annotation.vcf &



wait
