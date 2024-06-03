#module load Java/1.8.0_221

java -Xmx4g -jar ../utils/snpEff.jar -noStats -v -t GRCh38.86 ../data/${cohort}.snpeff_input.All_chr.vcf > ../out/snpeff_results/${cohort}_exons_mutations.snpeff_annotation.pass.vcf