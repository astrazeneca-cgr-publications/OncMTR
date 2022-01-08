## Workflow
```
- cd data/clinvar/vcf_GRCh38
  process_clinvar.py  


- cd filtered-clinvar-files
  prepare_snpeff_input.py 
	Input: full_clinvar_df.csv
	Output: full_clinvar_df.snpeff_input.vcf


- cd modules/snpeff_annotation
	./submit.snpEff.all_clinvar.sh
	Input: full_clinvar_df.snpeff_input.vcf
	Output: full_clinvar_df.snpeff_output.vcf


- cd data/clinvar/vcf_GRCh38/filtered-clinvar-files
  parse_clinvar_snpeff_output.py  
	Input: full_clinvar_df.snpeff_output.vcf
	Output: clinvar.All_annotations.tsv


  subset_clinvar_all_annotations.py
```





- Create ClinVar file without somatic variants
```
cat clinvar.vcf | grep -v '#' | grep -v ORIGIN=2 | grep -i -e CLNSIG=Pathogenic -e CLNSIG=Likely_pathogenic -e Uncertain_significance | sed -E 's/(.*\t.*\t.*\t.*\t.*\t.*\t.*\t).*CLNSIG=(.*)\.*/\1\2/' | sed 's/;.*//' > clinvar.pathogenic_subset.no_somatic.vcf
```

- Capture germline-only variants
```
cat clinvar.vcf | grep -e "ORIGIN=1$" | wc -l   # ORIGIN is the last field
cat clinvar.vcf | grep -e "ORIGIN=1;" | wc -l	# ORIGIN is followed by another field
```

