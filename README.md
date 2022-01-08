# OncMTR
Oncology-specific MTR score calculated from the gnomAD (v2) cohort:
- [MTR](#mtr)
- [OncMTR](#oncmtr)


<br><br>


MTR
===

### Calculate original and AB_MEDIAN-filtered versions of MTR

<br>

- **Pre-generated data**: Snpeff annotations on all exons for all possible mutations
```
out/snpeff_results/all_exons_mutations.snpeff_annotation.vcf
```

- Step 1. 
```
cd modules/snpeff_annotation
sbatch ./submit.snpEff.gnomad_mutations.sh

# Note - Expected input VCF (e.g. from gnomAD)]:
data/gnomad_GRCh38/snpeff_input.gnomad_exons_mutations.pass.vcf
```


- Step 2.
```
cd modules/snpeff_post_process
sbatch ./submit_parse_snpeff_output.sh 
```


- Step 3. (same dir as in step 2)
```
python collate_hits_per_transcript.py  
```


- Step 4. 
```
cd modules/mtr

./run_all.sh
```

or step by step

```
python -u calculate_mtr_scores.py -d gnomad -w 31 -s mtr -f  no
python -u calculate_mtr_scores.py -d gnomad -w 31 -s mtr-ab -f  yes

# submitted via
sbatch submit_mtr_calc.sh
```


<br>


OncMTR
===

### Compile OncMTR from MTR and MTR-AB results
- compile_oncmtr.py

```
# e.g.
python compile_oncmtr.py -d gnomad -w 31
```



### Annotate OncMTR signals with variants (e.g. clinvar/CGC) and create plots
- annotate_oncmtr_signals.py

```
# e.g.
python annotate_oncmtr_signals.py -d gnomad -w 31

# or with sbatch
sbatch submit_annotate_oncmtr_signals.sh
```



### Get genomic coordinates for each OncMTR score:
- merge_oncmtr_scores_with_genomic_coords.py

```
# e.g.
python merge_oncmtr_scores_with_genomic_coords.py -d gnomad -w 31

# or with sbatch
sbatch submit_merge_with_genomic_coords.sh
```

### Get variant level percentile scores
```
sbatch submit_get_pct_scores.sh
```




### Get OncMTR distribution (based on various metrics) in CGC Tier 1/2 vs rest of exome
```
python plot_cgc_geneset_distributions.py
```
