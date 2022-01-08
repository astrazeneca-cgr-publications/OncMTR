# Analyis Steps

## Compile OncMTR from MTR and MTR-AB results
- compile_oncmtr.py

```
# e.g.
python compile_oncmtr.py -d gnomad -w 31
```



## Annotate OncMTR signals with variants (e.g. clinvar/CGC) and create plots
- annotate_oncmtr_signals.py

```
# e.g.
python annotate_oncmtr_signals.py -d gnomad -w 31

# or with sbatch
sbatch submit_annotate_oncmtr_signals.sh
```



## Get genomic coordinates for each OncMTR score:
- merge_oncmtr_scores_with_genomic_coords.py

```
# e.g.
python merge_oncmtr_scores_with_genomic_coords.py -d gnomad -w 31

# or with sbatch
sbatch submit_merge_with_genomic_coords.sh
```

## Get variant level percentile scores
```
sbatch submit_get_pct_scores.sh
```




## Get OncMTR distribution (based on various metrics) in CGC Tier 1/2 vs rest of exome
```
python plot_cgc_geneset_distributions.py
```
