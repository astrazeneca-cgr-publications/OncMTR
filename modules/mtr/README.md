## Calculate MTR scores
```
python -u calculate_mtr_scores.py -d [dataset] -w [win_size] -s [out_suffix (optional)]
```
or via
```
sbatch submit_mtr_calc.sh
```



## Get variant-level MTR scores
```
sbatch submit_merge_with_genomic_coords.sh
```
