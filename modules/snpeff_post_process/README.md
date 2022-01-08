# Parse snpEff annotation output for various input sets 
```
sbatch ./submit_parse_snpeff_output.sh   # (~ 3hr)
```


# Create dictionary with all hits per transcript - to be used for MTR calculations
python collate_hits_per_transcript.py
# srun --mem-per-cpu=16G --cpus-per-task=1 python -u collate_hits_per_transcript.py



# ==== Genomic coordinates mapping ====
# - Map genomic to protein coordinates (Ensembl v92)
sbatch ./submit_map_coords.sh (~4hr)

# - Split files by chr
python split_coords_by_chr.py
