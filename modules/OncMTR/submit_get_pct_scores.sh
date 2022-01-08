#!/bin/bash
#SBATCH -o out.get_pct_scores
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH -t 24:00:0

python -u get_variant_level_percentile_oncmtr_scores.py -d gnomad -w 31
