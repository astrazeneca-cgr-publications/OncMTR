#!/bin/bash
#SBATCH -o out.get_pct_scores
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH -t 24:00:0

dataset=$1
win=$2
ab_filter=$3

echo "$dataset; $win; $ab_filter"
python -u get_variant_level_percentile_mtr_scores.py -d $dataset -w $win -f $ab_filter
