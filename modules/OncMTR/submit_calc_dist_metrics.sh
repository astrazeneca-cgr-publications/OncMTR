#!/bin/bash
#SBATCH -o out.calc_dist_metrics
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -t 24:00:0

python -u calc_oncmtr_dist_metrics.py -d gnomad -w 31 -s full_metrics
