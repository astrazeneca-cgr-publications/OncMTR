#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=16
#SBATCH -o out.gnomad31-OncMTR.merge_with_genomic_coords


python -u merge_oncmtr_scores_with_genomic_coords.py -d gnomad -w 31 