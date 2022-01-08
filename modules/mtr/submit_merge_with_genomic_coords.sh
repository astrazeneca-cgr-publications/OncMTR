#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=16
#SBATCH -o out.merge_with_genomic_coords

dataset=$1	# gnomad
win=$2
ab_filter=$3   # yes or no


echo "$dataset"
python -u merge_mtr_scores_with_genomic_coords.py -d $dataset -w $win -f $ab_filter

