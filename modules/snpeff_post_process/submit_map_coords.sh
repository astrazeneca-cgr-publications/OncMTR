#!/bin/bash
#SBATCH -o out.map_coords
#SBATCH -e err.map_coords
#SBATCH -t 8:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=64G

python -u  map_genomic_to_hgvs_from_snpeff_output.py
