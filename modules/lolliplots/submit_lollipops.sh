#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -o lollipops.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

python create_lollipop_for_all_enst_ids.py
