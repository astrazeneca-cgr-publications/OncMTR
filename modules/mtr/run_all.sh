#!/bin/bash

# - gnomAD
sbatch -o gnomad31_mtr_calc.out submit_mtr_calc.sh -d gnomad -w 31 
sbatch -o gnomad31_mtr_calc-ab_filtered.out submit_mtr_calc.sh -d gnomad -w 31 -f yes

