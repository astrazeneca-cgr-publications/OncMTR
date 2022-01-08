#!/bin/bash

# - gnomAD
sbatch -o out.gnomad_31.get_pct submit_get_pct_scores.sh gnomad 31 no
sbatch -o out.gnomad_31-ab_filtered.get_pct submit_get_pct_scores.sh gnomad 31 yes