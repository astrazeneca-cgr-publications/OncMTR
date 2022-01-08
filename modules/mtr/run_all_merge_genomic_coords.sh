#!/bin/bash

# - gnomAD
sbatch -o out.gnomad_31.merge_with_genomic_coords submit_merge_with_genomic_coords.sh gnomad 31 no
sbatch -o out.gnomad_31-ab_filtered.merge_with_genomic_coords submit_merge_with_genomic_coords.sh gnomad 31 yes
