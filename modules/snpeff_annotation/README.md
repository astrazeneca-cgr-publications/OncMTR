## Run snpEff to annotate different subsets of variants:

- All possible exon variants (artificially constructed based on Ensembl annotation)
```
sbatch ./submit.snpEff.all_mutations.sh
```

- All GnomAD PASS variants
```
sbatch ./submit.snpEff.gnomad_mutations.sh
```

- GnomAD PASS variants with AB_MEDIAN > 0.35
```
sbatch ./submit.snpEff.gnomad_mutations.ab_filtered.sh
```

- ClinVar pathogenic variants
```
sbatch ./submit.snpEff.clinvar_pathogenic.sh
```
