# Ancestry-specific MTR

Ancestry-specific MTR scores are calculated in the same way as regular MTR scores, but as an input they require lists of ancestry-specific genetic variants.

## Requirements

* Input: VCF file with variants present in a given ancestry (cohort)
* SnpEff 4.3t with GRCh38.86
* Python packages in `requirements.txt` (you can install them at once using `pip install -r requirements.txt`)

## Example usage

In this example, as an input we will use a VCF file with first 1000 variants present in the `NFE_440k` cohort (`cohort="NFE_440k"`): `data/snpeff_input/${cohort}.snpeff_input.All_chr.vcf`.

The whole pipeline consist of three steps: running SnpEff to annotate the variants, parsing the results, and calculating ancestry-specific MTR scores.

### Run SnpEff

SnpEff version 4.3t is used to annotate variants using GRCh38.86 database. You can download it from [here](https://pcingola.github.io/SnpEff/download/). Java 1.8.0 is required to run SnpEff.

The command to run SnpEff:

```
cohort="NFE_440k"
java -Xmx4g -jar ${path_to_snpEff_jar} -noStats -v -t GRCh38.86 data/snpeff_input/${cohort}.snpeff_input.All_chr.vcf > out/snpeff_results/${cohort}_exons_mutations.snpeff_annotation.pass.vcf
```

### Parse SnpEff results

```
cd MTR
python -u parse_snpeff_output.py ../out/snpeff_results/${cohort}_exons_mutations.snpeff_annotation.pass.vcf
python collate_hits_per_transcript.py ${cohort}
```

### Calculate MTR scores

Here, we will calculate MTR scores for a window of 31 codons. Lists of all possible missense and synonymous variants are located at `out/snpeff_processed_out` and are stored in the `pickle` format. For the purpose of this exercise, lists only for seven transcripts are given there. Please contact us if you need the full lists (as they are over 200 MB each). 

```
cd MTR
python -u calculate_mtr_scores.py ${cohort}_exons_mutations.snpeff_annotation.pass $cohort 31
```

The ancestry-specific MTR scores for individual transcripts will be saved in `out/MTR-${cohort}-win31`.

