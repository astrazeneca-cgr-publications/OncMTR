# Data preparation and pre-processing

`data_dir = ../../data/ensembl`

1. Extract sequences for each exon:
#### (call keep_canonical_transcripts.py within it and with input/output files accordingly if looking into canonical transcripts only)
```
./get_seqs_for_exons.sh   
```

- Output:
```
$data_dir/exon_sequences.fa
```

2. Create VCF for snpEFF.
python prepare_vcf_for_snpeff.py


3. Get dictionaries with ENST->ENSG and ENST->Uniprot mappings:
```
python map_enst_to_ensgenes_and_uniprot.py
```
