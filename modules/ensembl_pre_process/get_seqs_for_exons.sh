#!/bin/bash

# ================      Initialisation      ================
data_dir='../../data/ensembl'
seqs_dir="$data_dir/seq-extraction-factory"
rm -rf $seqs_dir; mkdir -p $seqs_dir
part_file_prefix=$seqs_dir/part_exons

merged_exons_bed="$data_dir/hsa_exons.GRCh38.92.merged.bed"
exon_2bit_intervals="$seqs_dir/hsa_exons.GRCh38.92.merged.2bit.intervals"



# =================         Functions        =================
function prepare_exon_intervals {
	echo "Retrieving exon intervals from gtf file into a BED file (converting 1-based to 0-based)..."
	zcat < $data_dir/Homo_sapiens.GRCh38.92.chr.gtf.gz | awk '{if($3 == "exon" && $1 != "MT") print $1"\t"$4-1"\t"$5"\t"$14"\t"$10"\t"$18}' | sed 's/[\";]//g' > $data_dir/hsa_exons.GRCh38.92.tsv


	# ---------- Canonical transcripts ----------
	# Keep canonical transcripts
	# -- Input: hsa_exons.GRCh38.92.tsv
	# -- Output: hsa_exons.GRCh38.92.canonical.bed         
	#	python keep_canonical_transcripts.py
	# Change input files accordingly ...
	# ___________________________________________


	# Get BED from tsv
	cat $data_dir/hsa_exons.GRCh38.92.tsv | awk '{print $1"\t"$2"\t"$3"\t"}' > $data_dir/hsa_exons.GRCh38.92.bed

	# Merge BED file with all transcripts
	cat $data_dir/hsa_exons.GRCh38.92.bed | sort -k1,1 -k2,2n | mergeBed > $merged_exons_bed
	
	cat $merged_exons_bed | awk '{print $1":"$2"-"$3}' > $exon_2bit_intervals
}


function split_exon_intervals_file_into_parts {
	echo "Splitting exon intervals file into smaller parts..."
	split -d -l 1000 $exon_2bit_intervals $part_file_prefix
}


function get_seqs_from_each_interval_file {

	human_ref_genome_2bit="$data_dir/homo_sapiens_GRCh38_FASTA/hsa38.2bit"

	cnt=1
	for f in ${part_file_prefix}*; do
		echo "Part: $cnt, $f"
		twoBitToFa $human_ref_genome_2bit -seqList=$f /dev/stdout | awk '/^>/ {printf("\n%s\n",$0); next;} { printf("%s",$0);}  END {printf("\n");}' | tail -n+2 > ${f}.fa

		cnt=$((cnt+1))
		rm $f
	done

	# cleanup 
	#rm ${part_file_prefix}*
}


function concat_seq_files {
	cat ${part_file_prefix}* > $data_dir/exon_sequences.fa 	
}


prepare_exon_intervals

split_exon_intervals_file_into_parts

get_seqs_from_each_interval_file

concat_seq_files
