import sys
import re

# -- Input: exon_sequences.fa
# -- Notes: Coordinates from input sequences file are 0-based.
# 	    Converting them to 1-based for VCF output.


def get_alternative_nts(ref_base):
	
	bases = ['A', 'C', 'G', 'T']
	alt_bases = [n for n in bases if n != ref_base]

	return alt_bases



def create_artef_entries_for_vcf(header, seq, out_vcf):

	chrom, start, end = re.split(':|-', header)
	#print(chrom, start, end)	

	# convert 0-based indexing to 1-based for VCF
	start = int(start) + 1
	end = int(end)

	idx = 0	
	for coord in range(start, end+1):

		ref_base = seq[idx]	
		idx += 1

		if ref_base == 'N':
			print("[Warning] Detected 'N' base in exon.") 
			continue

		alt_bases = get_alternative_nts(ref_base)

		for alt in alt_bases:
			tmp_out_line = chrom + '\t' + str(coord) + '\t' + '.' + '\t' + ref_base + '\t' + alt + '\t.\t.\t.\n'
			out_vcf.write(tmp_out_line)
			#print(tmp_out_line)
	
	


if __name__ == '__main__':

	data_dir = '../../data/ensembl/'

	input_fasta = data_dir + '/exon_sequences.fa' 
	output_vcf = data_dir + '/snpeff_input.all_exons_mutations.vcf'

	# ----------- Canonical (longest) transcripts only -----------
	#input_fasta = 'exon_sequences.canonical.fa' 
	#output_vcf = 'snpeff_input.all_exons_mutations.canonical.vcf'
	# ____________________________________________________________

	out_vcf = open(output_vcf, 'w')
	vcf_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	out_vcf.write(vcf_header)

	header = ''
	seq = ''

	cnt = 0
	with open(input_fasta, 'r') as seqs_fh:
		for line in seqs_fh:
			line = line.rstrip()

			if line.startswith('>'):
				header = line
				header = header.replace('>', '')
			else:
				seq = line
				create_artef_entries_for_vcf(header, seq, out_vcf)

			cnt += 1
			if cnt % 10000 == 0:
				print(cnt)

	out_vcf.close()
