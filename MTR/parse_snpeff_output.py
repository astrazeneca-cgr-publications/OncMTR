import os, sys
import re
import ntpath


def process_single_annotation_entry(ann, chrom, position):

	ann_vals = ann.split('|')
	annotations = [ ann_vals[1] ]
	gene_name = ann_vals[3]
	gene_id = ann_vals[4]
	transcript_id = ann_vals[6]
	transcript_biotype = ann_vals[7]
	hgvs_c = ann_vals[9]
	err_warn_msg = ann_vals[15]

	# remove version from transcript id
	transcript_id = re.sub(r'\..*', '', transcript_id)

	
	# Check for Error messages
	if 'ERROR' in err_warn_msg:
		print(err_warn_msg, ' - Skipping entry...')

	# Parse only synonymous and missense variants
	if len(set(valid_annot) & set(annotations)) == 0:
		return	

	# Extract position of mutation in the full coding reference sequence of this particular transcript
	if not hgvs_c.startswith('c.'):
		print("not hgvs_c.startswith('c.')")
		sys.exit(ann)

	coding_position = re.sub(r'c\.', '', hgvs_c)
	coding_position = re.sub(r'[AGCT].*', '', coding_position)

	new_line = transcript_id + '\t' + coding_position + '\n'
	

	if len(set(missense_annotations) & set(annotations)) >= 1:
		mis_fh.write(new_line)

	elif len(set(synonymous_annotations) & set(annotations)) >= 1:
		syn_fh.write(new_line)





if __name__ == '__main__':
		
	snpeff_results_file = sys.argv[1]

	snpeff_results_file_basename = ntpath.basename(snpeff_results_file)
	print("Input file basename:", snpeff_results_file_basename)

	snpeff_processed_out_dir = "../out/snpeff_processed_out"
	if not os.path.exists(snpeff_processed_out_dir):
		os.makedirs(snpeff_processed_out_dir)


	missense_annotations = ['missense_variant', 'missense_variant&splice_region_variant']
	synonymous_annotations = ['synonymous_variant', 'stop_retained_variant', 'splice_region_variant&stop_retained_variant', 'splice_region_variant&synonymous_variant']
	#missense_annotations = ['missense_variant', 'initiator_codon_variant']
	#synonymous_annotations = ['synonymous_variant', 'start_retained', 'stop_retained_variant']
	valid_annot = missense_annotations + synonymous_annotations
	

	mis_out_file = snpeff_processed_out_dir + '/' + snpeff_results_file_basename + '.missense.exon_coords'
	syn_out_file = snpeff_processed_out_dir + '/' + snpeff_results_file_basename + '.synonymous.exon_coords'
	
	mis_fh = open(mis_out_file, 'w')
	syn_fh = open(syn_out_file, 'w')


	cnt = 1
	with open(snpeff_results_file) as fh:
		for line in fh:
			if line.startswith('#'):
				continue

			line = line.rstrip()

			vals = line.split('\t')
			chrom, position = vals[0], vals[1]


			ann_field = vals[7].split(';')[0] # there may be some extra annotations ('LOF' or 'NMD') ";"-delimited, that we can discard
			ann_field_entries = ann_field.split(',')

			for ann in ann_field_entries:
				process_single_annotation_entry(ann, chrom, position)

			if cnt % 100000 == 0:
				print(cnt)
			cnt += 1


	mis_fh.close()
	syn_fh.close()
