import os, sys
import re


def process_single_annotation_entry(ann, chrom, position, ref, alt):

	global global_coord_map

	ann_vals = ann.split('|')
	try:
		annotations = [ ann_vals[1] ]
		gene_name = ann_vals[3]
		gene_id = ann_vals[4]
		transcript_id = ann_vals[6]
		transcript_biotype = ann_vals[7]
		hgvs_c = ann_vals[9]
		err_warn_msg = ann_vals[15]
	except:
		return


	# Check for Error messages
	if 'ERROR' in err_warn_msg:
		print(err_warn_msg, ' - Skipping entry...')



	# remove version from transcript id
	transcript_id = re.sub(r'\..*', '', transcript_id)

	# skip entries from structural variants
	if ':' in transcript_id:
		#print(transcript_id)
		return
	
	# Extract position of variant in the full coding reference sequence of this particular transcript
	coding_pattern = re.compile("c\.[0-9]")
	if not re.match(coding_pattern, hgvs_c):
		# skip gene upstream/downstream variants
		#print(hgvs_c)
		return
	elif re.search(r"\-|\+|\*", hgvs_c):
		# skip intron upstream/downstream variants
		#print(hgvs_c)
		return 
	


	coding_position = re.sub(r'c\.', '', hgvs_c)
	coding_position = re.sub(r'[AGCT].*', '', coding_position)
	# convert to amino-acid positions
	protein_position = str(int( (float(coding_position)-1)/3 ) + 1)



	if debug and (transcript_id) == debug_transcr:
		print('\n', chrom, position, ref, alt)
		print(annotations)
		print(transcript_id)
		print(hgvs_c)


	new_record = '=>'.join([chrom, position, ref, alt, transcript_id, protein_position])

	global_coord_map.add(new_record)



if __name__ == '__main__':
		
	global_coord_map = set()
	debug = False
	debug_transcr = 'ENST00000367028'  #'ENST00000335137'


	snpeff_results_file = '../../out/snpeff_results/all_exons_mutations.snpeff_annotation.vcf'
	print("Input file:", snpeff_results_file)


	genomic_coord_map_out_dir = "../../out/Genomic_to_protein_coords"
	if not os.path.exists(genomic_coord_map_out_dir):
		os.makedirs(genomic_coord_map_out_dir)



	cnt = 1
	with open(snpeff_results_file) as fh:
		for line in fh:
			if line.startswith('#'):
				continue

			line = line.rstrip()

			vals = line.split('\t')
			chrom, position, ref, alt = vals[0], vals[1], vals[3], vals[4]


			ann_field = vals[7].split(';')[0] # there may be some extra annotations ('LOF' or 'NMD') ";"-delimited, that we can discard
			ann_field_entries = ann_field.split(',')

			for ann in ann_field_entries:
				process_single_annotation_entry(ann, chrom, position, ref, alt)



			if cnt % 100000 == 0:
				print(cnt)
			if debug:
				if cnt > 200000:
					break
			cnt += 1


	if debug:
		print(global_coord_map)


	# Write global_coord_map into a custom-delimiter (=>) separated file, to then load with pandas and process it (split by transcript id, sort by p. coordinates)
	out_file = genomic_coord_map_out_dir + "/genome_coords_to_hgvs_p.tsv"
	if debug:
		out_file = genomic_coord_map_out_dir + "/" + debug_transcr + "_snpeff.dbg.res"

	print('> Writing global_coord_map to file...')
	out_fh = open(out_file, 'w')
	for record in global_coord_map:
		vals = record.split('=>')
		out_fh.write('\t'.join(vals) + '\n')

	print('Output file:', out_file)
