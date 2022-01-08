import sys

clinvar_dir = "../../data/clinvar/vcf_GRCh38"
out_fh = open(clinvar_dir + '/clinvar.pathogenic_subset.vcf', 'w')


with open(clinvar_dir + '/clinvar.vcf') as fh:
	for line in fh:
		if line.startswith('#'):
			out_fh.write(line)
			continue
		line = line.rstrip()

		vals = line.split('\t')
		#print(vals)
		
		info = vals[7]
		info_dict = dict(item.split("=") for item in info.split(";"))
		if 'CLNSIG' not in info_dict:
			continue
		else:
			#print(info_dict['CLNSIG'])
			simple_annot = ''
			if 'Pathogenic' in info_dict['CLNSIG']:
				simple_annot = 'Pathogenic'
			elif 'Likely_pathogenic' in info_dict['CLNSIG']:
				simple_annot = 'Likely_pathogenic'
			elif 'Uncertain_significance' in info_dict['CLNSIG']:
				simple_annot = 'Uncertain_significance'		
			else:
				continue

			out_fh.write('\t'.join(vals[:7] + [simple_annot]) + '\n')
	
out_fh.close()
