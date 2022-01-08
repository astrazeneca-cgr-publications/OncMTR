import re
import sys

out_fh = open('clinvar.pathogenic_subset.snpeff_annotation.coding.tsv', 'w')

cnt = 0
with open('clinvar.pathogenic_subset.snpeff_annotation.vcf') as fh:
	for line in fh:
		if line.startswith('#'):
			continue
		line = line.rstrip()
		vals = line.split('\t')

		info = vals[7]
		tmp_entries = info.split(',')

		clinvar_annot = tmp_entries[0].split(';')[0]
		#print('> Clinvar (simplified) annotation:', clinvar_annot, '\n')

		for entry in tmp_entries:
			entry_vals = entry.split('|')
			if len(entry_vals) < 16:
				continue

			#print(entry_vals)
			try:
				hgvs_p = entry_vals[10]
				#print(hgvs_p)

				# ignore entries without coding variants
				if hgvs_p == '':
					continue
				else:
					hgvs_p = hgvs_p.replace('p.', '')

				gene_name = entry_vals[3]
				ens_gene = entry_vals[4]
				enst_id = re.sub('\..*', '', entry_vals[6]) 
				
				#if gene_name == 'TP53':
				#	if hgvs_p != '':
				#		print(entry)
				#		sys.exit()
			except:
				continue

			new_entry = gene_name + '\t' + ens_gene + '\t' + enst_id + '\t' + hgvs_p + '\t' + clinvar_annot + '\n'
			out_fh.write(new_entry)

			cnt += 1
			if cnt % 10000 == 0:
				print(cnt)

	
out_fh.close()

