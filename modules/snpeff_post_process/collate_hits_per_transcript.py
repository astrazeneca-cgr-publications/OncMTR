import sys
import pickle
import pandas as pd
import numpy as np


def save_obj(obj, name ):
	with open(name + '.pkl', 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def store_hits_in_dict(input_file):

	print('\n> ' + input_file)
	
	target_dict = {}
	cnt = 1
	with open(input_file) as fh:
		for line in fh:
			line = line.rstrip()
			enst_id, position = line.split('\t')
			if '+' in position: 
				continue
# 				position = position.split('+')[0]
			if '-' in position: 
				continue
# 				position = position.split('-')[0]
			if '*' in position: 
				continue
# 				position = position.replace('*', '')
			if position == '': 
				continue
			position = int(position)

			if enst_id in target_dict:
				target_dict[enst_id].append(position)
			else:
				target_dict[enst_id] = [position]

			if cnt % 1000000 == 0:
				print(cnt)
			cnt += 1
	print('Total transcripts:', len(target_dict))

	# save dict to pickle file
	save_obj(target_dict, input_file)

	return target_dict
	

if __name__ == '__main__':
	
	snpeff_processed_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['snpeff_proc'].absolute_path
	)
	



# >>> gnomAD
print('> gnomAD')
gnomad_syn_input = snpeff_processed_dir + "/gnomad_exons_mutations.snpeff_annotation.pass.vcf.synonymous.exon_coords"
gnomad_syn_dict = store_hits_in_dict(gnomad_syn_input)
gnomad_mis_input = snpeff_processed_dir + "/gnomad_exons_mutations.snpeff_annotation.pass.vcf.missense.exon_coords"
gnomad_mis_dict = store_hits_in_dict(gnomad_mis_input)

# gnomAD: AB_MEDIAN > 0.30
print('> gnomAD (AB_MEDIAN > 0.30)')
gnomad_ab_syn_input = snpeff_processed_dir + "/gnomad_exons_mutations.snpeff_annotation.pass.ab_filtered.vcf.synonymous.exon_coords"
gnomad_ab_syn_dict = store_hits_in_dict(gnomad_ab_syn_input)
gnomad_ab_mis_input = snpeff_processed_dir + "/gnomad_exons_mutations.snpeff_annotation.pass.ab_filtered.vcf.missense.exon_coords"
gnomad_ab_mis_dict = store_hits_in_dict(gnomad_ab_mis_input)




# > all
print('> all')
all_syn_input = snpeff_processed_dir + "/all_exons_mutations.snpeff_annotation.vcf.synonymous.exon_coords"
all_syn_dict = store_hits_in_dict(all_syn_input)
all_mis_input = snpeff_processed_dir + "/all_exons_mutations.snpeff_annotation.vcf.missense.exon_coords"
all_mis_dict = store_hits_in_dict(all_mis_input)
