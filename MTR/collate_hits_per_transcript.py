import sys
import pickle


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
	
	snpeff_processed_dir = "../out/snpeff_processed_out"
	
	cohort = sys.argv[1] # NFE_440k
	
	syn_input = snpeff_processed_dir + "/" + cohort + "_exons_mutations.snpeff_annotation.pass.vcf.synonymous.exon_coords"
	syn_dict = store_hits_in_dict(syn_input)
	mis_input = snpeff_processed_dir + "/" + cohort + "_exons_mutations.snpeff_annotation.pass.vcf.missense.exon_coords"
	mis_dict = store_hits_in_dict(mis_input)
