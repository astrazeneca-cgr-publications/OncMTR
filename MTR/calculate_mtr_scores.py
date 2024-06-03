import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import sys
import numpy as np
import math
from mtr import MtrClass
import glob


def load_obj(name ):
	with open(name, 'rb') as f:        
		return pickle.load(f)


def concat_results_files(out_dir, dataset, win_size):

        concat_out_file = out_dir + '/All_MTR_scores-' +  dataset + '-win' + str(win_size) + '.tsv'
        print("\nConcatenating MTR score files into a single file:", concat_out_file)
        scores_dir = out_dir + '/scores'

        all_files = glob.glob(scores_dir + '/*.tsv')

        cnt = 0
        full_df = pd.DataFrame()
        for f in all_files:
                try:
                        tmp_df = pd.read_csv(f, sep='\t', header=None)
                except:
                        continue

                if full_df.shape[0] > 0:
                        full_df = pd.concat([full_df, tmp_df], axis=0)
                else:
                        full_df = tmp_df

                cnt += 1
                if cnt % 500 == 0:
                        print(cnt)


        full_df.columns = ['Gene', 'ENST_ID', 'MTR']
        print(full_df.shape)

        full_df.to_csv(concat_out_file, sep='\t', index=False)




if __name__ == '__main__':

	input_file_str = sys.argv[1] 	# e.g. UKB_exons_mutations.snpeff_annotation
	dataset = sys.argv[2]	# e.g. UKB
	win_size = int(sys.argv[3])	# e.g. 31


	# Output directory
	base_out_dir = '../out'
	out_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size)



	# Read ENST to HGNC mappings
	with open('../data/ensembl/ENST_to_HGNC_dict.pkl', 'rb') as fh:
		enst_to_hgnc = pickle.load(fh)

	snpeff_processed_dir = "../out/snpeff_processed_out"

	# Input dataset (e.g. UKB, gnomAD) -  synonymous
	input_file = snpeff_processed_dir + "/" + input_file_str + ".vcf.synonymous.exon_coords.pkl"
	print('Synonymous - loading file:', input_file)
	input_cohort_syn_dict = load_obj(input_file)
	print(dataset + ' - synonymous:', len(input_cohort_syn_dict))

	# Input dataset (e.g. UKB, gnomAD) - missense 
	input_file = snpeff_processed_dir + "/" + input_file_str + ".vcf.missense.exon_coords.pkl"
	print('Missense - loading file:', input_file)
	input_cohort_mis_dict = load_obj(input_file)
	print(dataset + ' - missense:', len(input_cohort_mis_dict))




	# all synonymous
	input_file = snpeff_processed_dir + "/all_exons_mutations.snpeff_annotation.vcf.synonymous.exon_coords.pkl"
	all_syn_dict = load_obj(input_file)
	print('All synonymous:', len(all_syn_dict))

	# all missense
	input_file = snpeff_processed_dir + "/all_exons_mutations.snpeff_annotation.vcf.missense.exon_coords.pkl"
	all_mis_dict = load_obj(input_file)
	print('All missense:', len(all_mis_dict))



	all_unique_enst_ids = list(set( set(input_cohort_syn_dict.keys()) | set(input_cohort_mis_dict.keys()) | set(all_syn_dict.keys()) | set(all_mis_dict.keys()) ))
	print('All ENST IDs:', len(all_unique_enst_ids))



	cnt = 1
	mtr_stats = []
	for enst_id in all_unique_enst_ids:
	#for enst_id in ['ENST00000342066']:
		gene_id = enst_to_hgnc[enst_id]
		#print(cnt, enst_id, gene_id)

		cnt += 1
		if cnt % 500 == 0:
			print('Cnt:', cnt)

		obs_syn = input_cohort_syn_dict.get(enst_id, [])
		obs_mis = input_cohort_mis_dict.get(enst_id, [])
		exp_syn = all_syn_dict.get(enst_id, [])
		exp_mis = all_mis_dict.get(enst_id, [])

		mtr_obj = MtrClass(win_size=win_size, min_variants_thres=3, out_dir=out_dir)
		mtr_stats_i = mtr_obj.run(exp_syn, exp_mis, obs_syn, obs_mis, enst_id, gene_id=gene_id)			

		# ommit results with return value: -1
		if not isinstance(mtr_stats_i, int):
			mtr_stats.append(mtr_stats_i)


	mtr_stats = pd.concat(mtr_stats)

	mtr_stats.to_csv(out_dir + '/' + dataset + '-MTR_stats.csv', index=False, sep='\t')


	# Concatenate MTR score files per transcript into a single file
	concat_results_files(out_dir, dataset, win_size)
