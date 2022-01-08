import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle
import glob
import pandas as pd
import sys
import numpy as np
import math
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from mtr import MtrClass
from statsmodels.stats.multitest import multipletests



def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', '1'):
        return True
    elif v.lower() in ('no', 'false', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def load_obj(dataset, filtering_str, var_type='synonymous'):

	snpeff_processed_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['snpeff_proc'].absolute_path
	)
	
	input_file = snpeff_processed_dir + "/" + dataset + "_exons_mutations.snpeff_annotation"  + filtering_str + ".vcf." + var_type + ".exon_coords.pkl"
	print('\n', var_type + ' - loading file:', input_file)


	with open(input_file, 'rb') as f:        
		input_cohort_dict = pickle.load(f)
		print(dataset + ' - ' + var_type + ':', len(input_cohort_dict))

	return input_cohort_dict


	
def concat_results_files(out_dir, dataset, win_size):

	concat_out_file = out_dir + '/All_MTR_scores-' +  dataset + '-win' + str(win_size) + '.tsv'
	print("\nConcatenating MTR score files into a single file:", concat_out_file)
	scores_dir = out_dir + '/scores'
	

	all_files = glob.glob(scores_dir + '/*.tsv')
	#print(all_files)

	cnt = 0
	full_df = pd.DataFrame()
	for f in all_files:
		try:
			tmp_df = pd.read_csv(f, sep='\t', header=None)
		except:
			continue

		if full_df.shape[0] > 0:
			full_df = pd.concat([full_df, tmp_df], axis=0)
			#print(full_df.shape)
		else:
			full_df = tmp_df
			#print(full_df.shape)

		cnt += 1
		if cnt % 500 == 0:
			print(cnt)


	full_df.columns = ['Gene', 'ENST_ID', 'MTR']
	#print(full_df.head())
	#print(full_df.info())
	print(full_df.shape)

	full_df.to_csv(concat_out_file, sep='\t', index=False)

	



if __name__ == '__main__':

	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-d", dest="dataset", help="Input dataset", choices=['gnomad', 'UKB', 'topmed', 'UKB_gnomad'], required=True)
	parser.add_argument("-w", dest="win_size", default=31, help="Window size to calculate MTR", required=True)
	parser.add_argument("-s", dest="out_suffix", default='', help="Custom output folder suffix", required=False)
	parser.add_argument("-f", dest="ab_filtered", type=str2bool, default=False, help="Filter out based on AB_Median {yes, true, 1 | no, false, 0}", required=False)

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)


	args = parser.parse_args()
	print(args, '\n')

	dataset = args.dataset
	win_size = args.win_size
	out_suffix = args.out_suffix
	ab_filtered = args.ab_filtered



	# Output directory
	base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['base_out'].absolute_path
	)

	out_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size)
	if ab_filtered:
		out_dir += '-AB-filtered'
	if out_suffix not in ['', 'None']:
		out_dir += '-' + out_suffix
	print('Output dir:', out_dir)


	filtering_str = ''
	if 'gnomad' in dataset and 'UKB' not in dataset:
		# PASS variants are available for the gnomad dataset only and are used by default instead of all variants
		filtering_str = '.pass'
	if ab_filtered:
		filtering_str += '.ab_filtered'


	# Read ENST to HGNC mappings
	enst_to_hgnc_fname = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['enst_hgnc'].absolute_path
	)
	with open(enst_to_hgnc_fname, 'rb') as fh:
		enst_to_hgnc = pickle.load(fh)





	# Input dataset (e.g. UKB, gnomAD) -  synonymous
	input_cohort_syn_dict = load_obj(dataset, filtering_str, var_type='synonymous')
	# Input dataset (e.g. UKB, gnomAD) - missense 
	input_cohort_mis_dict = load_obj(dataset, filtering_str, var_type='missense')




	# all synonymous
	all_syn_dict = load_obj('all', '', var_type='synonymous')
	# all missense
	all_mis_dict = load_obj('all', '', var_type='missense')



	all_unique_enst_ids = list(set( set(input_cohort_syn_dict.keys()) | set(input_cohort_mis_dict.keys()) | set(all_syn_dict.keys()) | set(all_mis_dict.keys()) ))
	print('All ENST IDs:', len(all_unique_enst_ids))




	# Calculate MTR per transcript
	cnt = 1
	mtr_stats = []
	fdr_stats = []
	for enst_id in all_unique_enst_ids:
	#for enst_id in ['ENST00000378609']:
		
		try:
			gene_id = enst_to_hgnc[enst_id]
		except:
			print('[Warning] No matching HGNC found for' + enst_id + ' - continuning with ENST_ID info only.')
			gene_id = enst_id
		#print(cnt, enst_id, gene_id)

		cnt += 1
# 		if cnt == 10:
# 			break
		if cnt % 1000 == 0:
			print(cnt)

		obs_syn = input_cohort_syn_dict.get(enst_id, [])
		obs_mis = input_cohort_mis_dict.get(enst_id, [])
		exp_syn = all_syn_dict.get(enst_id, [])
		exp_mis = all_mis_dict.get(enst_id, [])
	
		#print("enst_id\t", enst_id, "; obs_syn\t", obs_syn, "; obs_mis\t", obs_mis, "; exp_syn\t", exp_syn, "'; exp_mis\t", exp_mis)
	
	
		mtr_obj = MtrClass(win_size=win_size, min_variants_thres=3, out_dir=out_dir)
		mtr_stats_i, fdr_stats_i  = mtr_obj.run(exp_syn, exp_mis, obs_syn, obs_mis, enst_id, gene_id=gene_id, make_plot=True)

		# ommit results with return value: -1
		if not isinstance(mtr_stats_i, int):
			mtr_stats.append(mtr_stats_i)
			fdr_stats.append(fdr_stats_i)


	mtr_stats = pd.concat(mtr_stats)
	fdr_stats = pd.concat(fdr_stats)
	
	# save MTR-stats to file
	mtr_stats.to_csv(out_dir + '/' + dataset + '-win' + str(win_size) + '.MTR_stats.csv', index=False, sep='\t')
	
	# save FDR-stats to file
	fdr_stats = fdr_stats.fillna(1.0).assign(fdr = lambda x: np.round(multipletests(x.fdr.values.astype(float), method = 'fdr_bh')[1], 5))
	fdr_stats.to_csv(out_dir + '/' + dataset + '-win' + str(win_size) + '.fdr_stats.csv', index=False, sep='\t')
	
	
	fdr_stats_all = ( 
		fdr_stats[['fdr', 'enst']]
		.groupby(['enst'])
		.apply(lambda x: ', '.join(list(np.round(x.fdr.values, 5).astype(str))))
		.reset_index()
		.rename(columns = {0: 'FDR'})
		.merge(fdr_stats[['enst', 'gene_id']].drop_duplicates(), on = 'enst')
		.rename(columns = {'enst': 'ENST_ID', 'gene_id': 'Gene'})
		[['Gene', 'ENST_ID', 'FDR']]
	)
	
	fdr_stats_all.to_csv(out_dir + '/All_FDR_scores-' +  dataset + '-win' + str(win_size) + '.tsv', index=False, sep='\t')
	

	# Concatenate MTR score files per transcript into a single file
	concat_results_files(out_dir, dataset, win_size)
