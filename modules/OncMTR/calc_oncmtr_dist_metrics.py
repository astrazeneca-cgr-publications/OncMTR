import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.image as mimg
import pickle
import numpy as np
import pandas as pd
import random
import sys, os
import traceback
import glob
import ntpath
import re
from collections import Counter
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from scipy.integrate import simps
from scipy.stats import wasserstein_distance, fisher_exact
from sklearn.metrics import pairwise


sys.path.insert(0,'../')
from utils.resource_mapping import ResourceMapping
from utils.distance_utils import cross_entropy, get_cdist_metrics

sys.path.insert(0, '../lolliplots')
import lolliplot

pd.set_option('display.max_rows', 100)



def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 'True', '1'):
		return True
	elif v.lower() in ('no', 'false', 'False', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')
		


class OncMtrAnnotator:
	
	def __init__(self, dataset, win_size, geneset, out_file_suffix='', verbose=False):
		"""
			Init variables, input/output directories and read input mtr score dataframes
		"""
		self.dataset = dataset
		self.win_size = win_size
		self.geneset = geneset
		self.out_file_suffix = out_file_suffix
		self.verbose = verbose
		
		# Output directory
		self.base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
			.loc['base_out'].absolute_path
		)


		# OncMTR output dir
		self.oncmtr_dir = self.base_out_dir + '/OncMTR-' + self.dataset + '-win' + str(self.win_size)
		if not os.path.exists(self.oncmtr_dir):
			os.makedirs(self.oncmtr_dir)

		# Read OncMTR signal pairs data frame
		oncmtr_pairs_file = self.oncmtr_dir + '/OncMTR_signal_pairs_df-' + self.dataset + '-win' + str(self.win_size) + '.tsv'
		self.oncmtr_pairs_df = pd.read_csv(oncmtr_pairs_file, sep='\t')
		if self.verbose:
			print(self.oncmtr_pairs_df.head())





	def pad_with_zeros(self, mtr_scores, mtr_ab_scores):

		if len(mtr_scores) == len(mtr_ab_scores):
			return mtr_scores, mtr_ab_scores

		max_len = max(len(mtr_scores), len(mtr_ab_scores))
		mtr_scores = mtr_scores + ([0] * abs(max_len - len(mtr_scores)))
		mtr_ab_scores = mtr_ab_scores + ([0] * abs(max_len - len(mtr_ab_scores)))

		return mtr_scores, mtr_ab_scores

	
	
	def neutralise_na_positions(self, mtr_scores, mtr_ab_scores):

		if np.isnan(mtr_scores).any() or np.isnan(mtr_ab_scores).any():

			#print('Imputing NAs...')

			mtr_na_indexes = list(np.argwhere(np.isnan(mtr_scores))[:, 0])
			mtr_ab_na_indexes = list(np.argwhere(np.isnan(mtr_ab_scores))[:, 0])	
			#print('mtr_na_indexes:', mtr_na_indexes)
			#print('mtr_ab_na_indexes:', mtr_ab_na_indexes)

			all_na_indexes = list(set(mtr_na_indexes) | set(mtr_ab_na_indexes))

			for i in all_na_indexes:
				mtr_scores[i] = 1
				mtr_ab_scores[i] = 1	

			# DBG	
			#mtr_nas = [mtr_scores[v] for v in mtr_na_indexes]
			#mtr_ab_nas = [mtr_ab_scores[v] for v in mtr_ab_na_indexes]
			#print(mtr_nas)
			#print(mtr_ab_nas)

		return mtr_scores, mtr_ab_scores




	def normalise_signals_by_diff_sign(self, mtr_scores, mtr_ab_scores, drop_neg_diffs):
		
		tmp_df = pd.DataFrame({'mtr': mtr_scores, 'mtr_ab': mtr_ab_scores})
		
		# drop regions with positive selection, i.e. mtr - mtr_ab < 0
		if drop_neg_diffs:
			tmp_df['diff'] = tmp_df.mtr - tmp_df.mtr_ab

			diff_indexes = tmp_df.loc[ tmp_df['diff'] < 0, : ].index.values.tolist()
			tmp_df.loc[diff_indexes, 'mtr_ab'] = tmp_df.loc[diff_indexes, 'mtr']

		return tmp_df['mtr'], tmp_df['mtr_ab']







	def calc_distance_metrics(self, mtr_scores, mtr_ab_scores, ndigits=4, simple_set=False):

		dist_dict = {}		


		cdist_dict = get_cdist_metrics(mtr_scores, mtr_ab_scores)
		
		cur_mtr_range = np.max(mtr_scores) - np.min(mtr_scores)
		dist_dict['transcr_len'] = len(mtr_scores)

		# ----- EMD [ equals 0 for identical signals ] -----
		dist_dict['emd'] = wasserstein_distance(mtr_scores, mtr_ab_scores)

		# --- AUC ---
		dist_dict['trapz_area_diff'] = np.trapz(mtr_scores) - np.trapz(mtr_ab_scores)
		dist_dict['simps_area_diff'] = simps(mtr_scores) - simps(mtr_ab_scores)
		dist_dict['trapz_area_ratio'] = np.trapz(mtr_scores) / np.trapz(mtr_ab_scores)
		dist_dict['simps_area_ratio'] = simps(mtr_scores) / simps(mtr_ab_scores)

		# --- Metrics of point differences ---
		diff_mtr = np.array(mtr_scores) - np.array(mtr_ab_scores)
		dist_dict['std_diff'] = np.std(diff_mtr)
		dist_dict['cum_diff'] = sum(diff_mtr)
		dist_dict['norm_cum_diff'] = sum(diff_mtr) / cur_mtr_range
		dist_dict['mean_diff'] = np.mean(diff_mtr)

		# --- Covariance ---
		covar = np.cov(mtr_scores, mtr_ab_scores)
		dist_dict['covar'] = covar[0][1]

		# --- Cross-correlation ---
		dist_dict['cross_correlation'] = np.correlate(mtr_scores, mtr_ab_scores)[0]

		# --- Cross-entropy ---
		dist_dict['cross_entropy'] = cross_entropy(mtr_scores, mtr_ab_scores)
		

		# keep "ndigits" sifgnificant digits
		for k, _ in dist_dict.items():
			dist_dict[k] = round(dist_dict[k], ndigits)
		


		full_dist_dict = {**cdist_dict, **dist_dict}


		if simple_set:

			selected_metrics = ['cross_entropy', 'simps_area_diff', 'simps_area_ratio', 'cum_diff', 'norm_cum_diff', 'mean_diff', 'cross_correlation']

			full_dist_dict = {}
			for metric in selected_metrics:
				full_dist_dict[metric] = dist_dict[metric]


		return full_dist_dict




	
	def get_signals_distance(self, enst_id, ens_gene, hgnc_name, drop_neg_diffs=True):

		global cnt

		try: 
			mtr_scores = self.oncmtr_pairs_df.loc[ self.oncmtr_pairs_df.ENST_ID == enst_id, 'MTR'].values[0]
			mtr_scores = np.array( eval(mtr_scores.replace('nan', 'np.nan')) )
		
			mtr_ab_scores = self.oncmtr_pairs_df.loc[ self.oncmtr_pairs_df.ENST_ID == enst_id, 'MTR_ab'].values[0]
			mtr_ab_scores = np.array( eval(mtr_ab_scores.replace('nan', 'np.nan')) )

		except:
			print("[Warning] No OncMTR data available - Skipping ", hgnc_name, '('+enst_id+')')
			return None

			

		# Add padding with zeros for vectors of not equal length
		mtr_scores, mtr_ab_scores = self.pad_with_zeros(mtr_scores, mtr_ab_scores)


		
		# Neutralise NA positions with MTR=1 to avoid affecting the signal comparison
		mtr_scores, mtr_ab_scores = self.neutralise_na_positions(mtr_scores, mtr_ab_scores)


		# DBG
		#save_gnb1_raw_mtr_values(mtr_scores, mtr_ab_scores)




		# Remove positive Onc-MTRs (negative MTR/MTR-AB differences) for distance calculations and plotting
		mtr_scores, mtr_ab_scores = self.normalise_signals_by_diff_sign(mtr_scores, mtr_ab_scores, drop_neg_diffs)
		dist_dict = self.calc_distance_metrics(mtr_scores, mtr_ab_scores)

		
		
		if self.verbose:
			print('self.verbose:', self.verbose)
			print(enst_id)
			for k,v in sorted(dist_dict.items()):
				print(k + ': ', v)
			print('\n\n')




		dist_dict['ENSG'] = ens_gene
		dist_dict['HGNC-tmp'] = hgnc_name


		return dist_dict
	




	

# =================== CUSTOM FUNCTIONS ===================	

def save_gnb1_raw_mtr_values(mtr_scores, mtr_ab_scores):
	"""
		Get raw MTR score values
	"""

	mtr_scores_as_str = [str(i) for i in mtr_scores]
	mtr_ab_scores_as_str = [str(i) for i in mtr_ab_scores]

	fh1 = open('mtr_scores_GNB1.txt', 'w')
	fh1.write(', '.join(mtr_scores_as_str))

	fh2 = open('mtr_scores_GNB1.AB_MEDIAN_filtered.txt', 'w')
	fh2.write(', '.join(mtr_ab_scores_as_str))

	#print(len(mtr_scores))
	#print(len(mtr_ab_scores))




def read_input_args():
	
	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	# === Required ===
	parser.add_argument("-d", dest="dataset", help="\tInput dataset", choices=['gnomad', 'UKB', 'topmed'], required=True)
	parser.add_argument("-w", dest="win_size", default=31, help="\tWindow size to calculate MTR \t(default=31)\n\n", required=True)
	
	
	# === Optional ===
	parser.add_argument("-n", dest="drop_neg_diffs", type=str2bool, default=True, help="\tFilter out negative OncMTR scores (i.e. tolerant regions) when plotting \t(default=True)\n\n", required=False)
	
	parser.add_argument("-v", dest="verbose", type=str2bool, default=False, help="\tPrint verbose messages \t(default=False)\n\n", required=False)
	parser.add_argument("-s", dest="out_file_suffix", default='', help="\tCustom output file suffix \t(default='')\n\n", required=False)
	parser.add_argument("-t", dest="test_run", type=str2bool, default=False, help="\tTest run for a single transcript \t(default=False)\n\n", required=False)
	
	
	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	args = parser.parse_args()
	print(args, '\n')

	return args






		




if __name__ == '__main__':

	
	args = read_input_args()

	# Read required arguments
	dataset = args.dataset
	win_size = args.win_size
	
	# Read optional arguments
	drop_neg_diffs = args.drop_neg_diffs
	verbose = args.verbose
	out_file_suffix = args.out_file_suffix
	test_run = args.test_run


	
		
	# =============== External annotation / resources ==============
	# Get ENST to ENSG, HGNC and Uniprot Ids mapping
	res_map = ResourceMapping()
	
	
	
	

	
	# ---- Gene selection default parameters ----
	enst_ids = None 
	gene_list = None
	# --

	
	# Test-Run (bypass specified geneset)
	if test_run:
		enst_ids = ['ENST00000378609'] # ['ENST00000264033', 'ENST00000256196', 'ENST00000332029']
		gene_list = ['GNB1']
		geneset = 'test-run'




	if enst_ids is None and gene_list is None: 

		geneset = 'All-genes'

		# ---- Get all ENST IDs ----
		print("\nRunning differential MTR for all ENST IDs...")
		# Output directory
		base_out_dir = pd.read_csv('../../.conf', index_col = 0, sep = '\t').loc['base_out'].absolute_path

		# OncMTR output dir
		cur_oncmtr_dir = base_out_dir + '/OncMTR-' + dataset + '-win' + str(win_size)
		cur_oncmtr_df = pd.read_csv(cur_oncmtr_dir + '/All_OncMTR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t')
		enst_ids = list(cur_oncmtr_df.ENST_ID.values)
		print('All ENST IDs:', len(enst_ids))
	
	else:
		if enst_ids:
			print("\nRunning for ENST sub-list: ", len(enst_ids), " transcripts...")

		elif gene_list:
			# get ENST IDs from specified HGNC list
			print("\nRunning for gene sub-list: ", len(gene_list), " genes...")
			enst_ids = res_map.get_hgnc_to_enst_mapping(gene_list)




	# ********************* Main RUN *********************
	annotator = OncMtrAnnotator(dataset, win_size, geneset, verbose=verbose)

	
	full_ENST_dist_dict = {}

	cnt = 1
	for enst_id in enst_ids:
		print('>', enst_id)

		try:
			ens_gene = res_map.enst_to_ensg_dict[enst_id]
			hgnc_name = res_map.enst_to_hgnc_dict[enst_id]
		except:
			print('Could not find enst_id:', enst_id)
			continue


		cur_results = annotator.get_signals_distance(enst_id, ens_gene, hgnc_name, drop_neg_diffs=drop_neg_diffs)	
		if cur_results is None:
			continue 
		
		#print('\n\n=====>>>', cur_results, '\n\n')
		full_ENST_dist_dict[enst_id] = cur_results


		cnt += 1
		if cnt % 100 == 0:
			print('Counter:', cnt)




	# Save all distance metrics to then plot distributions
	enstid_dist_df = pd.DataFrame(full_ENST_dist_dict).T
	enstid_dist_df.reset_index(inplace=True)
	enstid_dist_df.columns.values[0] = 'ENST'
	enstid_dist_df.insert(0, 'HGNC', enstid_dist_df['HGNC-tmp'])
	enstid_dist_df.drop('HGNC-tmp', inplace=True, axis=1)

	dist_out_file = annotator.oncmtr_dir + '/enstid_dist_metrics.' + dataset + '.tsv'
	if out_file_suffix != '':
		 dist_out_file += '.' + out_file_suffix
	
	
	enstid_dist_df.to_csv(dist_out_file, sep='\t', index=False)
	print('\n- Distance metrics saved into ' + dist_out_file)
