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

from gene_datasets import get_gene_dataset

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
	
	def __init__(self, dataset, win_size, geneset, out_suffix='', verbose=False):
		"""
			Init variables, input/output directories and read input mtr score dataframes
		"""
		self.dataset = dataset
		self.win_size = win_size
		self.geneset = geneset
		self.out_suffix = out_suffix
		self.verbose = verbose
		
		# Output directory
		self.base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
			.loc['base_out'].absolute_path
		)


		# OncMTR output dir
		self.oncmtr_dir = self.base_out_dir + '/OncMTR-' + self.dataset + '-win' + str(self.win_size)
		if self.out_suffix not in ['', 'None']:
			self.oncmtr_dir += '-' + self.out_suffix
		if not os.path.exists(self.oncmtr_dir):
			os.makedirs(self.oncmtr_dir)

		# Read OncMTR signal pairs data frame
		oncmtr_pairs_file = self.oncmtr_dir + '/OncMTR_signal_pairs_df-' + self.dataset + '-win' + str(self.win_size) + '.tsv'
		self.oncmtr_pairs_df = pd.read_csv(oncmtr_pairs_file, sep='\t')
		if self.verbose:
			print(self.oncmtr_pairs_df.head())



			
		# Annotation results base dir
		self.res_base_dir = self.oncmtr_dir + '/Plots_and_Variant-annotations' 
		if not os.path.exists(self.res_base_dir):
			os.makedirs(self.res_base_dir)

		self.dataset_res_dir = self.res_base_dir + '/' + self.geneset
		if not os.path.exists(self.dataset_res_dir):
			os.makedirs(self.dataset_res_dir)

		# dataset-specific annotation out dir
		self.res_dir = self.dataset_res_dir + '/' + variant_annot_dataset + '_variants'
		if keep_significant_genes_only:
			self.res_dir += '-signif_hits'	
		else:
			self.res_dir += '-all_hits'	
		self.res_dir += '-top_perc' + str(top_percentile)
		if not os.path.exists(self.res_dir):
			os.makedirs(self.res_dir)
		
			
		self.signif_combined_figs_dir = self.res_dir + '/Significant_Diff-MTR_and_lolliplots'
		if not os.path.exists(self.signif_combined_figs_dir):
			os.makedirs(self.signif_combined_figs_dir)

		self.diff_mtr_figs_dir = self.res_dir + '/Diff-MTR_plots'
		if not os.path.exists(self.diff_mtr_figs_dir):
			os.makedirs(self.diff_mtr_figs_dir)
			
		self.oncmtr_figs_dir = self.res_dir + '/OncMTR_plots'
		if not os.path.exists(self.oncmtr_figs_dir):
			os.makedirs(self.oncmtr_figs_dir)
			


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
		
		if drop_neg_diffs:
			tmp_df['diff'] = tmp_df.mtr - tmp_df.mtr_ab

			diff_indexes = tmp_df.loc[ tmp_df['diff'] < 0, : ].index.values.tolist()
			tmp_df.loc[diff_indexes, 'mtr_ab'] = tmp_df.loc[diff_indexes, 'mtr']

		return tmp_df['mtr'], tmp_df['mtr_ab']




	def process_diff_mtrs_for_stat_test(self, enst, hgnc_name, mtr_scores, mtr_ab_scores, variant_annot_dataset, top_percentile):

		diff_mtr = mtr_ab_scores - mtr_scores
		# BETA - fix probable offest of codon indexes
		# [Solution]: No adjustment is required based on empirical results
		#diff_mtr = diff_mtr.iloc[1: ].reset_index(drop=True)  # index - 1
		#diff_mtr.index += 1  # index + 1


		# TODO: define better the 'negative'/'positive' sets [DONE: getting bottom 10-20% percentile of neg vals and all zeros for the rest]
		# 'negative' set: set with negative diff-MTRs
		# 'positive' set: set with 0 diff-MTRs (set drop_neg_diffs to False, to make sure > 0 diff-MTRs are not coerced into 0)
		neg_vals = diff_mtr[diff_mtr < 0]
		if neg_vals.shape[0] > 0:
			neg_vals = neg_vals.tolist()
		else:
			return -1

		diff_upper_thres = np.percentile(neg_vals, top_percentile)
		if diff_upper_thres == 0:
			return -1
		std_neg_vals = np.std(neg_vals)
		neg_range = max(neg_vals) - min(neg_vals)
		z_score_neg_range = neg_range / std_neg_vals
		avg_neg_val = np.mean(neg_vals)

		undiff_lower_thres = 0

		
		"""
		diff_upper_thres = np.percentile(diff_mtr, top_percentile)
		undiff_lower_thres = np.percentile(diff_mtr, 100 - top_percentile)
		median_diff = np.median(diff_mtr)	

		# Make sure I'm selecting codons that have mtr_ab - mtr != 0
		while top_percentile > 0:
			if diff_upper_thres == 0:
			# define lower percentiles if not
				top_percentile -= 1
				diff_upper_thres = np.percentile(diff_mtr, top_percentile)
			else:
				break
		else:
			return -1

		#print('Diff-MTR', top_percentile, '- top_percentile:', diff_upper_thres)
		#print('Non-Diff-MTR', str(100 - top_percentile), '- top_percentile:', undiff_lower_thres)
		#print('median:', median_diff)
		"""


		diff_codon_indexes = diff_mtr.loc[ diff_mtr <= diff_upper_thres].index.tolist()
		# BETA
		#diff_codon_indexes = diff_mtr.loc[ diff_mtr <= (avg_neg_val) ].index.tolist()
		
		#undiff_codon_indexes = diff_mtr.loc[ diff_mtr >= undiff_lower_thres].index.tolist()  # use when drop_neg_diffs = True
		undiff_codon_indexes = diff_mtr.loc[ diff_mtr == 0].index.tolist()	# Probably use that as default in all cases; [legacy: use when drop_neg_diffs = False]

		"""
		print('**', enst, '**')
		print('> diff_upper_thres:', diff_upper_thres)
		print('> std_neg_vals:', std_neg_vals)
		print('> neg_range:', neg_range)
		print('> z_score_neg_range:', z_score_neg_range)
		print('> avg_neg_val:', avg_neg_val)

		print('\nSelected neg vals:', diff_mtr.loc[diff_codon_indexes], '\n')
		print('\nAll neg vals:', diff_mtr.loc[ diff_mtr < 0], '\n')
		print('\ndiff_codon_indexes:', len(diff_codon_indexes))
		print(diff_codon_indexes)
		print('\nundiff_codon_indexes:', len(undiff_codon_indexes))
		print(undiff_codon_indexes)
		"""


		if len(diff_codon_indexes) > len(undiff_codon_indexes):
			sampled_undiff_codon_indexes = undiff_codon_indexes[:]
		else:
			sampled_undiff_codon_indexes = random.sample(undiff_codon_indexes, len(diff_codon_indexes))


		
		if variant_annot_dataset == 'oncology':
			df = pd.read_csv('../../data/cancer-hotspots/transcripts_hotspots_oncology.txt', sep='\t')
			df = df.loc[ df['Transcript_ID'] == enst, :]
			
			if df.empty:
				return -1
			
			df.reset_index(inplace=True)
			df.dropna(inplace=True)
			#print(df.head())
			
			hgvs_p_list = df['HGVSp'].tolist()
			
			hgvs_p_list = sorted([int(re.sub("[^0-9]", "", p)) for p in hgvs_p_list])
			#print('hgvs_p_list:', hgvs_p_list)
			#print('length:', len(hgvs_p_list))
			
		elif variant_annot_dataset == 'clinvar':
			df = clinvar_full_annot_df.loc[ clinvar_full_annot_df.ENST == enst_id, :].copy()
			if df.empty:
				return -1
		
			df.reset_index(inplace=True)

			hgvs_p_list = sorted([ int(re.sub("[^0-9]", "", p)) for p in df['HGVS_P'].tolist()] )

			#print('hgvs_p_list:', hgvs_p_list)
			#print('length:', len(hgvs_p_list))
		
		
		
		#print('diff_codon_indexes:', diff_codon_indexes)


		# - Get variant indexes with duplicates
		diff_codons_with_variants = [v for v in hgvs_p_list if v in diff_codon_indexes]
		undiff_codons_with_variants = [v for v in hgvs_p_list if v in undiff_codon_indexes]
		sampled_undiff_codons_with_variants = [v for v in hgvs_p_list if v in sampled_undiff_codon_indexes]


		# - Get unique variant indexes
		uniq_diff_codons_with_variants = list(set(hgvs_p_list) & set(diff_codon_indexes))
		uniq_undiff_codons_with_variants = list(set(hgvs_p_list) & set(undiff_codon_indexes))
		uniq_sampled_undiff_codons_with_variants = list(set(hgvs_p_list) & set(sampled_undiff_codon_indexes))

		# Compile Fisher 2x2 tables
		diff_mtr_counts = [len(diff_codons_with_variants), len(diff_codon_indexes) - len(uniq_diff_codons_with_variants)]
		undiff_mtr_counts = [len(undiff_codons_with_variants), len(undiff_codon_indexes) - len(uniq_undiff_codons_with_variants)]
		sampled_undiff_mtr_counts = [len(sampled_undiff_codons_with_variants), len(sampled_undiff_codon_indexes) - len(uniq_sampled_undiff_codons_with_variants)]
		
		
		full_fisher_table = [diff_mtr_counts, undiff_mtr_counts]
		sampled_fisher_table = [diff_mtr_counts, sampled_undiff_mtr_counts]


		full_fisher_res = fisher_exact(full_fisher_table)
		sampled_fisher_res = fisher_exact(sampled_fisher_table)

		full_odds_ratio, full_p_val = full_fisher_res
		sampled_odds_ratio, sampled_p_val = sampled_fisher_res


		# ommit plotting non-significant results (if 'keep_significant_genes_only' is set to True)
		if keep_significant_genes_only:
			if full_odds_ratio <= 1 or full_p_val >= 0.05:
				return -1
			else:
				signif_genes[hgnc_name] = 1
				signif_transcripts[enst] = 1
			
		print('\nHGNC:', hgnc_name, ' (' + enst + ')')
		print('- diff_mtr_counts:', diff_mtr_counts)
		print('- undiff_mtr_counts:', undiff_mtr_counts)
		print('- sampled_undiff_mtr_counts:', sampled_undiff_mtr_counts)

		print('\n> Full Fisher table:', full_fisher_res)
		print('> Sampled Fisher table:', sampled_fisher_res)
		print('----------------------------\n\n')


		fisher_result_str = '\nFisher\'s exact test - OR: ' + str(round(full_odds_ratio,4)) + '; p-val: ' + str(round(full_p_val, 5)) + '\n'
		#fisher_result_str += '(Sampled) Fisher test - OR: ' + str(sampled_odds_ratio) + '; p-val: ' + str(sampled_p_val) + '\n'
		
		fet_res_dict = {'FET Odds ratio (Full)': full_odds_ratio, 'FET P-value (Full)': full_p_val, 'FET Odds ratio (Sampled)': sampled_odds_ratio, 'FET P-value (Sampled)': sampled_p_val}

		
		return fet_res_dict, fisher_result_str




	def calc_distance_metrics(self, mtr_scores, mtr_ab_scores, ndigits=4, simple_set=True):

		dist_dict = {}		


		cdist_dict = get_cdist_metrics(mtr_scores, mtr_ab_scores)
		

		# ----- EMD [ equals 0 for identical signals ] -----
		dist_dict['emd'] = wasserstein_distance(mtr_scores, mtr_ab_scores)

		# --- AUC ---
		dist_dict['trapz_area_diff'] = np.trapz(mtr_scores) - np.trapz(mtr_ab_scores)
		dist_dict['simps_area_diff'] = simps(mtr_scores) - simps(mtr_ab_scores)
		#dist_dict['trapz_area_ratio'] = np.trapz(mtr_scores) / np.trapz(mtr_ab_scores)
		#dist_dict['simps_area_ratio'] = simps(mtr_scores) / simps(mtr_ab_scores)

		# --- Metrics of point differences ---
		diff_mtr = np.array(mtr_scores) - np.array(mtr_ab_scores)
		dist_dict['std_diff'] = np.std(diff_mtr)
		dist_dict['cum_diff'] = sum(diff_mtr)
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
			full_dist_dict = {'cross_entropy': dist_dict['cross_entropy'], 
					  'simps_area_diff': dist_dict['simps_area_diff']}

		return full_dist_dict




	
	def get_signals_distance(self, enst_id, ens_gene, hgnc_name, drop_neg_diffs=True, make_plots=True, top_percentile=None):

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


		mtr_scores_original, mtr_ab_scores_original = mtr_scores[:], mtr_ab_scores[:]


		mtr_scores, mtr_ab_scores = self.normalise_signals_by_diff_sign(mtr_scores, mtr_ab_scores, drop_neg_diffs)
		# BETA: exclude transcripts that have the exact same signals (zero euclidean distance)
		#if sum(mtr_scores - mtr_ab_scores) == 0:
		#	return
		
		is_significant = False
		fet_res_dict = {}
		fet_res_str = ''
		try:
			fet_res_dict, fet_res_str = self.process_diff_mtrs_for_stat_test(enst_id, hgnc_name, mtr_scores, mtr_ab_scores, variant_annot_dataset, top_percentile)
			
			if 'FET P-value (Full)' in fet_res_dict:
				if float(fet_res_dict['FET P-value (Full)']) < 0.05:
					is_significant = True
					print('\n-- Significant FET enrichment:', enst_id, hgnc_name, fet_res_dict['FET P-value (Full)'])

		except Exception as e:
			print('[Warning] Could not run the statistical test for variant enrichment in top Onc-MTR regions.\n')
			
			




		# Remove positive Onc-MTRs (negative MTR/MTR-AB differences) for distance calculations and plotting
		mtr_scores, mtr_ab_scores = self.normalise_signals_by_diff_sign(mtr_scores_original, mtr_ab_scores_original, drop_neg_diffs)
		dist_dict = self.calc_distance_metrics(mtr_scores, mtr_ab_scores)

		
		
		if self.verbose:
			print('self.verbose:', self.verbose)
			print(enst_id)
			for k,v in sorted(dist_dict.items()):
				print(k + ': ', v)
			print('\n\n')



		def plot_onc_mtr(fet_res_str):
			a = fig.add_subplot(2,1,1)
			oncmtr_scores = mtr_ab_scores - mtr_scores
			plt.plot(oncmtr_scores, linewidth=0.8, label='OncMTR', color='#3182bd')
			#plt.ylim(-0.1, max_mtr + 0.1)
			
			plt.legend(fontsize=8)
			plot_title = 'OncMTR: ' +  hgnc_name + ' - ' + enst_id + '\n'
			tmp_cnt = 0  
			for k,v in sorted(dist_dict.items()):
				plot_title += k + ': ' + str(v) + ';  '
				
				tmp_cnt += 1
				if tmp_cnt % 5 == 0:
					plot_title += '\n'
		
			plot_title += '\n' + str(fet_res_str)
			
			plt.title(plot_title, fontsize=7) 
			plt.grid(linewidth=0.5, linestyle='-')

			




		def plot_diff_mtr(fet_res_str):
			a = fig.add_subplot(2,1,1)
			plt.plot(mtr_ab_scores, linewidth=0.8, label='AB_MEDIAN > 0.3', color='#fd8d3c')
			plt.plot(mtr_scores, linewidth=0.8, label='All gnomAD PASS variants', color='#bdbdbd')
			#plt.ylim(-0.1, max_mtr + 0.1)
			
			plt.legend(fontsize=8)
			plot_title = hgnc_name + ' - ' + enst_id + '\n'
			tmp_cnt = 0  
			for k,v in sorted(dist_dict.items()):
				plot_title += k + ': ' + str(v) + ';  '
				
				tmp_cnt += 1
				if tmp_cnt % 5 == 0:
					plot_title += '\n'
		
			plot_title += '\n' + str(fet_res_str)
			
			plt.title(plot_title, fontsize=7) 
			plt.grid(linewidth=0.5, linestyle='-')

			

		signif_combined_fig_file = self.signif_combined_figs_dir + '/' + hgnc_name + '.' + enst_id + '.' + variant_annot_dataset + '-MTR_lolliplot.png'
		plot_has_been_created = False
		if os.path.exists(signif_combined_fig_file):
			plot_has_been_created = True


		if make_plots and not plot_has_been_created and is_significant:
			print("\n- Plotting two MTR signals (full and AB_MEDIAN-filtered) and lolliplot...")
			
			# Save OncMTR into oncmtr_figs_dir in a separate file
			fig = plt.figure(figsize=(12, 10))
			plot_onc_mtr(fet_res_str)
			fig.savefig(self.oncmtr_figs_dir + '/' + hgnc_name + '.' + enst_id + '.onc-MTR.pdf', bbox_inches='tight')
			plt.close()
			
			fig = plt.figure(figsize=(12, 10))
			plot_diff_mtr(fet_res_str)
			fig.savefig(self.diff_mtr_figs_dir + '/' + hgnc_name + '.' + enst_id + '.diff-MTR.pdf', bbox_inches='tight')
			plt.close()


			# Plot OncMTR signal and a lollipop plot
			fig = plt.figure(figsize=(12, 10))
			plot_onc_mtr(fet_res_str)

			# Lollipop plot
			try:
				if hgnc_name == 'BRAF':
					uniprot_id = 'P15056'
				#elif enst_id in enst_to_uniprot_dict:
				#	uniprot_id = enst_to_uniprot_dict[enst_id]
				else:
					uniprot_id = res_map.ensg_to_uniprot_dict[ens_gene]
					#print(ens_gene)
					#print(uniprot_id)
			except Exception as e:
				print('[Warning] No Uniprot id matching', enst_id, '\n\n')



			try:
				# TODO: 
				# catch failed runs that don't have a valid Uniprot ID when calling lollipops

				out_dir = self.res_dir + '/' + variant_annot_dataset + '-lolliplots'
				if not os.path.exists(out_dir):
					os.makedirs(out_dir)

				if variant_annot_dataset == 'oncology':
					out_img_file = lolliplot.onc_hotspots_lolliplot_for_transcript(enst_id, uniprot_id, out_dir, out_img_type, method=lolli_method, variant_annot_dataset=variant_annot_dataset, verbose=self.verbose)

				elif variant_annot_dataset == 'clinvar':
					df = clinvar_full_annot_df.loc[ clinvar_full_annot_df.ENST == enst_id, :].copy()
		

					out_img_file = lolliplot.clinvar_lolliplot_for_transcript(df, enst_id, uniprot_id, out_dir, out_img_type, method=lolli_method, variant_annot_dataset=variant_annot_dataset, verbose=self.verbose)
				print(out_img_file)


				lolli_img = mimg.imread(out_img_file)

				a = fig.add_subplot(2,1,2)
				plt.imshow(lolli_img, aspect=1)
				plt.axis('off')

			except Exception as e:
				print('[Warning] No lolliplot generated for', enst_id)
				print(e, '\n')
				

				
			fig.savefig(signif_combined_fig_file, dpi=500, bbox_inches='tight')
			plt.close()




		dist_dict['ENSG'] = ens_gene
		dist_dict['HGNC'] = hgnc_name

		fet_res_dict['ENSG'] = ens_gene
		fet_res_dict['HGNC'] = hgnc_name

		cur_results = {'distances': dist_dict, 'fet': fet_res_dict}

		return cur_results
	




	

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





def assess_enst_num_with_known_variants(enst_ids, variant_annot_dataset):
	"""
		Checking presence of known variant annotation for input ENST IDs/Genes
	"""

	# get total input HGNC genes
	input_hgnc_genes = {}
	for enst in enst_ids:
		try:
			input_hgnc_genes[ enst_to_hgnc_dict[enst] ] = 1
		except:
			continue


	if variant_annot_dataset == 'clinvar':
		sub_df = clinvar_full_annot_df.loc[ clinvar_full_annot_df.ENST.isin(enst_ids) ]

		uniq_present_hgnc_genes = sub_df['HGNC'].unique()
		uniq_present_enst_ids = sub_df['ENST'].unique()
	
	elif variant_annot_dataset == 'oncology':

		df = pd.read_csv('../../data/cancer-hotspots/transcripts_hotspots_oncology.txt', sep='\t')
		df.dropna(inplace=True)
		sub_df = df.loc[ df['Transcript_ID'].isin(enst_ids), :]
		
		uniq_present_hgnc_genes = sub_df['Hugo_Symbol'].unique()
		uniq_present_enst_ids = sub_df['Transcript_ID'].unique()
		

	print('- Present HGNC genes:', len(uniq_present_hgnc_genes))
	print('  Total input HGNC genes (non-accurate):', len(input_hgnc_genes.keys()))

	print('\n- Present enst_ids:', len(uniq_present_enst_ids))
	print('  Total input enst_ids:', len(enst_ids))






def read_clinvar_annot(drop_vars_of_uncert_signif=True):

	clinvar_full_annot_df = pd.read_csv('../../data/clinvar/vcf_GRCh38/filtered-clinvar-files/clinvar_pathogenic.mis_syn.non-somatic.tsv', sep='\t', low_memory=False)

	if drop_vars_of_uncert_signif:
		print('Removing variants of uncertain significance from ClinVar...')
		clinvar_full_annot_df = clinvar_full_annot_df.loc[ clinvar_full_annot_df.CLNSIG != 'Uncertain_significance' ]
		print(clinvar_full_annot_df.shape, '\n')		

	return clinvar_full_annot_df





	
def read_input_args():
	
	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	# === Required ===
	parser.add_argument("-d", dest="dataset", help="\tInput dataset", choices=['gnomad', 'UKB', 'topmed'], required=True)
	parser.add_argument("-w", dest="win_size", default=31, help="\tWindow size to calculate MTR \t(default=31)\n\n", required=True)
	
	
	# === Optional ===
	parser.add_argument("-a", dest="variant_annot_dataset", default='clinvar', choices=['clinvar', 'oncology'], help="\tVariant annotation dataset \t(default=clinvar)\n\n", required=False)
	parser.add_argument("-f", dest="keep_significant_genes_only", type=str2bool, default=False, help="\tKeep significant genes only [yes, true, True, 1 | no, false, False, 0] \t(default=False)\n\n", required=False)
	parser.add_argument("-t", dest="top_percentile", default=10, help="\tTop percentile of OncMTR regions within a gene to assess for variant enrichment \t(default=10)\n\n", required=False)
	parser.add_argument("-g", dest="geneset", default=None, help="\tGeneset to annotate with variants and examine enrichment with Fisher's exact test \t(default=None)\n\n", required=False)
	

	parser.add_argument("-m", dest="make_plots", type=str2bool, default=True, help="\tMake plots for Diff MTR & lolliplots \t(default=True)\n\n", required=False)
	parser.add_argument("-i", dest="out_img_type", default='png', help="\tPlot output image type \t(default=png)\n\n", required=False)
	parser.add_argument("-l", dest="lolli_method", default='lollipops', choices=['lollipops', 'genvisr'], help="\tMethod to use for plotting lolliplots \t(default=lollipops)\n\n", required=False)
	
	parser.add_argument("-n", dest="drop_neg_diffs", type=str2bool, default=True, help="\tFilter out negative OncMTR scores (i.e. tolerant regions) when plotting \t(default=True)\n\n", required=False)
	parser.add_argument("-u", dest="drop_vars_of_uncert_signif", type=str2bool, default=False, help="\tFilter out variants of uncertain significance from ClinVar when plotting a lolliplot \t(default=False)\n\n", required=False)
	
	parser.add_argument("-v", dest="verbose", type=str2bool, default=False, help="\tPrint verbose messages \t(default=False)\n\n", required=False)
	parser.add_argument("-s", dest="out_suffix", default='', help="\tCustom output folder suffix \t(default='')\n\n", required=False)
	
	
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
	variant_annot_dataset = args.variant_annot_dataset
	keep_significant_genes_only = args.keep_significant_genes_only
	top_percentile = args.top_percentile
	geneset = args.geneset

	make_plots = args.make_plots
	out_img_type = args.out_img_type
	lolli_method = args.lolli_method
	drop_neg_diffs = args.drop_neg_diffs
	drop_vars_of_uncert_signif = args.drop_vars_of_uncert_signif
	verbose = args.verbose
	out_suffix = args.out_suffix


	
		
	
	
	# =============== External annotation / resources ==============
	# Get ENST to ENSG, HGNC and Uniprot Ids mapping
	res_map = ResourceMapping()
	
	
	
	
	
	


	
	# ---- Gene selection default parameters ----
	enst_ids = None 
	gene_list = None

	signif_genes = {}
	signif_transcripts = {}
	# --

	
	

	#geneset = 'All_prioritised_oncogenes'
	#geneset = 'Shortlist_prioritised_oncogenes'
	#geneset = 'Manually_inspected_prioritised_oncogenes'
	#geneset = 'AM_Leukemia-genes'
	#geneset = 'CL_Leukemia-genes'
	#geneset = 'Rasopathy-genes'
	#geneset = 'Hematological-gene-examples'


	ids = get_gene_dataset(geneset)
	enst_ids = ids['enst_ids']
	gene_list = ids['gene_list']
	

	# DEBUG (bypass specified geneset)
	#enst_ids = ['ENST00000378609'] # ['ENST00000264033', 'ENST00000256196', 'ENST00000332029']
	#gene_list = ['GNB1']
	#geneset = 'debug'




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




	# === Read full ClinVar annotation ===
	if variant_annot_dataset == 'clinvar':
		print('Reading clinvar variant annotation...')
		clinvar_full_annot_df = read_clinvar_annot(drop_vars_of_uncert_signif)





	# ********************* Main RUN *********************
	annotator = OncMtrAnnotator(dataset, win_size, geneset, verbose=verbose)

	
	full_ENST_dist_dict = {}
	fet_res_dict_per_ENST = {}

	cnt = 1
	for enst_id in enst_ids:
		print('\n\n\n\n=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=\n-', enst_id)

		try:
			ens_gene = res_map.enst_to_ensg_dict[enst_id]
			hgnc_name = res_map.enst_to_hgnc_dict[enst_id]
		except:
			print('Could not find enst_id:', enst_id)
			continue


		cur_results = annotator.get_signals_distance(enst_id, ens_gene, hgnc_name, drop_neg_diffs=drop_neg_diffs, make_plots=make_plots, top_percentile=top_percentile)	
		if cur_results is None:
			continue 
		
		#print('\n\n=====>>>', cur_results, '\n\n')

		full_ENST_dist_dict[enst_id] = cur_results['distances']
		if 'FET P-value (Full)' in cur_results['fet'].keys():
			fet_res_dict_per_ENST[enst_id] = cur_results['fet']	


		cnt += 1
		if cnt % 100 == 0:
			print('Counter:', cnt)





	# Save all distance metrics to then plot distributions
	enstid_dist_df = pd.DataFrame(full_ENST_dist_dict).T
	enstid_dist_df.reset_index(inplace=True)
	enstid_dist_df.columns.values[0] = 'ENST'
	dist_out_file = annotator.dataset_res_dir + '/enstid_dist_metrics.' + dataset + '.tsv'
	enstid_dist_df.to_csv(dist_out_file, sep='\t', index=False)
	print('\n- Distance metrics saved into ' + dist_out_file)



	# Save all valid FET results into a data frame
	try:
		enstid_fet_df = pd.DataFrame(fet_res_dict_per_ENST).T
		enstid_fet_df.reset_index(inplace=True)
		enstid_fet_df.columns.values[0] = 'ENST'
		enstid_fet_df.sort_values(by='FET P-value (Full)', ascending=True, inplace=True)
		fet_out_file = annotator.res_dir + '/enstid_FET_results.' + dataset + '.' + variant_annot_dataset + '.tsv'
		enstid_fet_df.to_csv(fet_out_file, sep='\t', index=False)
		print("\n- Fisher's exact test results saved into " + fet_out_file)

	except Exception as e:
		print(e, '\n')
	
	





	# =============== PRINT RUN STATS ===============
	print('\n\n==== RUN STATS ====')
	assess_enst_num_with_known_variants(enst_ids, variant_annot_dataset)

	if keep_significant_genes_only:
		print('\n>> Significant HGNC genes:', len(signif_genes.keys()))
		print('>> Significant transcripts:', len(signif_transcripts.keys()))
	# ===============================================
