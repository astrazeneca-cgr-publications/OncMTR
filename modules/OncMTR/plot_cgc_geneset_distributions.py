import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import sys, os
from stats_util import is_outlier




def plot_metric_distr(metric, out_dir, outliers_filter):
	
	cols_to_keep = ['ENSG', metric]

	tmp_cgc_tier1_df = cgc_tier1_metrics_df[cols_to_keep].dropna()
	tmp_cgc_tier2_df = cgc_tier2_metrics_df[cols_to_keep].dropna()
	tmp_rest_df = rest_metrics_df[cols_to_keep].dropna()

	
	uniq_cgc_tier1_df = tmp_cgc_tier1_df.groupby(['ENSG']).max()
	uniq_cgc_tier2_df = tmp_cgc_tier2_df.groupby(['ENSG']).max()
	uniq_rest_df = tmp_rest_df.groupby(['ENSG']).max()
	#print('Tier1:', uniq_cgc_tier1_df.shape)
	#print('Tier2:', uniq_cgc_tier2_df.shape)
	#print('Rest:', uniq_rest_df.shape)
	
	random_sample_df = pd.concat([uniq_cgc_tier1_df, uniq_cgc_tier2_df, uniq_rest_df], axis=0).sample(n=uniq_cgc_tier1_df.shape[0])
	#print('Random sample:', random_sample_df.shape)
	

	uniq_cgc_tier1 = uniq_cgc_tier1_df[metric].copy()	
	uniq_cgc_tier2 = uniq_cgc_tier2_df[metric].copy()	
	random_sample = random_sample_df[metric].copy()
	uniq_rest = uniq_rest_df[metric].copy()


	if outliers_filter == 'no_outliers':
		uniq_cgc_tier1 = uniq_cgc_tier1[ ~is_outlier(uniq_cgc_tier1, 3.5)]
		uniq_cgc_tier2 = uniq_cgc_tier2[ ~is_outlier(uniq_cgc_tier2, 3.5)]
		random_sample = random_sample[ ~is_outlier(random_sample, 3.5)]
		uniq_rest = uniq_rest[ ~is_outlier(uniq_rest, 3.5)]
	elif outliers_filter == 'only_outliers':
		uniq_cgc_tier1 = uniq_cgc_tier1[ is_outlier(uniq_cgc_tier1, 3.5)]
		uniq_cgc_tier2 = uniq_cgc_tier2[ is_outlier(uniq_cgc_tier2, 3.5)]
		random_sample = random_sample[ is_outlier(random_sample, 3.5)]
		uniq_rest = uniq_rest[ is_outlier(uniq_rest, 3.5)]


	sns.set(font_scale=1.5)
	sns.set_style("ticks")
	fig, ax = plt.subplots(figsize=(15, 8))
	sns.kdeplot(uniq_cgc_tier1, 
				color = '#31a354', 
				#hist_kws={'edgecolor':'black'},
				shade=True,
				linewidth= 2,
				label='CGC_tier1 ('+str(len(uniq_cgc_tier1))+')' )
 

	sns.kdeplot(uniq_cgc_tier2, 
				color = '#3182bd', 
				#hist_kws={'edgecolor':'black'},
				shade=True,
				linewidth= 2,
				label='CGC_tier2 ('+str(len(uniq_cgc_tier2))+')')
				
	sns.kdeplot(random_sample, 
				color = '#feb24c', 
				#hist_kws={'edgecolor':'black'},
				shade=True,
				linewidth= 2,
				label='random-sample ('+str(len(random_sample))+')')
				
	sns.kdeplot(uniq_rest, 
				color = 'black', 
				#hist_kws={'edgecolor':'black'},
				shade=True,
				linewidth= 2,
				label='Rest ('+str(len(uniq_rest))+')')
				
	# Calculate Mann-Whittney U p-values between distributions of different gene sets 
	mannu_tier1_vs_rest_pval = '{:.2e}'.format(float(mannwhitneyu(uniq_cgc_tier1, uniq_rest).pvalue))
	mannu_tier2_vs_rest_pval = '{:.2e}'.format(float(mannwhitneyu(uniq_cgc_tier2, uniq_rest).pvalue))
	mannu_tier1_vs_tier2_pval = '{:.2e}'.format(float(mannwhitneyu(uniq_cgc_tier1, uniq_cgc_tier2).pvalue))
	mannu_tier1_vs_random_pval = '{:.2e}'.format(float(mannwhitneyu(uniq_cgc_tier1, random_sample).pvalue))


	print('uniq_cgc_tier1:', np.median(uniq_cgc_tier1))
	print('uniq_cgc_tier2:', np.median(uniq_cgc_tier2))
	print('uniq_rest:', np.median(uniq_rest))

	
	plt.title('Mann-Whittney U tests:\n' + '\u2022 CGC_tier1 vs Rest = ' + mannu_tier1_vs_rest_pval + '\n' \
			+ '\u2022 CGC_tier2 vs Rest = ' + mannu_tier2_vs_rest_pval + '\n'\
			+ '\u2022 CGC_tier1 vs random_sample = ' + mannu_tier1_vs_random_pval + '\n'\
			+ '\u2022 CGC_tier1 vs CGC_tier2 = ' + mannu_tier1_vs_tier2_pval, \
			fontsize=20)
	plt.legend(fontsize=20)

	out_dir += '/' + outliers_filter
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	out_fig_file = metric + '_distr.' + outliers_filter + '.pdf'
	fig.savefig(out_dir + '/' + out_fig_file, bbox_inches='tight')
	
	




	
if __name__ == '__main__':
	
	input_dir = 'OncMTR-gnomad-win31' 
	#outliers_filter = 'no_outliers' # no_outliers, only_outliers, None
	#outliers_filter = 'only_outliers'
	#outliers_filter = 'all'


	out_dir = '../../out/' + input_dir + '/Distance-distribution-figs'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)	


	metrics_df = pd.read_csv('../../out/' + input_dir + '/enstid_dist_metrics.gnomad.tsv.simple_set', index_col=None, sep='\t')
	print(metrics_df.head())
	print(metrics_df.shape)


	cgc_tier1_df = pd.read_csv('../../data/Cancer_Gene_Census/CGC_tier1.tsv.processed.csv', header=None)
	cgc_tier1_genes = cgc_tier1_df.iloc[ :, 1].values.tolist()
	print(len(cgc_tier1_genes))


	cgc_tier2_df = pd.read_csv('../../data/Cancer_Gene_Census/CGC_tier2.tsv.processed.csv', header=None)
	cgc_tier2_genes = cgc_tier2_df.iloc[ :, 1].values.tolist()
	print(len(cgc_tier2_genes))



	cgc_tier1_metrics_df = metrics_df.loc[ metrics_df.ENSG.isin(cgc_tier1_genes), : ]
	cgc_tier2_metrics_df = metrics_df.loc[ metrics_df.ENSG.isin(cgc_tier2_genes), : ]
	rest_metrics_df = metrics_df.loc[ ~metrics_df.ENSG.isin(cgc_tier1_genes + cgc_tier2_genes), : ]
	print(cgc_tier1_metrics_df.shape)
	print(cgc_tier2_metrics_df.shape)
	print(rest_metrics_df.shape)
	

	for outliers_filter in ['all', 'no_outliers', 'only_outliers']:
		print('>>', outliers_filter)
	
		for metric in metrics_df.columns:
			if metric in ['ENST', 'ENSG', 'HGNC', 'transcr_len']:
				continue
			else:
				print('\n', metric)
			
			plot_metric_distr(metric, out_dir, outliers_filter)
