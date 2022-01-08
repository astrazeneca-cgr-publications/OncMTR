import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import numpy as np
import math
import pickle
from scipy.stats import binom_test



class MtrClass():

	def __init__(self, win_size=31, min_variants_thres=3, out_dir='MTR-out'):
		self.win_size = int(win_size) #31
		self.min_variants_thres = int(min_variants_thres) #3
		self.offset = int(self.win_size / 2)

		self.out_dir = out_dir
		self.mtr_scores_dir = out_dir + '/scores'
		self.mtr_plots_dir = out_dir + '/pdfs' 	
		
		self.fdr_scores_dir = out_dir + '/fdr_scores'
		self.fdr_plots_dir = out_dir + '/fdr_pdfs'
		
		if not os.path.exists(self.out_dir):
			os.makedirs(self.out_dir)

		if not os.path.exists(self.mtr_scores_dir):
			os.makedirs(self.mtr_scores_dir)
		
		if not os.path.exists(self.fdr_scores_dir):
			os.makedirs(self.fdr_scores_dir)
		
		if not os.path.exists(self.mtr_plots_dir):
			os.makedirs(self.mtr_plots_dir)
		
		if not os.path.exists(self.fdr_plots_dir):
			os.makedirs(self.fdr_plots_dir)
			
		self.valid_data = True # will become False if missing data are detected
	

	def mtr_formula(self, obs_mis, obs_syn, exp_mis, exp_syn, base):

		#if obs_mis + obs_syn == 0:
		#	print("obs_mis + obs_syn = 0")
		#if exp_mis + exp_syn == 0:
		#	print("exp_mis + exp_syn = 0")


		if (obs_mis + obs_syn) < self.min_variants_thres:
			#print('(obs_mis + obs_syn) < self.min_variants_thres -- Returning Nan...')
			return np.nan, np.nan, np.nan

		#print('> exp_mis:', exp_mis)
		#print('> exp_syn:', exp_syn)
		#print('> obs_mis:', obs_mis)
		#print('> obs_syn:', obs_syn)

		try:
			mtr_obs = (obs_mis / (obs_mis + obs_syn))
			mtr_exp = (exp_mis / (exp_mis + exp_syn))
			mtr = mtr_obs / mtr_exp

			return np.round(mtr_obs, 5), np.round(mtr_exp, 5), np.round(mtr, 5)

		except:
			return np.nan, np.nan, np.nan



	def fdr_formula(self, obs_mis, obs_syn, exp_mis, exp_syn, base):

		if (obs_mis + obs_syn) < self.min_variants_thres:
			return np.nan

		try:
			x = obs_mis
			n = obs_mis + obs_syn
			p = exp_mis / (exp_mis + exp_syn)
			fdr = binom_test(x, n, p, alternative='two-sided')
			
			return np.round(fdr, 5)
		except:
			return np.nan

	


	def roundup(self, x):
		return math.ceil(x)




	def read_exp_obs_data(self, enst_id, exp_syn, exp_mis, obs_syn, obs_mis):

		self.is_reverse = False
		try:
			if exp_mis[0] > exp_mis[-1]:
				self.is_reverse = True
		except:
			# > Abort: expected transcript has zero hits - insufficient data/not well-annotated transcript
			print('[Error] ' + enst_id + ' - Non valid data for MTR calculation...', exp_mis)
			self.valid_data = False
			return -1


		exp_syn_df = pd.DataFrame(pd.Series(exp_syn).value_counts())
		exp_syn_df.columns = ['exp_syn']
		#print(exp_syn_df.head())

		exp_mis_df = pd.DataFrame(pd.Series(exp_mis).value_counts())
		exp_mis_df.columns = ['exp_mis']
		#print(exp_mis_df.head())

		obs_syn_df = pd.DataFrame(pd.Series(obs_syn).value_counts())
		obs_syn_df.columns = ['obs_syn']
		#print(obs_syn_df.head())

		obs_mis_df = pd.DataFrame(pd.Series(obs_mis).value_counts())
		obs_mis_df.columns = ['obs_mis']
		#print(obs_mis_df.head())


		# merge hits from exp. syn/mis and obs. syn/mis into a single data frame
		full_df = pd.merge(obs_syn_df, obs_mis_df, left_index=True, right_index=True, how='outer')
		full_df = pd.merge(full_df, exp_syn_df, left_index=True, right_index=True, how='outer')
		full_df = pd.merge(full_df, exp_mis_df, left_index=True, right_index=True, how='outer')
		#print(full_df.head())
		#print(full_df.tail())



		# fill-in missing coding sequence positions with zeros
		max_position = max(full_df.index)
		valid_rows = full_df.index.tolist()
		all_rows = list(range(1, max_position + 1))

		zero_rows = list(set(all_rows) - set(valid_rows))
		zero_df = pd.DataFrame(0, index=zero_rows, columns=['zero_pos'])

		full_df = pd.merge(full_df, zero_df, left_index=True, right_index=True, how='outer')
		full_df.fillna(0, inplace=True)
		full_df.drop(['zero_pos'], axis=1, inplace=True)
		#print('\n', full_df.head(60))
		#print(full_df.shape)

		# count hits per codon	
		full_df['codon'] = full_df.index.values / 3
		full_df['codon'] = full_df['codon'].apply(self.roundup)
		#print('\n', full_df.head(20))
		#print(full_df.shape)

		self.codon_df = full_df.groupby('codon').sum()
		#print('\nCodon_df:', self.codon_df.head(20))
		#print(self.codon_df.shape)



	def calc_full_mtr_signal(self, enst_id):
		valid_mtr = True

		self.mtr_vector = []
		self.mtr_stats = []
		
		self.fdr_vector = []
		self.fdr_stats = []

		max_position = max(self.codon_df.index)
		for base in range(1, max_position + 1):
			
			left = base - self.offset
			right = base + self.offset

			#print('base:', base, ', left:', left, ', right:', right)
			#print(self.codon_df.loc[left:right])

			
			collapsed_df = self.codon_df.loc[left:right].sum(axis=0)
			#print(collapsed_df)
			obs_syn = collapsed_df.loc['obs_syn']
			obs_mis = collapsed_df.loc['obs_mis']
			exp_syn = collapsed_df.loc['exp_syn']
			exp_mis = collapsed_df.loc['exp_mis']

			
			mtr_obs, mtr_exp, tmp_mtr = self.mtr_formula(obs_mis, obs_syn, exp_mis, exp_syn, base)
			tmp_fdr = self.fdr_formula(obs_mis, obs_syn, exp_mis, exp_syn, base)

			#print(tmp_mtr)
			self.mtr_vector.append(tmp_mtr)
			self.mtr_stats.append([base, obs_mis, obs_syn, exp_mis, exp_syn, mtr_obs, mtr_exp, tmp_mtr])
			
			self.fdr_vector.append(tmp_fdr)
			self.fdr_stats.append([base, obs_mis, obs_syn, exp_mis, exp_syn, tmp_fdr])


		# Check if MTR vector contains only NaNs
		if np.isnan(self.mtr_vector).all() or all(v==0 for v in self.mtr_vector):
			valid_mtr = False
			print(enst_id, ': Invalid MTR vector - all NANs or all zeros')
			return valid_mtr

		self.mtr_stats = pd.DataFrame(self.mtr_stats, columns = ['position', 'obs_mis', 'obs_syn', 'exp_mis', 'exp_syn', 'mtr_obs', 'mtr_exp', 'mtr'])
		self.fdr_stats = pd.DataFrame(self.fdr_stats, columns = ['position', 'obs_mis', 'obs_syn', 'exp_mis', 'exp_syn', 'fdr'])

		return valid_mtr



	def save_mtr(self, enst_id, gene_id='ENSG'):

		with open(self.mtr_scores_dir + '/' + gene_id + '.' + enst_id + '.MTR_score.tsv', 'w') as fout:
			mtr_list_str = ', '.join([str(e) for e in self.mtr_vector])

			fout.write(gene_id+'\t'+enst_id+'\t'+mtr_list_str+'\n')
			#pickle.dump(self.mtr_vector, fout)

	def save_fdr(self, enst_id, gene_id='ENSG'):

		with open(self.fdr_scores_dir + '/' + gene_id + '.' + enst_id + '.fdr_score.tsv', 'w') as fout:
			fdr_list_str = ', '.join([str(e) for e in self.fdr_vector])

			fout.write(gene_id+'\t'+enst_id+'\t'+fdr_list_str+'\n')
			#pickle.dump(self.fdr_vector, fout)
			

	def plot_mtr(self, enst_id, gene_id='ENSG'):
		#if self.is_reverse:
		#	self.mtr_vector = self.mtr_vector[::-1]
		try:
			max_mtr = np.nanmax(self.mtr_vector)
		except:
			pass

		# > Do not plot transcripts with all values NaN
		if math.isnan(max_mtr): 
			return -1  

		try:
			fig = plt.figure(figsize=(8, 5))
			plt.plot(self.mtr_vector)
			plt.ylim(-0.1, max_mtr + 0.1)
			plt.grid(linewidth=0.5, linestyle='-')

			fig.savefig(self.mtr_plots_dir + '/' + gene_id + '.' + enst_id + '_Plot.pdf')	
			plt.close()
		except:
			return -1
		

	def plot_fdr(self, enst_id, gene_id='ENSG'):
		#if self.is_reverse:
		#	self.fdr_vector = self.fdr_vector[::-1]
		try:
			max_fdr = np.nanmax(self.fdr_vector)
		except:
			pass

		# > Do not plot transcripts with all values NaN
		if math.isnan(max_fdr): 
			return -1  


		fig = plt.figure(figsize=(8, 5))
		plt.plot(self.fdr_vector)
		plt.ylim(-0.1, max_fdr + 0.1)
		plt.grid(linewidth=0.5, linestyle='-')

		fig.savefig(self.fdr_plots_dir + '/' + gene_id + '.' + enst_id + '_Plot.pdf')	
		plt.close()
		
		
		
	def run(self, exp_syn, exp_mis, obs_syn, obs_mis, enst_id, gene_id='ENSG', make_plot=True):

		self.read_exp_obs_data(enst_id, exp_syn, exp_mis, obs_syn, obs_mis)
		if not self.valid_data:
			return -1, -1

		valid_mtr = self.calc_full_mtr_signal(enst_id)
		if not valid_mtr:
			return -1, -1

		self.save_mtr(enst_id, gene_id)
		if make_plot:
			self.plot_mtr(enst_id, gene_id)

# 		self.save_fdr(enst_id, gene_id)
# 		if make_plot:
# 			self.plot_fdr(enst_id, gene_id)
		

		self.mtr_stats = self.mtr_stats.assign(enst = enst_id).assign(gene_id = gene_id)
		self.fdr_stats = self.fdr_stats.assign(enst = enst_id).assign(gene_id = gene_id)

		return self.mtr_stats, self.fdr_stats
			


if __name__ == '__main__':

	win_size = 31
	min_variants_thres = 3	

	# --------- GNB1 example ---------
	gnb1_dict = {}
	with open('gnb1_example.txt') as fh:
		for line in fh:
			line = line.rstrip()
			vals = line.split('\t')
			attr, value = vals[0], vals[1]

			if '[' in value:
				value = eval(value)

			gnb1_dict[attr] = value

	gene_id = gnb1_dict['gene_id']
	enst_id = gnb1_dict['enst_id']
	exp_syn = gnb1_dict['exp_syn']
	exp_mis = gnb1_dict['exp_mis']
	obs_syn = gnb1_dict['obs_syn']
	obs_mis = gnb1_dict['obs_mis']
	#print(gnb1_dict)
	# ---------------------------------

	mtr_obj = MtrClass(win_size=31, min_variants_thres=3, out_dir='MTR-example')
	mtr_obj.run(exp_syn, exp_mis, obs_syn, obs_mis, enst_id, gene_id)
