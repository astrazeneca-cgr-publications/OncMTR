import pandas as pd
import multiprocessing as mp
import subprocess
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import sys
import os



def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', '1'):
		return True
	elif v.lower() in ('no', 'false', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')




# Get percentile MTR scores
if __name__ == "__main__":


	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-d", dest="dataset", help="Input dataset", choices=['gnomad', 'UKB_gnomad', 'UKB', 'topmed'], required=True)
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


	# Input directory
	base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['base_out'].absolute_path
	)

	input_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size)
	if ab_filtered:
		input_dir += '-AB-filtered'
	if out_suffix not in ['', 'None']:
		input_dir += '-' + out_suffix
	print('Input dir:', input_dir)

	variant_level_dir = input_dir + '/variant-level-scores'


	print('Reading (tmp) full vcf file...')
	# debug: set nrows=100000
	#full_df = pd.read_csv(variant_level_dir + '/MTR-' +  dataset + '-win' + str(win_size) + '.full.tmp', sep='\t', nrows=100000)
	full_df = pd.read_csv(variant_level_dir + '/MTR-' +  dataset + '-win' + str(win_size) + '.full.tmp', sep='\t', low_memory=False)
	print(full_df.info())

	# cleanup redundant header rows
	full_df = full_df.loc[ full_df.chrom != 'chrom', :]
	print(full_df.head())
	print(full_df.shape)
	
	"""
	print(full_df.tail())
	full_df.iloc[ 10000:19999, 0] = '11'
	full_df.iloc[ 20000:29999, 0] = '21'
	full_df.iloc[ 30000:39999, 0] = '10'
	full_df.iloc[ 40000:49999, 0] = '1'
	full_df.iloc[ 50000:59999, 0] = '2'
	full_df.iloc[ 60000:69999, 0] = '9'
	full_df.iloc[ 70000:79999, 0] = 'Y'
	full_df.iloc[ 80000:89999, 0] = 'X'
	full_df.iloc[ 90000:99999, 0] = '4'
	print(full_df.tail())
	print(full_df.chrom.unique())
	"""
			
	full_df_nona = full_df[['chrom', 'pos', 'ref', 'alt', 'enst_id', 'MTR']].copy().dropna()
	print('No NAs:', full_df_nona.shape)
	full_df_nona.index = full_df_nona['chrom'].astype(str) + '_' + full_df_nona['pos'].astype(str) + '_' + full_df_nona['ref'] + '_' + full_df_nona['alt'] + '_' + full_df_nona['enst_id']
	full_df_nona.drop(['chrom', 'pos', 'ref', 'alt', 'enst_id'], axis=1, inplace=True)

	print(full_df_nona.head())
	print(full_df_nona.tail())
	print(full_df_nona.shape)

	print('\n- Calculating percentiles...')
	full_df_nona['MTR'] = full_df_nona['MTR'].astype(float)
	pct_df = full_df_nona.rank(pct=True)
	print(pct_df.head())
	print(pct_df.tail())

	#cleanup
	del full_df_nona

	full_df.index = full_df['chrom'].astype(str) + '_' + full_df['pos'].astype(str) + '_' + full_df['ref'] + '_' + full_df['alt'] + '_' + full_df['enst_id']
	full_pct_df = full_df.merge(pct_df, left_index=True, right_index=True)

	#cleanup
	del pct_df

	full_pct_df.reset_index(drop=True, inplace=True)
	full_pct_df.rename(columns={"MTR_x": "MTR", "MTR_y": "MTR_percentile"}, inplace=True)


	full_pct_df = full_pct_df[['chrom', 'pos', 'ref', 'alt', 'enst_id', 'gene_id', 'protein_pos', 'MTR', 'MTR_percentile', 'FDR']]

	print('\n- Sorting results by chrom and position...')
	full_pct_df['chrom'] = full_pct_df['chrom'].apply(lambda x: x if x in ['X', 'Y'] else int(x))
	full_pct_df['pos'] = full_pct_df['pos'].astype(int)
	full_pct_df.sort_values(by=['chrom', 'pos'], inplace=True, ascending=True)
	
	
	out_file = variant_level_dir + '/MTR-' +  dataset + '-win' + str(win_size) + '.full.vcf'
	print('\n- Saving to:', out_file)
	full_pct_df.to_csv(out_file, sep='\t', index=False)


	print('> bgzip ...')
	cmd = 'bgzip -f ' + out_file
	print(cmd)
	proc = subprocess.Popen(cmd, shell=True,
		    stdout=subprocess.PIPE,
		    stderr=subprocess.PIPE)
	(output, err) = proc.communicate()
	proc.wait()


	print('> tabix ...')
	cmd = 'tabix -f -p vcf ' + out_file + '.gz'
	print(cmd)
	proc = subprocess.Popen(cmd, shell=True,
		    stdout=subprocess.PIPE,
		    stderr=subprocess.PIPE)
	(output, err) = proc.communicate()
	proc.wait()

