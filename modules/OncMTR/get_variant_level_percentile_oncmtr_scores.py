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




# Get percentile OncMTR scores
if __name__ == "__main__":


	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-d", dest="dataset", help="Input dataset", choices=['gnomad', 'UKB', 'topmed'], required=True)
	parser.add_argument("-w", dest="win_size", default=31, help="Window size to calculate OncMTR", required=True)
	parser.add_argument("-s", dest="out_suffix", default='', help="Custom output folder suffix", required=False)

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	args = parser.parse_args()
	print(args, '\n')

	dataset = args.dataset
	win_size = args.win_size
	out_suffix = args.out_suffix


	# Input directory
	base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['base_out'].absolute_path
	)

	input_dir = base_out_dir + '/OncMTR-' + dataset + '-win' + str(win_size)
	if out_suffix not in ['', 'None']:
		input_dir += '-' + out_suffix
	print('Input dir:', input_dir)

	variant_level_dir = input_dir + '/variant-level-scores'


	print('Reading (tmp) full vcf file...')
	# debug: set nrows=100000
	#full_df = pd.read_csv(variant_level_dir + '/OncMTR-' +  dataset + '-win' + str(win_size) + '.full.tmp', sep='\t', nrows=100000)
	full_df = pd.read_csv(variant_level_dir + '/OncMTR-' +  dataset + '-win' + str(win_size) + '.full.tmp', sep='\t', low_memory=False)
	print(full_df.head())
	print(full_df.shape)


			
	full_df_nona = full_df[['chrom', 'pos', 'ref', 'alt', 'enst_id', 'OncMTR']].copy().dropna()
	print('No NAs:', full_df_nona.shape)
	full_df_nona.index = full_df_nona['chrom'].astype(str) + '_' + full_df_nona['pos'].astype(str) + '_' + full_df_nona['ref'] + '_' + full_df_nona['alt'] + '_' + full_df_nona['enst_id']
	full_df_nona.drop(['chrom', 'pos', 'ref', 'alt', 'enst_id'], axis=1, inplace=True)

	print(full_df_nona.head())
	print(full_df_nona.tail())

	print('\n- Calculating percentiles...')
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
	full_pct_df.rename(columns={"OncMTR_x": "OncMTR", "OncMTR_y": "OncMTR_percentile"}, inplace=True)


	full_pct_df = full_pct_df[['chrom', 'pos', 'ref', 'alt', 'enst_id', 'gene_id', 'protein_pos', 'OncMTR', 'OncMTR_percentile', 'MTR_FDR', 'MTR_ab_FDR']]

	print('\n- Sorting results by chrom and position...')
	full_pct_df.sort_values(by=['chrom', 'pos'], inplace=True, ascending=True)
	
	
	out_file = variant_level_dir + '/OncMTR-' +  dataset + '-win' + str(win_size) + '.full.vcf'
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

