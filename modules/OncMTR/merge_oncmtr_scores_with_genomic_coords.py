import pandas as pd
import numpy as np
import multiprocessing as mp
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import subprocess
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
		


def merge_oncmtr_with_coords_for_chrom(chrom):

	print('Reading genomic coords to hgvs_p for chrom', chrom, '...\n')
	df = pd.read_csv(base_out_dir + '/Genomic_to_protein_coords/genome_coords_to_hgvs_p.chr' + str(chrom) + '.tsv', sep='\t', header=None, low_memory=False)

	df.columns = ['chrom', 'pos', 'ref', 'alt', 'ENST_ID', 'protein_pos']
	#print(df.head())

	print('Sorting genomic coords df for chrom', chrom, '...\n')
	genomic_sorted_df = df.sort_values(by=['chrom', 'pos', 'ENST_ID', 'protein_pos'], ascending=True)

	# Edit protein_pos so that start codon (ATG) is at index 0
	genomic_sorted_df.protein_pos -= 1
	
	
	print('Merging OncMTR scores with genomic coords for chrom', chrom, '...\n')
	full_chrom_df = pd.DataFrame()
	cnt = 1
	for enst_id in genomic_sorted_df.ENST_ID.unique():
		
		if not enst_id in oncmtr_df.ENST_ID.values:
			continue
			
		cur_oncmtr_entry = oncmtr_df.loc[ oncmtr_df.ENST_ID == enst_id ]
		cur_mtr = cur_oncmtr_entry['OncMTR'].values[0].split(', ')
		cur_gene = cur_oncmtr_entry['Gene'].values[0]
		cur_oncmtr_df = pd.DataFrame({'gene_id': cur_gene, 'OncMTR': cur_mtr})	
		
		
		cur_genom_df = genomic_sorted_df.loc[ genomic_sorted_df.ENST_ID == enst_id]
		if cur_genom_df.shape[0] != (cur_oncmtr_df.shape[0] * 9):
			pass
			#print('[Warning] Inconsistent OncMTR and genomic coords list sizes for', enst_id)
			#print('genomic coords:', cur_genom_df.shape[0], ' -  mtr (x9):', cur_oncmtr_df.shape[0] * 9, '\n')
			

		cur_merged_df = cur_genom_df.merge(cur_oncmtr_df, left_on='protein_pos', right_index=True)

		# get FDR scores
		cur_mtr_fdr = np.array(eval(mtr_fdr_df.loc[ mtr_fdr_df.ENST_ID == enst_id, 'FDR'].values[0].replace('nan', 'np.nan')))
		cur_mtr_ab_fdr = np.array(eval(mtr_fdr_ab_df.loc[ mtr_fdr_ab_df.ENST_ID == enst_id, 'FDR'].values[0].replace('nan', 'np.nan')))
		cur_fdrs_df = pd.DataFrame({'MTR_FDR': cur_mtr_fdr, 'MTR_ab_FDR': cur_mtr_ab_fdr})	

		# merge with rest
		cur_merged_df = cur_merged_df.merge(cur_fdrs_df, left_on='protein_pos', right_index=True)

		
		if full_chrom_df.shape[0] > 0:
			full_chrom_df = pd.concat([full_chrom_df, cur_merged_df], axis=0)
		else:
			full_chrom_df = cur_merged_df

		#print(cnt, enst_id)
		if cnt % 200 == 0:
			print('\n\n>>>', chrom, '-', cnt, '\n')
			#break
		cnt += 1
		
			

	full_chrom_df.reset_index(inplace=True, drop=True)
	full_chrom_df.rename(columns={'ENST_ID': 'enst_id'}, inplace=True)
	full_chrom_df = full_chrom_df[['chrom', 'pos', 'ref', 'alt', 'enst_id', 'gene_id', 'protein_pos', 'OncMTR', 'MTR_FDR', 'MTR_ab_FDR']]

	# remove start codon coords
	full_chrom_df = full_chrom_df.loc[ full_chrom_df.protein_pos != 0]
	full_chrom_df.to_csv(out_dir + '/OncMTR-' + dataset + '-win' + str(win_size) + '.chr' + str(chrom) + '.tsv', index=False, header=False, sep='\t')





# Merge OncMTR scores genomic-to-protein coordinate mappings
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



	# Input OncMTR directory
	base_out_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['base_out'].absolute_path
	)

	input_dir = base_out_dir + '/OncMTR-' + dataset + '-win' + str(win_size)
	if out_suffix not in ['', 'None']:
		input_dir += '-' + out_suffix
	print('Input dir:', input_dir)

	# OncMTR scores
	input_mtr_full_file = 'All_OncMTR_scores-' + dataset + '-win' + str(win_size) + '.tsv'
	oncmtr_df = pd.read_csv(input_dir + '/' + input_mtr_full_file, sep='\t')


	# Read MTR / MTR-AB FDR scores from base dirs 
	mtr_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size)
	mtr_fdr_df = pd.read_csv(mtr_dir + '/All_FDR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t')

	mtr_ab_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size) + '-AB-filtered'
	mtr_fdr_ab_df = pd.read_csv(mtr_ab_dir + '/All_FDR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t')


	# Output directory
	out_dir = input_dir + '/variant-level-scores'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	print('Output dir:', out_dir)


	

	# Debug run
	#merge_oncmtr_with_coords_for_chrom(21)
	#sys.exit()
	all_chroms = list(range(1,23)) + ['X', 'Y']
	

	# ======= Main Run =======
	pool = mp.Pool(16)
	result = pool.map(merge_oncmtr_with_coords_for_chrom, all_chroms)

	pool.close()
	pool.join()
    # ========================
	


	print('> Merge all files in a single VCF ...')
	vcf_columns = ['chrom', 'pos', 'ref', 'alt', 'enst_id', 'gene_id', 'protein_pos', 'OncMTR', 'MTR_FDR', 'MTR_ab_FDR']
	with open(out_dir + '/vcf_header.tmp', 'w') as out_fh:
		out_fh.write('\t'.join(vcf_columns) + '\n')


	# add header for full VCF
	input_vcfs_str = out_dir + '/vcf_header.tmp '
	for chrom in all_chroms:
		print(chrom)
		input_vcfs_str += out_dir + '/OncMTR-' + dataset + '-win' + str(win_size) + '.chr' + str(chrom) + '.tsv '
	
	full_vcf_out_file = out_dir + '/OncMTR-' + dataset + '-win' + str(win_size) + '.full.tmp'
	
	cmd = 'cat ' + input_vcfs_str + ' > ' + full_vcf_out_file
	#print(cmd)
	proc = subprocess.Popen(cmd, shell=True,
							stdout=subprocess.PIPE,
							stderr=subprocess.PIPE)
	(output, err) = proc.communicate()  

	# wait for cmd to finish
	p_status = proc.wait()


	# cleanup
	os.system('rm ' + out_dir + '/vcf_header.tmp')
