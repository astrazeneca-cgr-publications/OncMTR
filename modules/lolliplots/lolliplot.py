import pickle
import subprocess
import sys, os
import pandas as pd
import pathlib



def onc_hotspots_lolliplot_for_transcript(enst, uniprot_id, out_dir, out_img_type, method='lollipops', variant_annot_dataset='clinvar', verbose=False):

	df = pd.read_csv('../../data/cancer-hotspots/transcripts_hotspots_oncology.txt', sep='\t')
	df = df.loc[ df['Transcript_ID'] == enst, :]
	df.reset_index(inplace=True)
	df.dropna(inplace=True)
	if verbose:
		print(df.head())

	
	hgvs_p_list = df['HGVSp'].tolist()
	cur_hgnc = df.loc[0, 'Hugo_Symbol']

	
	if method == 'genvisr':
		genvisr_data = pd.DataFrame({'gene': cur_hgnc, 'amino_acid_change': hgvs_p_list, 'transcript_name': enst})
		#print(genvisr_data.head())

		tmp_dir = './tmp'
		if not os.path.exists(tmp_dir):
			os.makedirs(tmp_dir)
		tmp_data_file = tmp_dir + '/' + enst + '.aa.tsv'

		genvisr_data.to_csv(tmp_data_file, index=False, sep='\t')

		out_img_file = run_genvisr(tmp_data_file, enst, cur_hgnc, out_dir, variant_annot_dataset, out_img_type)

	elif method == 'lollipops':		

		#print(uniprot_id, enst, cur_hgnc)
		#print(hgvs_p_list)

		out_img_file = run_lollipops_cmd(uniprot_id, enst, cur_hgnc, hgvs_p_list, out_dir, variant_annot_dataset)

	return out_img_file







def clinvar_lolliplot_for_transcript(df, enst, uniprot_id, out_dir, out_img_type, method='lollipops', variant_annot_dataset='clinvar', verbose=False):

	uncertain_signif_color = '#969696'

	# Set different color to variants of 'Uncertain_significance'
	df.reset_index(inplace=True)

	if verbose:
		print(df.head())
		print(df.shape)
		print(df.columns)
		print(enst)
	

	hgvs_p_list = ('p.' + df['HGVS_P']).tolist()
	#uniprot_id = df['UNIPROT_ID'].unique().tolist()[0]
	cur_hgnc = df.loc[0, 'HGNC']

	if verbose:
		print('hgvs_p_list:', hgvs_p_list)
		print('uniprot_id:', uniprot_id)
		print('cur_hgnc:', cur_hgnc)
	
	
	if method == 'genvisr':
		genvisr_data = pd.DataFrame({'gene': cur_hgnc, 'amino_acid_change': hgvs_p_list, 'transcript_name': enst, 'SnpEff_Annotation': df['SNPEFF_ANNOT'].tolist()})
		#print(genvisr_data.head())

		tmp_dir = './tmp'
		if not os.path.exists(tmp_dir):
			os.makedirs(tmp_dir)
		tmp_data_file = tmp_dir + '/' + enst + '.aa.tsv'

		genvisr_data.to_csv(tmp_data_file, index=False, sep='\t')

		out_img_file = run_genvisr(tmp_data_file, enst, cur_hgnc, out_dir, variant_annot_dataset, out_img_type)

	elif method == 'lollipops':		
		# Add special colouring to variants of uncertain significance
		df['HGVS_P_LOLLI'] = df['HGVS_P'].copy()
		df.loc[ df.CLNSIG.str.contains('Uncertain_significance'), 'HGVS_P_LOLLI'] += '#969696'
		hgvs_p_lolli_list = df['HGVS_P_LOLLI'].tolist()

		out_img_file = run_lollipops_cmd(uniprot_id, enst, cur_hgnc, hgvs_p_lolli_list, out_dir, variant_annot_dataset)
		
	return out_img_file
	

			





def run_lollipops_cmd(uniprot_id, enst, hgnc, hgvs_p_list, out_dir, variant_annot_dataset, dpi=600, verbose=False):

	base_dir = '../../utils/lollipops-v1.5.1-linux64/'
	#cmd = base_dir + "lollipops -legend -f=" + base_dir + "/arial.ttf -dpi=300"  # with legend (deprecated as it gets truncated)
	cmd = base_dir + "lollipops -f=" + base_dir + "/arial.ttf "

	# Entry ids without matching UniprotKB/Swiss-Prot ID are passed on with their gene name.
	# In that case I should ommit the '-U' option prior to the gene name.


	out_img_file = out_dir + '/' + hgnc + '_' + enst + '.' + variant_annot_dataset + '.lolliplot.png'

	cmd += ' -o=' + out_img_file

	cmd += ' -U ' + uniprot_id
	cmd += ' -dpi=' + str(dpi) + ' '
	cmd += ' '.join(hgvs_p_list)


	if len(hgvs_p_list) == 0:
		return 

	try:
		p = subprocess.Popen(cmd, shell=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE)
	except Exception as e:
		print('Error in lollipops cmd', e)
		sys.exit()

	stdout, stderr = p.communicate()
	#stdout = str(stdout, "utf-8")
	#print('STDOUT:\n', stdout)
	stderr = str(stderr, "utf-8")
	if verbose:
		print('\n\n-------\nLolliplot output:\n', cmd, '\n', uniprot_id, '\n', stderr)

	return out_img_file






def run_genvisr(tmp_data_file, enst, hgnc, out_dir, variant_annot_dataset, out_img_type):
	
	attempt = 1
	max_attempts = 5
	success = False

	out_img_file = out_dir + '/' + hgnc + '_' + enst + '.' + variant_annot_dataset + '.GenVisR.' + out_img_type
	
	cur_dir = str(pathlib.Path(__file__).parent.absolute())
	print("Lolliplot - current working dir:", cur_dir)

	cmd = "Rscript " + cur_dir + "/run_genvisr.R " + tmp_data_file + " " + enst + " " + hgnc + " " + out_img_file + " " + out_img_type
	print("\n ", cmd)

	while not success and attempt <= max_attempts:
		# TODO -- this does not currently capture errors in the R script
		# need read print output from R script to see if it has been succesful
		try:
			print("Calling BioMart API - attempt " + str(attempt) + "...")
			os.system(cmd)
			success = True
			print('Success.')
		except:	
			print('Current attemtp failed. Retrying...')
			attempt += 1
			sleep(1)
			
	os.remove(tmp_data_file)
	
	if not success:
		return -1
	else:
		return out_img_file








if __name__ == '__main__':

	method = 'lollipops' # 'lollipops' # genvisr
	variant_annot_dataset = 'clinvar' #'oncology' # clinvar
	out_dir = '.'
	out_img_type = 'png'
	enst = 'ENST00000378609'


	with open('../../data/ensembl/ENST_to_Uniprot_dict.pkl' ,'rb') as fh:
		enst_to_uniprot_dict = pickle.load(fh)

	try:
			if enst == 'ENST00000288602':
				uniprot_id = 'P15056'
			else:
				uniprot_id = enst_to_uniprot_dict[enst]
	except Expection as e:
		print('[Error] Cannot find Uniprot id for transcript', enst, '\n')
		print(e, '\n\n')
		sys.exit()



		
	if variant_annot_dataset == 'oncology':

		onc_hotspots_lolliplot_for_transcript(enst, uniprot_id, out_dir, out_img_type, method=method, variant_annot_dataset=variant_annot_dataset)


	elif variant_annot_dataset == 'clinvar':
		# get subset of HGVS_p variants for transcript
		full_annot_df = pd.read_csv('../../data/clinvar/vcf_GRCh38/filtered-clinvar-files/clinvar_pathogenic.mis_syn.non-somatic.tsv', sep='\t', low_memory=False)
		df = full_annot_df.loc[ full_annot_df.ENST == enst, :].copy()

		clinvar_lolliplot_for_transcript(df, enst, uniprot_id, out_dir, out_img_type, method=method, variant_annot_dataset=variant_annot_dataset)
