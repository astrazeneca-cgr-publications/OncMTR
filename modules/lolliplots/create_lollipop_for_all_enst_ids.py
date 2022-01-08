import pickle
import subprocess
import sys, os


def get_aa_mutations_per_transcript():

	pathog_lolli_entries = {}
	uncert_lolli_entries = {}

	lolli_entry_to_enst_id_dict = {}
	
	# Compile lists of amino-acid mutations per transcript
	cnt = 0
	with open('../../data/clinvar/vcf_GRCh38/clinvar.pathogenic_subset.snpeff_annotation.coding.tsv') as fh:
		for line in fh:
			line = line.rstrip()
			gene_name, ens_gene, enst_id, hgvs_p, clinvar_annot = line.split('\t')

			#if gene_name != 'TP53':
			#	continue
			#else:
			#	print('TP53 hit')

			if enst_id != 'ENST00000378609':
				continue
			
			# get unique ID per transcript: either Uniprot_ID (if it exist) or gene name
			if enst_id not in enst_to_uniprot_ids_dict:
				entry_id = gene_name
			else:
				entry_id = enst_to_uniprot_ids_dict[enst_id] #uniprot_kb ID
			lolli_entry_to_enst_id_dict[entry_id] = enst_id
			

			if clinvar_annot == 'Pathogenic' or clinvar_annot == 'Likely_pathogenic':
				if entry_id not in pathog_lolli_entries:
					pathog_lolli_entries[entry_id] = [hgvs_p]
				else:
					pathog_lolli_entries[entry_id].append(hgvs_p)
			elif clinvar_annot == 'Uncertain_significance':
				if entry_id not in uncert_lolli_entries:
					uncert_lolli_entries[entry_id] = [hgvs_p]
				else:
					uncert_lolli_entries[entry_id].append(hgvs_p)
	

	return pathog_lolli_entries, uncert_lolli_entries, lolli_entry_to_enst_id_dict



			
			
def compile_lolliplot_cmd(entry_id, enst_id, hgvs_p_list, out_dir):

	base_dir = '../../utils/lollipops-v1.5.1-linux64/'
	cmd = base_dir + "lollipops -legend -f=" + base_dir + "/arial.ttf -dpi=300"

	# Entry ids without matching UniprotKB/Swiss-Prot ID are passed on with their gene name.
	# In that case I should ommit the '-U' option prior to the gene name.
	cmd += ' -o=' + out_dir + '/' + enst_id + '.png '
	if enst_id not in enst_to_uniprot_ids_dict:
		cmd += ' ' + entry_id + ' '
	else:
		cmd += ' -U ' + entry_id + ' '
	cmd += ' '.join(hgvs_p_list)
	print(cmd)
	sys.exit()

	return cmd
	

def make_lolliplots_for_transcript_list(lolli_entries, lolli_entry_to_enst_id_dict, out_dir):

	fh_err = open('ambiguous_lolliplot_entries.err', 'w')

	cnt = 0
	for entry_id, hgvs_p_list in lolli_entries.items():


		if len(hgvs_p_list) == 0:
			return -1

		hgvs_p_list = list(set(hgvs_p_list))
		#print(entry_id)
		#print(hgvs_p_list)

		enst_id = lolli_entry_to_enst_id_dict[entry_id]
	
		cmd = compile_lolliplot_cmd(entry_id, enst_id, hgvs_p_list, out_dir)
	
		if entry_id != 'P62873':
			continue

		try:
			p = subprocess.Popen(cmd, shell=True,
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE)
		except:
			print('Error in lolliplot cmd')
			sys.exit()

		stdout, stderr = p.communicate()
		#stdout = str(stdout, "utf-8")
		#print('STDOUT:\n', stdout)
		stderr = str(stderr, "utf-8")
		print('\n\n-------\nSTDERR:\n', cmd, '\n', entry_id, '\n', stderr)

		if 'Uniprot returned' in stderr:
			fh_err.write('\n--------\n' + stderr + '\n')
	
		cnt += 1
		if cnt % 100 == 0:
			print(cnt)

	fh_err.close()



if __name__ == '__main__':

	plots_dir = '../../out/lollipop-plots'
	pathogenic_plots = plots_dir + '/pathogenic'
	uncert_signif_plots = plots_dir + '/uncertain_signif'

	if not os.path.exists(plots_dir):
		os.makedirs(plots_dir)
	if not os.path.exists(pathogenic_plots):
		os.makedirs(pathogenic_plots)
	if not os.path.exists(uncert_signif_plots):
		os.makedirs(uncert_signif_plots)


	with open('../../data/ensembl/ENST_to_Uniprot_dict.pkl', 'rb') as fh:
		enst_to_uniprot_ids_dict = pickle.load(fh)


	pathog_lolli_entries, uncert_lolli_entries, lolli_entry_to_enst_id_dict = get_aa_mutations_per_transcript()

	# create lollipop plots for pathogenic variants
	make_lolliplots_for_transcript_list(pathog_lolli_entries, lolli_entry_to_enst_id_dict, pathogenic_plots)
