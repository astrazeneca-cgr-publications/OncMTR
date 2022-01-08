import pickle
import sys


data_dir = '../../data/ensembl/'

#  **************    ENST ---> ENSG    **************
enst_to_ensg_dict = {}

processed_ensembl_gtf_annot_file = data_dir + '/hsa_exons.GRCh38.92.tsv'
with open(processed_ensembl_gtf_annot_file) as fh:
	for line in fh:
		vals = line.split('\t')

		enst_id = vals[3]
		ens_gene = vals[4]

		if enst_id not in enst_to_ensg_dict:
			enst_to_ensg_dict[enst_id] = ens_gene

with open(data_dir + '/ENST_to_ENSG_dict.pkl', 'wb') as fout:
	pickle.dump(enst_to_ensg_dict, fout)
print('ENST to ENSG dict size:', len(enst_to_ensg_dict.keys()))






#  **************    ENST ---> UNIPROT IDs   **************
enst_to_uniprot_dict = {}
ensg_to_uniprot_dict = {}

with open(data_dir + '/enst_to_uniprot_ids.biomart.txt') as fh:
	for line in fh:
		line = line.rstrip()
		vals = line.split('\t')
		#print(vals)

		if len(vals) < 3:
			continue

		ensg = vals[0]
		enst_id = vals[1]
		uniprot_id = vals[2]

		enst_to_uniprot_dict[enst_id] = uniprot_id
		ensg_to_uniprot_dict[ensg] = uniprot_id


with open(data_dir + '/ENST_to_Uniprot_dict.pkl', 'wb') as fout:
	pickle.dump(enst_to_uniprot_dict, fout)
print('ENST to Uniprot dict size:', len(enst_to_uniprot_dict.keys()))


with open(data_dir + '/ENSG_to_Uniprot_dict.pkl', 'wb') as fout:
	pickle.dump(ensg_to_uniprot_dict, fout)
print('ENSG to Uniprot dict size:', len(ensg_to_uniprot_dict.keys()))


