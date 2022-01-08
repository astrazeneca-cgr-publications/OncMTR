import pandas as pd
import pickle



class ResourceMapping:

	def __init__(self):
		# Read base directory
		self.base_dir = ( 
		pd.read_csv('../../.conf', index_col = 0, sep = '\t')
			.loc['base_dir'].absolute_path
		)
		self.data_dir = self.base_dir + '/data/ensembl'
				
		# Read HGNC to ENST and ENSG mappings data frame
		self.hgnc_df = pd.read_csv(self.data_dir + '/HGNC_to_ensembl_mappings.tsv', sep='\t')
		
		# Initialise ENST/ENSG/HGNC mappings
		self.read_enst_to_ensg_and_hgnc_maps()




	def read_enst_to_ensg_and_hgnc_maps(self):

		with open(self.data_dir + '/ENST_to_ENSG_dict.pkl', 'rb') as fh:
			self.enst_to_ensg_dict = pickle.load(fh)

		with open(self.data_dir + '/ENST_to_HGNC_dict.pkl', 'rb') as fh:
			self.enst_to_hgnc_dict = pickle.load(fh)

		with open(self.data_dir + '/ENST_to_Uniprot_dict.pkl', 'rb') as fh:
			self.enst_to_uniprot_dict = pickle.load(fh)

		with open(self.data_dir + '/ENSG_to_Uniprot_dict.pkl', 'rb') as fh:
			self.ensg_to_uniprot_dict = pickle.load(fh)
			
			
			
	def get_hgnc_to_enst_mapping(self, gene_list):

		tmp_df = self.hgnc_df.loc[ self.hgnc_df.HGNC.isin(gene_list) ]
		#print(tmp_df.head())

		enst_ids = tmp_df['ENST'].tolist()
		enst_ids = [v.split(',') for v in enst_ids]
		enst_ids = [enst_ids[i][j] for i in range(len(enst_ids)) for j in range(len(enst_ids[i])) ]

		return enst_ids
