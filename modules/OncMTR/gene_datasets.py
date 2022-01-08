import pandas as pd


def get_gene_dataset(geneset):
  
	ids = {'enst_ids': None, 'gene_list': None}	

	# ===================== ENST ID lists =====================
	if geneset == 'All_prioritised_oncogenes':
		enst_ids = pd.read_csv('../../data/cancer-hotspots/enst_of_onc_interest.txt', header=None)
		ids['enst_ids'] = enst_ids.iloc[:, 0].tolist()


	elif geneset == 'Shortlist_prioritised_oncogenes':
		enst_ids = pd.read_csv('../../data/cancer-hotspots/shortlist_prioritised.enst_of_onc_interest.txt', header=None)
		ids['enst_ids'] = enst_ids.iloc[:, 0].tolist()

	
	elif geneset == 'Manually_inspected_prioritised_oncogenes':
		enst_ids = pd.read_csv('../../data/cancer-hotspots/manually_inspected_prioritised.enst_of_onc_interest.txt', header=None)
		ids['enst_ids'] = enst_ids.iloc[:, 0].tolist()

	
	elif geneset == 'AM_Leukemia-genes':
		# Acute Myleoid Leukemia
		#gene_list = ['FLT3', 'NPM1', 'DNMT3A', 'IDH2', 'IDH1', 'TET2', 'RUNX1', 'NRAS', 'CEBPA', 'PTPN11', 'U2AF1', 'KRAS', 'SMC1A', 'TP53', 'WT1', 'KIT']
		enst_ids = pd.read_csv('../../data/cancer-hotspots/AML_transcripts.txt', header=None)
		ids['enst_ids'] = enst_ids.iloc[:, 0].tolist()

	elif geneset == 'CL_Leukemia-genes':
		# Chronic lymphocytic leukemia (CLL)
		#gene_list = ['SF3B1', 'NOTCH1', 'POT1', 'ATM', 'CHD2', 'PCLO']
		enst_ids = pd.read_csv('../../data/cancer-hotspots/CLL_transcripts.txt', header=None)
		ids['enst_ids'] = enst_ids.iloc[:, 0].tolist()


	
	# ===================== Gene lists (HGNC) =====================
	elif geneset == 'Rasopathy-genes': 
		gene_list = ['BRAF', 'CBL', 'HRAS', 'KRAS', 'LZTR1', 'MAP2K1', 'MAP2K2', 'NF1', 'NRAS', 'PPP1CB', 'PTPN11', 'RAF1', 'RIT1', 'SHOC2', 'SOS1', 'SOS2', 'SPRED1', 'A2ML1', 'RASA2']
		ids['gene_list'] = gene_list

	
	elif geneset == 'Hematological-gene-examples':
		gene_list = ['UA2F1', 'SRFS2']  # To be Fixed: non-HGNC names provided
		ids['gene_list'] = gene_list
	


	return ids
