import sys
from os import system
import pandas as pd
import pickle

tmp_out_file = 'tmp.hgnc_to_ensembl_mappings'

system("""zcat < Homo_sapiens.GRCh38.92.chr.gtf.gz | awk '{if($3 == "exon" && $1 != "MT") print $20"\t"$14"\t"$10}' | awk '!seen[$0]++ {print $0}' | sed 's/[\";]//g' > """ + tmp_out_file)


df = pd.read_csv(tmp_out_file, sep='\t')
df.columns = ['HGNC', 'ENST', 'ENSG']
print(df.head())
print(df.shape)

# save enst to hgnc mappings
enst_to_hgnc_dict = dict(zip(df['ENST'], df['HGNC']))
with open('./ENST_to_HGNC_dict.pkl', 'wb') as fout:         
	pickle.dump(enst_to_hgnc_dict, fout)



# ENST_to_ENSG_dict.pkl

enst_df = df.groupby(['HGNC'])['ENST'].apply(lambda x: ','.join(x.unique())).reset_index()
ensg_df = df.groupby(['HGNC'])['ENSG'].apply(lambda x: ','.join(x.unique())).reset_index()

print('ENST:')
print(enst_df.head(10))
print(enst_df.shape)

print('ENSG:')
print(ensg_df.head(10))
print(ensg_df.shape)

print('Gene-ENST-ENSG total df:')
total_df = enst_df.merge(ensg_df, left_on='HGNC', right_on='HGNC')
print(total_df.head(10))
print(total_df.shape)

# filter (some of the) non-gene entries, i.e. those with multiple ENSG values
total_df = total_df.loc[ ~total_df.ENSG.str.contains(',') ]

total_df.to_csv('HGNC_to_ensembl_mappings.tsv', sep='\t', index=False)
#system("rm " + tmp_out_file)

