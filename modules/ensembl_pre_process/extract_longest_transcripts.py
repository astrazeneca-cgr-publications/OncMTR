import pandas as pd
import sys

data_dir = '../../data/ensembl'

df = pd.read_csv(data_dir + '/hsa_exons.GRCh38.92.tsv', sep='\t', header=None, low_memory=False)
# 1  11868  12227  ENST00000456328  ENSG00000223972	1
df.columns = ['chr', 'start', 'end', 'transcript_id', 'gene_id', 'exon_num']
df['length'] = df.end - df.start + 1
print(df.head())
print(df.shape)


# > Get total length of exons per transcript
sum_df = df.groupby(['transcript_id', 'gene_id'])[['transcript_id', 'gene_id', 'length']].sum()
sum_df.reset_index(inplace=True)
print(sum_df.head())
print(sum_df.shape)

print('----------------------------')

## DEBUG
#print(sum_df.loc[ sum_df.gene_id == 'ENSG00000012048'])
#tt = sum_df.loc[ sum_df.gene_id == 'ENSG00000012048']
#print(tt.loc[ tt.groupby(['gene_id'])['length'].idxmax() ] )
#print( tt.groupby(['gene_id'])['length'].idxmax())


# > Keep transcript with maximum length
print(sum_df.info())
longest_df = sum_df.loc[ sum_df.groupby(['gene_id'])['length'].idxmax() ]
print(longest_df.head())
print(longest_df.shape)

#longest_df.drop('length', axis=1, inplace=True)
longest_df.columns = ['ENSG', 'longest_ENST', 'length']
print(longest_df.head())
print(longest_df.shape)

longest_df.to_csv(data_dir + '/hsa_exons.GRCh38.92.longest_transcript.tsv', sep='\t', index=False)



""" 
# Deprecated
longest_transcripts = max_exon_df['transcript_id']
longest_df = df.loc[ df.transcript_id.isin(longest_transcripts), :]
print(longest_df.head())
print(longest_df.shape)

# > Save full df with longest transcripts
longest_df.to_csv(data_dir + '/hsa_exons.GRCh38.92.longest_transcript.tsv', sep='\t', index=False)


# > Create BED file with longest transcripts
longest_bed = longest_df[['chr', 'start', 'end']].copy()
longest_bed.to_csv(data_dir + '/hsa_exons.GRCh38.92.longest_transcript.bed', sep='\t', header=False, index=False)
"""
