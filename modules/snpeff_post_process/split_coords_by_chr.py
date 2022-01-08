import sys
import os


out_dir = '../../out/Genomic_to_protein_coords'

file_handles = {}
for chrom in list(range(1,23)) + ['X', 'Y']:
	file_handles[str(chrom)] = open(out_dir + '/genome_coords_to_hgvs_p.chr' + str(chrom) + '.tsv', 'w')
	print(chrom)

cnt = 0
with open(out_dir + '/genome_coords_to_hgvs_p.tsv') as fh:
	for line in fh:
		vals = line.rstrip().split("\t")
		cur_chrom = vals[0]

		file_handles[cur_chrom].write('\t'.join(vals) + '\n')	
