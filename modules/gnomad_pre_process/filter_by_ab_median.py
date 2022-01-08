import sys

input_file = sys.argv[1] #"gnomad.exomes.r2.0.2.sites.processed.snp.uniq.pass.vcf"

ab_median_cutoff = 0.30

out_file = input_file.replace('.vcf', '.ab_filtered.vcf')
fout = open(out_file, 'w')


cnt = 0
with open(input_file) as fh:
	for line in fh:
		cnt += 1
		if line.startswith('#'):
			fout.write(line)
			continue

		line = line.rstrip()
		vals = line.split('\t')

		info = vals[7]
		info_vals = info.split(';')
		ab_mean_field = [f for f in info_vals if 'AB_MEDIAN' in f]
		ab_median_val = ab_mean_field[0].split('=')[1]
		if ab_median_val == '.':
			continue
		else:
			ab_median_val = float(ab_median_val)
			#print(ab_median_val)

		if ab_median_val > ab_median_cutoff:
			fout.write(line + '\n')


		if cnt % 100000 == 0:
			print(cnt)
fout.close()
