import sys, os
import pandas as pd
from io import StringIO



def process_clinvar_into_data_frame():

	input_file = 'clinvar.vcf'
	processed_file_str = ''

	cnt = 0
	with open(input_file) as fh:
		for line in fh:
			if line.startswith('#'):
				continue		
			line = line.rstrip()


			vals = line.split('\t')

			header = '\t'.join(vals[:5])


			info_dict = {}
			info_vals = vals[7].split(';')
			for val in info_vals:
				tmp_key, tmp_val = val.split('=')
				
				# process MC field
				mc_val_list = []
				if tmp_key == 'MC':
					mc_vals = tmp_val.split(',')
					for mv in mc_vals:
						_, tmp_conseq = mv.split('|')
						mc_val_list.append(tmp_conseq)

					info_dict[tmp_key] = ','.join(mc_val_list)
				else:
					info_dict[tmp_key] = tmp_val


			info_str = ''
			fields_to_retain = ['RS', 'CLNSIG', 'ORIGIN', 'CLNHGVS', 'CLNVC', 'MC', 'CLNDN']


			for f in fields_to_retain:
				field_val = 'NA'
				if f in info_dict:
					field_val = info_dict[f]

				info_str += '\t' + field_val

			new_entry = header + info_str
				

			processed_file_str += new_entry + '\n'


			if DEBUG:
				if cnt > 1000:
					break
			cnt += 1
			if cnt % 10000 == 0:
				print(cnt)


	processed_df_str = StringIO(processed_file_str)
	full_df = pd.read_csv(processed_df_str, sep='\t', header=None, low_memory=False, keep_default_na=False)


	full_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT'] + fields_to_retain

	return full_df




if __name__ == '__main__':

	DEBUG = False

	out_dir = 'filtered-clinvar-files'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	full_clinvar_df = process_clinvar_into_data_frame()
	print('full_clinvar_df:', full_clinvar_df.shape)
	full_clinvar_df.to_csv(out_dir + '/full_clinvar_df.csv', sep='\t', index=False)	
	
	#INFO=<ID=ORIGIN,Number=.,Type=String,Description="Allele origin. One or more of the following values may be added: 
	# 0 - unknown;
	# 1 - germline;
	# 2 - somatic;
	# 4 - inherited;
	# 8 - paternal;
	# 16 - maternal;
	# 32 - de-novo;
	# 64 - biparental;
	# 128 - uniparental;
	# 256 - not-tested;
	# 512 - tested-inconclusive;
	# 1073741824 - other">


	# Save [germline] only variants
	germline_df = full_clinvar_df.loc[ full_clinvar_df.ORIGIN == '1', :]
	print('germline:', germline_df.shape)
	germline_df.to_csv(out_dir + '/germline.clinvar_df.csv', sep='\t', index=False)	
	
	# Save all [non-somatic] variants
	non_somatic_df = full_clinvar_df.loc[ full_clinvar_df.ORIGIN != '2', :]
	print('non_somatic_df:', non_somatic_df.shape)
	non_somatic_df.to_csv(out_dir + '/non_somatic.clinvar_df.csv', sep='\t', index=False)	

	# Save all somatic variants
	somatic_df = full_clinvar_df.loc[ full_clinvar_df.ORIGIN == '2', :]
	print('somatic_df:', somatic_df.shape)
	somatic_df.to_csv(out_dir + '/somatic.clinvar_df.csv', sep='\t', index=False)	
