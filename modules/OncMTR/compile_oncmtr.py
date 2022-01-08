import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import numpy as np
import sys
import os



def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', '1'):
        return True
    elif v.lower() in ('no', 'false', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def calc_onc_mtr(x, debug=False, debug_enst_id='ENST00000378609'):
    
        global cnt

        #print(x)
        gene = x['Gene']
        enst_id = x['ENST_ID']
        
        if debug and enst_id != debug_enst_id:
                return
        
        mtr = x['MTR'].replace('nan', 'np.nan')
        mtr_ab = x['MTR_ab'].replace('nan', 'np.nan')
        
        mtr = np.array( eval(mtr) )
        mtr_ab = np.array( eval(mtr_ab) )
        #print(type(mtr))
        #print(type(mtr_ab))
        
        cur_oncmtr = mtr_ab - mtr
        #print(cur_oncmtr)
        
        
        cur_oncmtr_str = ', '.join([str(np.round(x, 5)) for x in cur_oncmtr])
        new_oncmtr_entry = '\t'.join([gene, enst_id, cur_oncmtr_str])
        #print(new_oncmtr_entry)
        

        cnt += 1
        if cnt % 1000 == 0:
            print(cnt)

        if debug and enst_id == debug_enst_id:
            plt.figure(figsize=(10,4))
            plt.plot(mtr)

            #plt.figure(figsize=(10,4))
            plt.plot(mtr_ab)

            plt.figure(figsize=(10,4))
            plt.plot(cur_oncmtr)
            sys.exit()
        
        return new_oncmtr_entry
        




if __name__ == '__main__':

	# Read input arguments
	parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("-d", dest="dataset", help="Input dataset", choices=['gnomad', 'UKB', 'topmed'], required=True)
	parser.add_argument("-w", dest="win_size", default=31, help="Window size to calculate MTR", required=True)
	parser.add_argument("-s", dest="out_suffix", default='', help="Custom output folder suffix", required=False)

	if len(sys.argv)==1:
		parser.print_help(sys.stderr)
		sys.exit(1)


	args = parser.parse_args()
	print(args, '\n')

	dataset = args.dataset
	win_size = args.win_size
	out_suffix = args.out_suffix

	#dataset = 'gnomad'
	#win_size = 31
	#out_suffix = ''


	# Output directory
	base_out_dir = ( 
	pd.read_csv('../../.conf', index_col = 0, sep = '\t')
		.loc['base_out'].absolute_path
	)

	# MTR dir
	print('\n- Reading MTR scores...')
	mtr_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size)
	if out_suffix not in ['', 'None']:
		mtr_dir += '-' + out_suffix
	print('MTR dir:', mtr_dir)
	mtr_df = pd.read_csv(mtr_dir + '/All_MTR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t')
	print(mtr_df.head())
	print(mtr_df.shape)


	# MTR AB-filtered dir
	print('\n- Reading MTR AB-filtered scores...')
	mtr_ab_dir = base_out_dir + '/MTR-' + dataset + '-win' + str(win_size) + '-AB-filtered'
	if out_suffix not in ['', 'None']:
		mtr_ab_dir += '-' + out_suffix
	print('MTR AB-filtered dir:', mtr_ab_dir)
	mtr_ab_df = pd.read_csv(mtr_ab_dir + '/All_MTR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t')
	mtr_ab_df.columns = ['Gene', 'ENST_ID', 'MTR_ab']
	print(mtr_ab_df.head())
	print(mtr_ab_df.shape)



	# OncMTR output dir
	oncmtr_dir = base_out_dir + '/OncMTR-' + dataset + '-win' + str(win_size)
	if out_suffix not in ['', 'None']:
		oncmtr_dir += '-' + out_suffix
	if not os.path.exists(oncmtr_dir):
		os.makedirs(oncmtr_dir)


	
	print('\n- Creating merged MTR / MTR AB-filtered df...')
	oncmtr_df = mtr_df.merge(mtr_ab_df, left_on=['Gene', 'ENST_ID'], right_on=['Gene', 'ENST_ID'])
	print(oncmtr_df.head())
	print(oncmtr_df.shape)

	oncmtr_df.to_csv(oncmtr_dir + '/OncMTR_signal_pairs_df-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t', index=False)


	cnt = 0
	# Calculate OncMTR
	print('\n- Calculating OncMTR per transcript...')
	res = oncmtr_df.apply(lambda x: calc_onc_mtr(x), axis=1)


	# Save OncMTR scores into a data frame
	oncmtr_full_str = '\n'.join(res)
	final_oncmtr_df = pd.read_csv(StringIO(oncmtr_full_str), sep="\t", header=None)
	final_oncmtr_df.columns = ['Gene', 'ENST_ID', 'OncMTR']
	final_oncmtr_df.head()
	final_oncmtr_df.to_csv(oncmtr_dir + '/All_OncMTR_scores-' + dataset + '-win' + str(win_size) + '.tsv', sep='\t', index=False)

