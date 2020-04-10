import pandas as pd 
import numpy as np 
import os
import sys
import glob
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
NUM_CORES = 2
RMIN = 0.0
RMAX = 1.0
NBIN = 1000

def read_one_hist_df(hist_fn):
	return(pd.read_csv(hist_fn, sep = '\t', header = 0))

def average_histogram_across_all_ct(all_ct_hist_dir, out_dir, num_chromHMM_state):
	ct_hist_fn_list = glob.glob(all_ct_hist_dir + "/E*/prop_hist.txt.gz")
	ct_hist_df_list = map(read_one_hist_df, ct_hist_fn_list)
	avg_hist_df = pd.concat(ct_hist_df_list).groupby(level = 0).mean()
	save_fn = os.path.join(out_dir, 'avgerage_posterior_histogram.txt.gz')
	avg_hist_df.to_csv(save_fn, sep = '\t', header = True, index = False, compression = 'gzip')

def main():
	if len(sys.argv) != 4:
		usage()
	all_ct_hist_dir = sys.argv[1]
	helper.check_dir_exist(all_ct_hist_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[3])
	print "Done getting command line arguments"
	average_histogram_across_all_ct(all_ct_hist_dir, out_dir, num_chromHMM_state)
	print "Done!"

def usage():
	print "python average_histogram_across_ct.py"
	print "all_ct_hist_dir: where each ct has a folder, and inside each folder is the information about histogram of that cell type"
	print "out_dir: where data of average histogram is stored"
	print "num_chromHMM_state"
	exit(1)

main()