import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper
import multiprocessing as mp
NUM_CORES = 4
RMIN = 0.0
RMAX = 1.0
NBIN = 1000
def calculate_hist_one_column(column, rmin, rmax, nbin):
	hist = np.histogram(column, range = (rmin, rmax), bins = nbin)
	return hist[0] # hist is a tuple of 2: [0] --> counts of elements in each bin, [1] --> the bounds of bins, this item always has length one longer than the first

def calculate_hist_one_genome_pos(avg_state_fn, num_chromHMM_state):
	avg_state_df = pd.read_csv(avg_state_fn, header = 0, index_col = 0, sep = '\t') # read one file of average state segmentation for one position on the genome
	avg_state_df.reset_index(inplace = True)
	avg_state_df = avg_state_df[map(lambda x: 'state_' + str(x + 1), range(num_chromHMM_state))] # get rid of the index which is just the indices of the genomic positions
	hist_df = avg_state_df.apply(lambda x: calculate_hist_one_column(x, RMIN, RMAX, NBIN), axis = 0) # apply function to each column  --> df: rows: each bins in the range that we specify, columns: each state
	return hist_df 

def calculate_hist_one_process(avg_state_dir, genome_pos_fn_list, out_dir, num_chromHMM_state, process_index):
	result_df_list = []
	for avg_state_fn in genome_pos_fn_list:
		this_pos_hist_df = calculate_hist_one_genome_pos(avg_state_fn, num_chromHMM_state)
		result_df_list.append(this_pos_hist_df) # add this position's data into the list of histograms results
	process_sum_hist_df = pd.concat(result_df_list).groupby(level = 0).sum() # sum the histogram data for all the position, so that we have the total count for this process
	save_fn = os.path.join(out_dir, 'hist_draft_' + str(process_index) + ".txt.gz")
	process_sum_hist_df.to_csv(save_fn, sep = '\t', header = True, index = False, compression = 'gzip')

def aggregate_data_all_processes(out_dir, num_chromHMM_state):
	processes_fn_list = map(lambda x: os.path.join(out_dir, 'hist_draft_' + str(x) +  ".txt.gz"), range(NUM_CORES))
	hist_df_list = map(lambda x: pd.read_csv(x, header = 0, sep = '\t'), processes_fn_list) # list of df of results from each process
	sum_hist_df = pd.concat(hist_df_list).groupby(level = 0).sum() # sum the df from all processes
	sum_hist_df['all_state'] = sum_hist_df.sum(axis = 1) # row sum --> data for all states
	count_save_fn = os.path.join(out_dir, 'count_hist.txt.gz')
	sum_hist_df.to_csv(count_save_fn, sep = '\t', header = True, index = False, compression = 'gzip')
	total_pos = sum_hist_df.sum(axis = 0) # column sum --> for eahc state, we got the total number of positions in each state. It should be the case that this array has all values being equal
	print "Num unique total_pos (should be 1): " ,len(np.unique(total_pos))
	assert len(np.unique(total_pos)) == 2, "The number of genomic pos in each state is not consistent" # the all_state will be the sum of all the other states
	sum_hist_df = sum_hist_df.div(total_pos, axis = 1) # divide each row bey total_pos
	prop_save_fn = os.path.join(out_dir, 'prop_hist.txt.gz')
	sum_hist_df.to_csv(prop_save_fn, sep = '\t', header = True, index = False, compression = 'gzip')
	return


def calculate_hist_parallel(avg_state_dir, out_dir, num_chromHMM_state):
	genome_pos_fn_list = glob.glob(avg_state_dir + "/*_avg_all_ct.txt.gz")
	partition_genPos_list = helper.partition_file_list(genome_pos_fn_list, NUM_CORES) # partiton that list of files so that we can give them to processes
	processes = [mp.Process(target = calculate_hist_one_process, args = (avg_state_dir, partition_genPos_list[i], out_dir, num_chromHMM_state, i)) for i in range(NUM_CORES)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print "Process " + str(i) + " is finished!"
	print "Done calculating histogram for each proceses"
	# now aggregate data from all processes
	aggregate_data_all_processes(out_dir, num_chromHMM_state)
	print "Done calculating the data for all cell type within this cell group"

def main():
	if len(sys.argv) != 4:
		usage()
	avg_state_dir = sys.argv[1]
	helper.check_dir_exist(avg_state_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[3])
	print "Done getting command line arguments"
	calculate_hist_parallel(avg_state_dir, out_dir, num_chromHMM_state)


def usage():
	print "python count_histogram_pos_prob_one_cg.py"
	print "avg_state_dir: folder representative_data/avg_segmentation of the cell group"
	print "out_dir: where output of the representative chromHMM tracks should be stored. Inside this folder there would be two folders: avg_segmentation, segmentation_igv_format"
	print "num_chromHMM_state"
	exit(1)
main()