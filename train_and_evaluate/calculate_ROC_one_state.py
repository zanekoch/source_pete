import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper

def count_num_segment_in_state_and_not(pred_df, num_chromHMM_state, ct):
	# ct is used to extract the true segmentation data 
	# for each state calculate the number of positions, within this pred_df, that are not actually assigned to each of the state
	# pandas series with indices E<state_index>
	# one series will report the number of segments that are assigned to the state in the entire file
	# one series will report the number of segments that are not assigned to the state in the entire file
	# first count the number of times that each state is assigned across this file pred_df
	state_segment_count = pred_df[ct].value_counts() # index: E<state_index>
	num_total_segment = pred_df.shape[0] # number of rows --> num positions
	for state_index in range(num_chromHMM_state):
		state_str = 'E' + str(state_index + 1)
		if state_str not in state_segment_count.index:
			state_segment_count[state_str] = 0
	
	count_segment_not_state = num_total_segment - state_segment_count
	return state_segment_count, count_segment_not_state


def get_pred_segment_one_fn(segment_fn, pred_fn, reverse_lower_bound_list, ct, num_chromHMM_state, result_colnames, result_rownames):
	pred_df = pd.read_csv(pred_fn, header = 0, index_col = 0, sep = '\t')
	segment_df = pd.read_csv(segment_fn, sep = '\t', header = 0)
	pred_df[ct] = segment_df[ct] # get the segmentation of the ct that we are looking at
	# pred_df now has 'state_1' ... 'state_18' --> posterior probabilities, based on model prediction, and 'ct' --> true segmentation based on chromHMM
	# count_segment_not_state : pandas series such that: --> E<state_index, one_based> --> number of positions that are actually (truth)  not assigned to each of the state, based on the data of pred_df
	(count_segment_in_state, count_segment_not_state) = count_num_segment_in_state_and_not(pred_df, num_chromHMM_state, ct) # count the number of segments that we observe across this file, solely based on the true segmentation (regardless of prosterior score thresholds), that are in each of the state and that are also not in each of the states. Each entry in the series is indexed 'E'<state_index> --> total numbers across the entire file
	# CREATE A RESULT_DF
	result_df = pd.DataFrame(columns = result_colnames, index = result_rownames) # create a df where each state has 4 fields
	for bin_index in range(len(reverse_lower_bound_list)): # for each bin of score, loop through each of the state and get the places where the posterior probabilities for this state are below the upper threshold or above the lower threshold for this bin
		# get lower bound of score of bins in descending order
		low_score_bin = reverse_lower_bound_list[bin_index]
		this_bin_report = pd.Series([]) # an empty pandas series that we will report later with data related to all states
		for state_index in range(num_chromHMM_state): # zero based indices
			one_based_state_index = state_index + 1
			state_str = 'state_' + str(one_based_state_index)
			this_bin_state_df = pred_df[pred_df[state_str] >= low_score_bin] # get positions where the posterior score of this state is at least the low_score_bin
			this_state_true_pos = np.sum(this_bin_state_df[ct] == ('E' + str(one_based_state_index)))
			num_seg_pred_state_above_thres = this_bin_state_df.shape[0] # number of positions that are assigned  to this state because the posterior probabilities is higher than the bin thresholds
			this_state_false_pos = num_seg_pred_state_above_thres - this_state_true_pos
			num_seg_not_state = count_segment_not_state['E' + str(one_based_state_index)]
			this_bin_report['true_pos_num_S' + str(one_based_state_index)] = this_state_true_pos
			this_bin_report['false_pos_num_S' + str(one_based_state_index)] = this_state_false_pos
			this_bin_report['num_seg_true_S' + str(one_based_state_index)] = count_segment_in_state['E' + str(one_based_state_index)] # number of segments in this file that are assigned to this state, by truth, and we do not care about the score thresholds of posterior prediction here
			this_bin_report['num_seg_true_not_S' + str(one_based_state_index)] = count_segment_not_state['E' + str(one_based_state_index)] # number of segments in this file that are not assigned, by truth, to the state, we do not care about posterior predictions here
			this_bin_report['num_seg_above_thres_S' + str(one_based_state_index)] = num_seg_pred_state_above_thres # count the number of segments that score above the threshold for this state
		# now we filled in enough report of all the state for this posterior score bin
		result_df.loc[result_rownames[bin_index]] =  this_bin_report
	return result_df

def get_tp_fp_data_all_regions(true_segment_dir, pred_dir, reverse_lower_bound_list, ct, num_chromHMM_state):
	pred_fn_list = glob.glob(pred_dir + "/*.txt.gz")
	chrom_region_list = map(lambda x: (x.split('/')[-1]).split('.')[0][:-9], pred_fn_list) # convert from 'chr7_10_avg_pred.txt.gz' to chr7_10
	true_segment_fn_list = map(lambda x: os.path.join(true_segment_dir, x + '_combined_segment.bed.gz'), chrom_region_list) # get a list of file name ex: /path/to/true/segment/chr2_1_combined_segment.bed.gz 
	# format the result table
	# column names
	result_colnames = []
	for state_index in range(num_chromHMM_state):
		this_state_colnames = ['true_pos_num_S' + str(state_index + 1), 'false_pos_num_S' + str(state_index + 1), 'num_seg_true_S' + str(state_index + 1), 'num_seg_true_not_S' + str(state_index + 1), 'num_seg_above_thres_S' + str(state_index + 1)]
		# true_pos_num_S1: within in each bin, the number of positions that are predicted rightly as S1
		# num_seg_pred_S1: number of positions that have S1 scores in the bin --> which is then predicted as S1
		# false_pos_num_S1: number of positions that have S1 scores in the bin but are aactually not assigned to the state S1
		result_colnames += this_state_colnames
	# row names
	result_rownames = (map(lambda x: "score_greater_" + str(x), reverse_lower_bound_list))
	# declare the results
	total_result_df = pd.DataFrame(0, columns = result_colnames, index = result_rownames)
	# call function get_pred_segment_one_fn
	for region_index, region in enumerate(chrom_region_list):
		this_region_segment_fn = true_segment_fn_list[region_index]
		this_region_pred_fn = pred_fn_list[region_index]
		this_region_df = get_pred_segment_one_fn(this_region_segment_fn, this_region_pred_fn, reverse_lower_bound_list, ct, num_chromHMM_state, result_colnames, result_rownames)
		total_result_df = total_result_df.add(this_region_df)
	return total_result_df

def get_score_bins (num_score_bins):
	bin_size = 1.0 / float(num_score_bins)
	upper_bound_list = map(lambda x: bin_size * x, range(1, num_score_bins + 1)) # if num_score_bins: 10 then it should be [0.1, 0.2, .., 1]
	lower_bound_list = map(lambda x: bin_size * x, range(num_score_bins)) # [0, 0.1, ..., 0.9]
	lower_bound_list.reverse() # reversed now
	return (lower_bound_list, upper_bound_list)

def calculate_tpr_fpr(total_tp_fp_df, num_chromHMM_state, save_fn):
	# total_tp_fp_df is the output of function get_tp_fp_data_all_regions
	for state_index in range(num_chromHMM_state):
		one_based_state_index = str(state_index + 1)
		total_tp_fp_df['tpr_S' + (one_based_state_index)] = total_tp_fp_df['true_pos_num_S' + one_based_state_index] / total_tp_fp_df['num_seg_true_S' + one_based_state_index]
		total_tp_fp_df['fpr_S' + one_based_state_index] = total_tp_fp_df['false_pos_num_S' + one_based_state_index] / total_tp_fp_df['num_seg_true_not_S' + one_based_state_index]
		total_tp_fp_df['precision_S' + one_based_state_index] = total_tp_fp_df['true_pos_num_S' + one_based_state_index] / total_tp_fp_df['num_seg_above_thres_S' + one_based_state_index]
		total_tp_fp_df['recall_S' + one_based_state_index] = total_tp_fp_df['tpr_S' + one_based_state_index]
	total_tp_fp_df.to_csv(save_fn, sep = '\t', header = True, index = True, compression = 'gzip')
		
def main():
	if len(sys.argv) != 7:
		usage()
	pred_dir = sys.argv[1]
	helper.check_dir_exist(pred_dir)
	true_segment_dir = sys.argv[2]
	helper.check_dir_exist(true_segment_dir)
	ct = sys.argv[3]
	outDir = sys.argv[4]
	helper.make_dir(outDir)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[5])
	num_score_bins = helper.get_command_line_integer(sys.argv[6])
	print "Done getting command line arguments"
	# first get the upper bounds for the score bins
	(reverse_lower_bound_list, upper_bound_score_list) = get_score_bins(num_score_bins)
	print "Get the bounds of posterior probabilities that we will set for each of the bin"
	# get the count of true positives and false positives, etc. across all regions in the genome
	total_tp_fp_df = get_tp_fp_data_all_regions(true_segment_dir, pred_dir, reverse_lower_bound_list, ct, num_chromHMM_state)
	print "Done processing all the files corresponding to all the regions in the genome"
	# calculate tpr and fpr values for each of the state
	save_fn = os.path.join(outDir, 'tpr_fpr_all_states.txt.gz')
	calculate_tpr_fpr(total_tp_fp_df, num_chromHMM_state, save_fn)
	print "Done calculating true positive rates and false positive rates in all bins"

def usage():
	print "python calculate_ROC_one_state.py"
	print "pred_dir: dir of the different positions on the genome and the probabilities that each of positions are in any particular states."
	print "true_segment_dir: dir of the segmentations for different regions on the genome"
	print "ct : cell type whose segmentation we are trying to evaluate prediction for"
	print "outDir: directory where the data of ROC curves for each of the states are stored"
	print "num_chromHMM_state"
	print "num_score_bins: number of bins to divide [0,1] into"
	exit(1)

main()