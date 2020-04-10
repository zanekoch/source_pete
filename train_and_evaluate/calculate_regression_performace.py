import numpy as np 
import pandas as pd 
import os
import sys
import helper
import glob
def check_results_one_genome_window(reg_res_fn, segment_fn, check_ct, num_chromHMM_state):
	segment_df = pd.read_csv(segment_fn, sep = '\t', header = 0)
	truth_segment = segment_df[check_ct] # get only the segmentation data of the check ct
	reg_res_df = pd.read_csv(reg_res_fn, sep = '\t', header = 0, index_col = 0)
	predict_segment = reg_res_df.idxmax(axis = 1) # list of 'state_<index>' of the state assigned at each position
	predict_segment = predict_segment.apply(lambda x: 'E' + x.split('_')[1]) # convert from 'state_18' to E18
	assert len(truth_segment) == len(predict_segment), "The number of genomic positions in prediction file: " + reg_res_fn + " DOES NOT MATCh the number of genomic positions in segment_fn: " + segment_fn
	# declare a confusion matrix
	count_df = pd.DataFrame({'truth' : list(truth_segment), 'pred' : list(predict_segment)}) 
	count_df = count_df.groupby(['truth', 'pred']).size() # count the number times that we observe a particular pair of true-predicted states
	count_df = count_df.reset_index() # this is necessary to turn the results from previous line to a proper data frame
	confusion_df = count_df.pivot(index = 'truth', columns = 'pred', values = 0) # convert from a count table to a confusion matrix, entries are count of times we observe any pair of true-predicted states
	confusion_df = confusion_df.fillna(0) # if we dont observe any pairs of true-predicted states, we will fill it with a count of 0
	state_annot_list = map(lambda x: 'E' + str(x+1), range(num_chromHMM_state))
	# check to fill out missing rows first with zeros
	(nrow, ncol) = confusion_df.shape
	if nrow != num_chromHMM_state: 
		not_included_states = np.setdiff1d(np.array(state_annot_list), confusion_df.index)
		for state in not_included_states: # this state is not present in this region in the truthful segmentation
			confusion_df.loc[state,:] = np.zeros(ncol)
	# check to fill out missing columns with zeros
	(nrow, ncol) = confusion_df.shape
	if ncol != num_chromHMM_state:
		not_included_states = np.setdiff1d(np.array(state_annot_list), confusion_df.columns)
		for state in not_included_states:
			confusion_df[state] = np.zeros(nrow)
	# now we have a complete num_chromHMM_state * num_chromHMM_state matrix
	confusion_df = confusion_df.reindex(state_annot_list) # reorder the rows so that they are fron state 1 to state num_chromHMM_state
	confusion_df = confusion_df[state_annot_list] # reorder the columns so that they are fron state 1 to state num_chromHMM_state
	# rename indices and columns
	confusion_df.index = map(lambda x: 'true_' + x, confusion_df.index)
	confusion_df.columns = map(lambda x: 'pred_' + x, confusion_df.columns) # confusion matrix: indicates, for each state of truth (original) segmentation, the percentage that the segmentation gets redistributed to the prediction states
	return confusion_df

def check_results_all_genome_windows(reg_res_dir, all_ct_segment_dir, check_ct, output_fn, num_chromHMM_state):
	truth_fn_list = glob.glob(all_ct_segment_dir + "/*_combined_segment.bed.gz")
	gen_loc_list = map(lambda x: (x.split('/')[-1]).split('_combined_segment.bed.gz')[0], truth_fn_list)
	reg_predict_fn_list = map(lambda x: os.path.join(reg_res_dir, x + '_avg_pred.txt.gz'), gen_loc_list) # get the list of respective prediction files, ordered similarly to truth_fn_list
	colnames_confusion_df = map(lambda x: 'pred_E' + str(x+1), range(num_chromHMM_state))
	rownames_confusion_df = map(lambda x: 'true_E' + str(x+1), range(num_chromHMM_state))
	confusion_df = pd.DataFrame(0, columns = colnames_confusion_df, index = rownames_confusion_df) # confusion matrix: indicates, for each state of truth (original) segmentation, the percentage that the segmentation gets redistributed to the prediction states
	for pos_index in range(len(reg_predict_fn_list)):
		reg_res_fn = reg_predict_fn_list[pos_index]
		segment_fn = truth_fn_list[pos_index]
		this_pos_confusion_df = check_results_one_genome_window(reg_res_fn, segment_fn, check_ct, num_chromHMM_state)
		confusion_df = confusion_df.add(this_pos_confusion_df) # add the confusion matrices from multiple files corresponding to multiple regions on the genome
	# now we have calculated all the number of times that any pairs of true_predicted states appeared, we will convert all our data into the frequencies form
	num_truth_segment_per_state = confusion_df.sum(axis = 1) # row sum
	confusion_df = confusion_df.divide(num_truth_segment_per_state, axis = 0) # divide each row by num_truth_segment_per_state --> frequencies of assigning each true state to a predicted state
	# write results to file
	confusion_df.to_csv(output_fn, sep = '\t', header = True, index = True)

def main():
	if len(sys.argv) != 6:
		usage()
	reg_res_dir = sys.argv[1]
	helper.check_dir_exist(reg_res_dir)
	all_ct_segment_dir = sys.argv[2]
	helper.check_dir_exist(all_ct_segment_dir)
	output_fn = sys.argv[3]
	helper.create_folder_for_file(output_fn)
	check_ct = sys.argv[4]
	try:
		num_chromHMM_state = int(sys.argv[5])
	except:
		print "num_chromHMM_state IS NOT VALID: " + sys.argv[5]
		usage()
	print "Done getting command line arguments"
	check_results_all_genome_windows(reg_res_dir, all_ct_segment_dir, check_ct, output_fn, num_chromHMM_state)

def usage():
	print "python calculate_regression_performance.py"
	print "reg_res_dir: directory of regression results for different positions: rows: positions, columns: states --> prob of each position being in each state"
	print "all_ct_segment_dir: fn where the truth segmentation of all cell types are stored"
	print "output_fn"
	print "check_ct: the ct that we are trying to calculate the performance of segmentation prediction"
	print "num_chromHMM_state: self-explanatory"
	exit(1)

main()