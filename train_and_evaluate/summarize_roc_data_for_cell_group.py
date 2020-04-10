import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper
import pickle
ROC_TERMS_LIST = ['tpr', 'fpr', 'precision', 'recall']
def get_columns_names_to_fetch(num_chromHMM_model):
	columns_to_get = []	
	for state_index in range(num_chromHMM_model):
		colnames_this_state = map(lambda x: x + "_" + 'S' + str(state_index + 1), ROC_TERMS_LIST) # [u'tpr_S18', u'fpr_S18', u'precision_S18', u'recall_S18']
		columns_to_get += colnames_this_state
	return columns_to_get

def get_roc_data_one_cell_type(fn, num_chromHMM_model, num_score_bins):
	df = pd.read_csv(fn, sep = '\t', header = 0, index_col = 0)
	columns_to_get = get_columns_names_to_fetch(num_chromHMM_model) # column names that we want to get from each roc file
	try: 
		df = df[columns_to_get] # get only the columns that we are interested in
		assert df.shape[0] == num_score_bins, "Number of score bins in files: " + fn + " IS NOT EQUAL TO NUMBER OF SCORE BINS " + str(num_score_bins)
		# assert that the number of rows is the same as the number of score bins in each cell type's roc data file
	except: 
		print "Data in file: " + fn + " IS NOT IN VALID FORM"
		exit(1)
	return (df)

def calculate_summary_staistics_across_ct(cg_dir, out_dir, num_chromHMM_model, num_score_bins, ct_list):
	ct_roc_fn_list = map(lambda x: os.path.join(cg_dir, 'val_' + x, 'roc', 'tpr_fpr_all_states.txt.gz'), ct_list)
	# get files like /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood/baseline/val_E034/roc/tpr_fpr_all_states.txt.gz for all cell types
	map(helper.check_file_exist, ct_roc_fn_list) # make sure that we have all the roc data files that we need
	roc_df_list = map(lambda x: get_roc_data_one_cell_type(x, num_chromHMM_model, num_score_bins), ct_roc_fn_list)
	columns_to_get = get_columns_names_to_fetch(num_chromHMM_model) # column names that we got from each roc file
	sum_df = pd.DataFrame()
	for colname in columns_to_get: # loop through each columns that we are interested in
		this_col_data_all_ct = map(lambda x: x[colname], roc_df_list) 
		this_col_data_all_ct = np.vstack(this_col_data_all_ct) # a matrix where rows are ccell types, columns are different points on the roc/ precision-recall curves
		this_col_avg =  np.nanmean(this_col_data_all_ct, axis = 0) # a list of average roc data for all points on the roc curves
		this_col_min = np.nanmin(this_col_data_all_ct, axis = 0) # min across cell types
		this_col_max = np.nanmax(this_col_data_all_ct, axis = 0)
		this_col_std = np.nanstd(this_col_data_all_ct, axis = 0) # std deviation of this data point across the cell types
		sum_df['avg_' + colname] = this_col_avg
		sum_df['min_' + colname] = this_col_min
		sum_df['max_' + colname] = this_col_max
		sum_df['std_' + colname] = this_col_std
	output_fn = os.path.join(out_dir, 'summary_all_ct_roc.pickle')
	outF = open(output_fn, 'ab')
	pickle.dump(sum_df, output_fn) # save to a pickle file


def main():
	if len(sys.argv) != 6:
		usage()
	cg_dir = sys.argv[1]
	helper.check_dir_exist(cg_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	num_chromHMM_model = helper.get_command_line_integer(sys.argv[3])
	num_score_bins = helper.get_command_line_integer(sys.argv[4])
	cell_type_list_fn = sys.argv[5]
	ct_list = helper.get_list_from_line_seperated_file(cell_type_list_fn)
	helper.check_file_exist(cell_type_list_fn)
	print "Done getting command line arguments"
	calculate_summary_staistics_across_ct(cg_dir, out_dir, num_chromHMM_model, num_score_bins, ct_list)
	print "Done!"

def usage():
	print "python summarize_roc_data_for_cell_group.py"
	print "cg_dir: directory where the validation data for all the cell types are stored for this cell group"
	print "out_dir: where the summazied roc data for across the cell types inside this cell group is stored"
	print "num_chromHMM_model"
	print "num_score_bins"
	print "cell_type_list_fn: fn where each cell types are stored separate lines"
	exit(1)
main()