'''
This function will only train the cell types for one validation cell type. We need to call a bash script to create jobs for do multiple training jobs for multiple validation cell types
'''

import os
import sys
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
from subprocess import call # to call other python script like from command line arguments
import pandas as pd
import glob
import numpy as np

def get_all_train_ct_list(all_ct_list_fn, validate_ct):
	ct_list =  list(pd.read_csv(all_ct_list_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group. Note the [0] is necessary for us to get the first column
	# get the index of the validation ct and then remove it from the list, so that we can focus this training process on ct other than the validation ct
	val_ct_i = ct_list.index(validate_ct)
	assert val_ct_i != -1, "the validation ct is not present in the list of all cell type of the group that we are trying to train on"
	ct_list = ct_list[:val_ct_i] + ct_list[(val_ct_i + 1):] # skip the validation ct
	return ct_list

def put_one_result_file_to_df(fn):
	return pd.read_csv(fn, header = 0, sep = '\t')

def average_multiple_result_files(result_fn_list, chrom_index, validate_outDir):
	result_df_list = map(put_one_result_file_to_df, result_fn_list) # read all the files data and put them into  a data frame
	avg_df = pd.concat(result_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	(num_bins_this_chrom, num_chromHMM_state) = avg_df.shape # this works based on the assumption that the length of chromosomes in all the cell types are similar. This has been checked using code inthe download_data folder. See the code and run for yourself if you want to confirm that fact. each bin is 200bp long
	num_window_this_chrom = int(np.ceil(float(num_bins_this_chrom) / float(helper.NUM_BIN_PER_WINDOW))) # number of windows that we will divied this chromosome into
	for window_index in range(num_window_this_chrom - 1): 
		this_window_output_fn = os.path.join(validate_outDir, "chr" + str(chrom_index) + "_" + str(window_index) + "_avg_pred.txt.gz")
		this_window_start_bin_index = window_index  * helper.NUM_BIN_PER_WINDOW # each bin is a row on the avg_df. We want to cut a window of NUM_BIN_PER_WINDOW bins for each window
		this_window_end_bin_index = (window_index + 1) * helper.NUM_BIN_PER_WINDOW - 1
		this_window_df = avg_df.loc[this_window_start_bin_index:this_window_end_bin_index] # only get the bins in the range of indices that we want. Note that the indexing system in pandas is different from normal python. When I specifiy [0:2] it gives me 0,1,2 indices, not just 0,1 as usual
		this_window_df.to_csv(this_window_output_fn, sep = '\t', header = True, compression = 'gzip')
	# now do the last window which is probably not a full window
	last_window_index = num_window_this_chrom - 1
	last_window_output_fn = os.path.join(validate_outDir, "chr" + str(chrom_index) + "_" + str(last_window_index) + "_avg_pred.txt.gz")
	last_window_start_bin_index = last_window_index * helper.NUM_BIN_PER_WINDOW
	last_window_df = avg_df.loc[last_window_start_bin_index:]
	if last_window_df.shape[0] > 0: # if there are rows of data in this window
		last_window_df.to_csv(last_window_output_fn, sep = '\t', header = True, compression = 'gzip') # save to file when there is data. which there will definitely will
	return

def averaging_predictions_to_validate_one_ct(validate_ct_dir, validate_ct, num_pred_ct):
	# num_pred_ct: number of ct whose predictions we use to average out and get the predictions for the validate cell type
	# validate_ct_dir = the directory where the data that are associated with the ct used for validation  cell type
	pred_dir_list = glob.glob(validate_ct_dir + "/pred_*")
	pred_ct_list = map(lambda x: (x.split('/')[-1]).split('_')[-1], pred_dir_list) # path/to/pred_E034 --> E034
	assert len(pred_dir_list) == num_pred_ct, 'Number of pred_dir_list is not the same as number of specificed ct used to predict the model'
	# get the folder where the results of averaging across different prediction cell types will be stored
	validate_outDir = os.path.join(validate_ct_dir, 'average_predictions')
	helper.make_dir(validate_outDir)
	for chrom_index in helper.CHROMOSOME_LIST: # loop through each genomic window and then get the avrage result across different predictions for all positions in this window
		this_window_pred_fn_list = map(lambda x: os.path.join(x, 'chr' + str(chrom_index) +  "_pred_out.txt.gz"), pred_dir_list)
		# calculate the average prediction results across different prediction cell types for this window, and save the results
		average_multiple_result_files(this_window_pred_fn_list, chrom_index, validate_outDir) # this function will take the average of predictions of multiple cell types, for each genomic bin(200bp), and then divide the averaged data for each chromosome into 
		print "Done averaging results for chromosome: " + str(chrom_index)
	return 



def call_cross_validation_functions(validate_ct, ct_list, outDir, train_sampled_data_fn, all_ct_posterior_folder, num_chromHMM_state):
	code_file_fn = '/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/posterior_based/train_predict_chromHMM_posterior.py'
	val_outDir = os.path.join(outDir, 'val_' + validate_ct)
	print "Running validation of ct: " + validate_ct
	for train_ct_i, train_ct in enumerate(ct_list): # this ct will be used be the response variable for training the data.
		predictor_ct_list = ct_list[:train_ct_i] + ct_list[(train_ct_i + 1):] # leave out the response ct. All the remaining ones will be used as predictor and will be passed into the program
		num_predictor_ct = len(predictor_ct_list)			
		this_predict_outDir = os.path.join(val_outDir, 'pred_' + train_ct)
		helper.make_dir(this_predict_outDir)
		command = ['python', code_file_fn, train_sampled_data_fn, all_ct_posterior_folder, this_predict_outDir, train_ct, str(num_chromHMM_state), str(num_predictor_ct)] + predictor_ct_list
		print "Within, running predicting cell type: " + train_ct
		# call(command)
	print "Averaging results from different predictions for this validation"
	averaging_predictions_to_validate_one_ct(validate_ct_dir = val_outDir, validate_ct = validate_ct, num_pred_ct = len(ct_list))
	print ""
	print ""


def main():
	if len(sys.argv) != 7:
		usage()
	train_sampled_data_fn = sys.argv[1]
	helper.check_file_exist(train_sampled_data_fn)
	outDir = sys.argv[2]
	helper.make_dir(outDir)
	all_ct_posterior_folder = sys.argv[3]
	helper.check_dir_exist(all_ct_posterior_folder)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[4]) 
	validate_ct = sys.argv[5]
	all_ct_list_fn = sys.argv[6]
	print "Done getting command line arguments"
	# get all cell types
	ct_list = get_all_train_ct_list(all_ct_list_fn, validate_ct)
	print ct_list
	# call all cell types
	call_cross_validation_functions(validate_ct, ct_list, outDir, train_sampled_data_fn, all_ct_posterior_folder, num_chromHMM_state)


def usage():
	print "python run_cross_validation_pipeline.py"
	print "train_sampled_data_fn: where the segmentation training data of 10% genome are stored, and also information about all the cell types that we are running the pipeline on"
	print "outDir: Where all the output data will be stored in the most structured form"
	print "all_ct_segment_fn: where whole-genome segmentation data for all 127 cell types are stored"
	print "num_chromHMM_state: Num chromHMM states that we will predict the genomic positions upon"
	print "validate_ct: ct that we use to validate the training process. We will leave out this ct when we try to train the model useing the rest of the cell types as predictors and response variables"
	print "all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list"
	exit(1)

main()