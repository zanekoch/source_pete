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
	ct_list =  list(pd.read_csv(all_ct_list_fn, sep = '\n', header = None)[0]) # -->  a list with each entry being the cell type in this cell group
	# get the index of the validation ct and then remove it from the list, so that we can focus this training process on ct other than the validation ct
	val_ct_i = ct_list.index(validate_ct)
	assert val_ct_i != -1, "the validation ct is not present in the list of all cell type of the group that we are trying to train on"
	ct_list = ct_list[:val_ct_i] + ct_list[(val_ct_i + 1):] # skip the validation ct
	return ct_list

def get_genomic_positions_list(all_ct_segment_folder):
	gen_pos_segment_fn_list = glob.glob(all_ct_segment_folder + '/chr*')
	gen_pos_list = map(lambda x: (x.split('/')[-1]).split('_combined_segment.bed.gz')[0], gen_pos_segment_fn_list) # get the list of all genomic positions available: [chr9_11, chr9_12, etc.]
	return gen_pos_list

def put_one_result_file_to_df(fn):
	return pd.read_csv(fn, header = 0, sep = '\t')

def average_multiple_result_files(result_fn_list, output_fn):
	result_df_list = map(put_one_result_file_to_df, result_fn_list) # read all the files data and put them into  a data frame
	avg_df = pd.concat(result_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	avg_df.to_csv(output_fn, compression = 'gzip', header = True, index = True, sep = '\t') # save and compression to file
	return

def averaging_predictions_to_validate_one_ct(validate_ct_dir, validate_ct, num_pred_ct, gen_pos_list):
	# num_pred_ct: number of ct whose predictions we use to average out and get the predictions for the validate cell type
	# validate_ct_dir = the directory where the data that are associated with the ct used for validation  cell type
	pred_dir_list = glob.glob(validate_ct_dir + "/pred_*")
	pred_ct_list = map(lambda x: (x.split('/')[-1]).split('_')[-1], pred_dir_list) # path/to/pred_E034 --> E034
	assert len(pred_dir_list) == num_pred_ct, 'Number of pred_dir_list is not the same as number of specificed ct used to predict the model'
	# get the folder where the results of averaging across different prediction cell types will be stored
	validate_outDir = os.path.join(validate_ct_dir, 'average_predictions')
	helper.make_dir(validate_outDir)
	for gene_window in gen_pos_list: # loop through each genomic window and then get the avrage result across different predictions for all positions in this window
		this_window_output_fn = os.path.join(validate_outDir, gene_window + "_avg_pred.txt.gz")
		this_window_pred_fn_list = map(lambda x: os.path.join(x, gene_window + "_pred_out.txt.gz"), pred_dir_list)
		# calculate the average prediction results across different prediction cell types for this window, and save the results
		average_multiple_result_files(this_window_pred_fn_list, this_window_output_fn)
	return 



def call_cross_validation_functions(validate_ct, ct_list, outDir, train_sampled_data_fn, all_ct_segment_folder, num_chromHMM_state, gen_pos_list, train_mode):
	code_file_fn = '/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/train_multinomial_logistic_regression.py'
	val_outDir = os.path.join(outDir, 'val_' + validate_ct)
	print "Running validation of ct: " + validate_ct
	for train_ct_i, train_ct in enumerate(ct_list): # this ct will be used be the response variable for training the data.
		predictor_ct_list = ct_list[:train_ct_i] + ct_list[(train_ct_i + 1):] # leave out the response ct. All the remaining ones will be used as predictor and will be passed into the program
		num_predictor_ct = len(predictor_ct_list)			
		this_predict_outDir = os.path.join(val_outDir, 'pred_' + train_ct)
		helper.make_dir(this_predict_outDir)
		command = ['python', code_file_fn, train_sampled_data_fn, all_ct_segment_folder, this_predict_outDir, train_ct, str(num_chromHMM_state), str(num_predictor_ct), train_mode] + predictor_ct_list
		print "Within, running predicting cell type: " + train_ct
		call(command)
	print "Averaging results from different predictions for this validation"
	averaging_predictions_to_validate_one_ct(validate_ct_dir = val_outDir, validate_ct = validate_ct, num_pred_ct = len(ct_list), gen_pos_list = gen_pos_list)
	print ""
	print ""


def main():
	if len(sys.argv) != 8:
		usage()
	train_sampled_data_fn = sys.argv[1]
	helper.check_file_exist(train_sampled_data_fn)
	outDir = sys.argv[2]
	helper.make_dir(outDir)
	all_ct_segment_folder = sys.argv[3]
	helper.check_dir_exist(all_ct_segment_folder)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[4]) 
	validate_ct = sys.argv[5]
	train_mode = sys.argv[6]
	all_ct_list_fn = sys.argv[7]
	print "Done getting command line arguments"
	# get the list of all genomic positions used to segment the genome for our model training (we exclude chromosome Y in all analysis)
	gen_pos_list = get_genomic_positions_list(all_ct_segment_folder)
	# get all cell types
	ct_list = get_all_train_ct_list(all_ct_list_fn, validate_ct)
	# call all cell types
	call_cross_validation_functions(validate_ct, ct_list, outDir, train_sampled_data_fn, all_ct_segment_folder, num_chromHMM_state, gen_pos_list, train_mode)


def usage():
	print "python run_cross_validation_pipeline.py"
	print "train_sampled_data_fn: where the segmentation training data of 10% genome are stored, and also information about all the cell types that we are running the pipeline on"
	print "outDir: Where all the output data will be stored in the most structured form"
	print "all_ct_segment_fn: where whole-genome segmentation data for all 127 cell types are stored"
	print "num_chromHMM_state: Num chromHMM states that we will predict the genomic positions upon"
	print "validate_ct: ct that we use to validate the training process. We will leave out this ct when we try to train the model useing the rest of the cell types as predictors and response variables"
	print "train_mode: ['baseline', 'multi_logistic', 'ml_chromHMM_pos']"
	print "all_ct_list_fn: each line of this file is the name of a cell type that is associated with the cell group that we are running the pipeline on. An example is /u/home/h/havu73/project-ernst/diff_pete/roadmap/blood.list"
	exit(1)

main()