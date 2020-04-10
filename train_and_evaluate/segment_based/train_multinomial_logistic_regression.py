import pandas as pd 
import numpy as np 
from sklearn.linear_model import LogisticRegression
import os
import sys
import time
import glob
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
import multiprocessing as mp
TRAIN_MODE_LIST = ['baseline', 'multi_logistic', 'ml_chromHMM_pos']
def get_X_colnames (train_cell_types, num_chromHMM_state):
	'''
	train_cell_types = ['E034', 'E037']
	num_chromHMM_state = 2
	--> ['E034_1', 'E034_2', 'E037_1', 'E037_2']
	'''
	results = []
	for ct in train_cell_types:
		results += map(lambda x: ct + "_" + str(x + 1), range(num_chromHMM_state))
	return results

def get_state_number(state_annot):
	# 'E18' --> 18 (integers)
	return int(state_annot[1:])

def get_binary_state_assignment(state_train_segment, num_chromHMM_state): 
	'''
	a series of state assignment with indices ['E034', 'E037'] --> [1, 2]
	num_chromHMM_state = 2
	--> a series of 0/1 with indices ['E034_1', 'E034_2', 'E037_1', 'E037_2'] --> [1,0,0,1]
	'''
	results = [0] * (len(state_train_segment) * num_chromHMM_state)
	for ct_index, train_segment in enumerate(state_train_segment): 
		index_to_put_one = (int(train_segment) - 1) + ct_index * num_chromHMM_state
		results[index_to_put_one] = 1
	return pd.Series(results)

def get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_segment_fn, train_mode):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X) and response (Y). Predictor will be binarized with ct_state combinations as columns. Response will be just state labels 1 --> num_chromHMM_state. Rows of these two data frames correspond to different positions on the genome.
	train_segment_df = pd.read_csv(train_segment_fn, sep = '\t', header = 0)
	Xtrain_segment_df = train_segment_df[train_cell_types] # only get the data for the cell types that we will use to train the model 
	Y_df = train_segment_df[response_ct] # get the state train_segmentation for the response variable
	Y_df = Y_df.apply(get_state_number)
	Xtrain_segment_df = Xtrain_segment_df.applymap(get_state_number) # get the state train_segmentation from 'E18' --> 18 (integers)
	Xtrain_segment_df = Xtrain_segment_df.apply(lambda x: get_binary_state_assignment(x, num_chromHMM_state), axis = 1) # apply function row-wise, change from the state train_segmentation to binarized of state
	Xtrain_segment_df.columns = get_X_colnames(train_cell_types, num_chromHMM_state)
	return Xtrain_segment_df, Y_df # X is binarized, Y is just state label 1 --> num_chromHMM_state

def get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X). Predictor will be binarized with ct_state combinations as columns. Rows of this data frame correspond to different positions on the genome.
	segment_df = pd.read_csv(segment_fn, sep = '\t', header = 0)
	segment_df = segment_df[train_cell_types] # only get the data of the cell types that we need as predictors
	segment_df = segment_df.applymap(get_state_number)# get the state train_segmentation from 'E18' --> 18 (integers)
	segment_df = segment_df.apply(lambda x: get_binary_state_assignment(x, num_chromHMM_state), axis = 1) # apply function row-wise, change from the state train_segmentation to binarized of state
	segment_df.columns = get_X_colnames(train_cell_types, num_chromHMM_state)
	return segment_df

def train_multinomial_logistic_regression(X_df, Y_df, num_chromHMM_state):
	# give the Xtrain_segment_df and Y_df obtained from get_XY_segmentation_data --> train a logistic regression object
	regression_machine = LogisticRegression(random_state = 0, solver = 'lbfgs', multi_class = 'multinomial').fit(X_df, Y_df)
	return regression_machine 

def train_model(X_df, Y_df, num_chromHMM_state, train_mode):
	if train_mode == 'baseline':
		return 
	if train_mode == 'multi_logistic':
		return train_multinomial_logistic_regression(X_df, Y_df, num_chromHMM_state)

def predict_baseline_segmentation(predictor_df, num_chromHMM_state):
	# predictor_df : <ct>_state<state_index> --> each columns corresponds to 0/1 based on segmentation of that specific cell type
	# rows: genomic positions
	# return response_df: rows: genomic pos, columns: state<state_index> --> probability that each pos is in each state
	num_train_ct = len(predictor_df.columns) / num_chromHMM_state 
	train_df_list = []
	state_colnames = map(lambda x: "state_" + str(x+1), range(num_chromHMM_state)) # ex: [state_1, state_2, ..., state_18]
	for ct_index in range(num_train_ct): # get data for each ct that we use to train the model
		start_colname_index = ct_index * num_chromHMM_state # starting from <ct>_state1
		end_colname_index = (ct_index + 1) * num_chromHMM_state # till <ct>_state<num_chromHMM_state>
		this_ct_train_df = predictor_df[predictor_df.columns[start_colname_index:end_colname_index]] # get the data corresponding to states in this particular cell type
		this_ct_train_df.columns = state_colnames
		train_df_list.append(this_ct_train_df)
	# Now we have obtained the training data for each of the cell types, we will average out the one-hot encoding of the state assignment across cell types. That's the basis of the baseline training where basically the state assigned to each position is the state that have the highest count of cell types where this position is assigned to. 
	response_df = pd.concat(train_df_list).groupby(level = 0).mean() # get the average across all the df. So what we get is a df : rows: genomic positions, columns: states, each entry is the average of the respective cells in all the input dfa
	return response_df

def predict_regression_segmentation(predictor_df, regression_machine, num_chromHMM_state):
	response_df = regression_machine.predict_proba(predictor_df) # --> 2D: rows: positions (observations), columns: states (types) --> probability that each obs is of each type
	response_df = pd.DataFrame(response_df) # convert to a dataframe 
	# 3. Turn the results into readable format, then write to file. 
	response_df.columns = map(lambda x: "state_" + str(x + 1), range(num_chromHMM_state))
	return response_df

def predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine, train_mode):
	# based on the machine created through training, predict the segmentation corresponding to one specific window on the genome, specified in segment_fn. And print out, for each position, and for each chromHMM state, the probability that the region fall in to the state.
	# 1. Get the data of predictor cell types into the right format
	predictor_df = get_predictorX_segmentation_data(train_cell_types, num_chromHMM_state, segment_fn) # rows: predictor cell type - state combinations, columns: positions inside a  window on the genome
	# 2. Do the prediction job. Different model has different prediction functions
	if train_mode == 'baseline':
		response_df = predict_baseline_segmentation(predictor_df, num_chromHMM_state)
	elif train_mode == 'multi_logistic':
		response_df = predict_regression_segmentation(predictor_df, regression_machine, num_chromHMM_state)
	response_df.to_csv(output_fn, header = True, index = False, compression = 'gzip', sep = '\t')
	print "Done producing file: " + output_fn


def one_job_run_predict_segmentation(segment_fn_list, output_fn_list, train_cell_types, response_ct, num_chromHMM_state, regression_machine, train_mode):
	'''
	segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).  
	Each element corresponds to a region on the genome
	'''
	for (window_index, segment_fn) in enumerate(segment_fn_list):
		output_fn = output_fn_list[window_index]
		predict_segmentation_one_genomic_window(segment_fn, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine, train_mode)
	return 

def partition_file_list(file_list, num_cores):
	results = [] # list of lists of file names
	num_files_per_core = len(file_list) / num_cores
	for core_i in range(num_cores):
		if core_i < (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core : (core_i + 1) * num_files_per_core]
		elif core_i == (num_cores - 1):
			this_core_files = file_list[core_i * num_files_per_core :]
		results.append(this_core_files)
	return results

def predict_segmentation (all_ct_segment_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state, train_mode):
	# 1. Get list of segmentation files corresponding to different windows on the genome.
	segment_fn_list = glob.glob(all_ct_segment_folder + "/*.bed.gz")
	genome_pos_list = map(lambda x: (x.split('/')[-1]).split('_combined_segment.bed.gz')[0],  segment_fn_list) # from /path/to/chr9_14_combined_segment.bed.gz --> chr9_14
	output_fn_list = map(lambda x: os.path.join(predict_outDir, x + "_pred_out.txt.gz"), genome_pos_list) # get the output file names corresponding to different regions on the genome
	# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
	num_cores = 4
	partition_segment_fn_list = partition_file_list(segment_fn_list, num_cores)
	partition_output_fn_list = partition_file_list(output_fn_list, num_cores)
	processes = [mp.Process(target = one_job_run_predict_segmentation, args = (partition_segment_fn_list[i], partition_output_fn_list[i], train_cell_types, response_ct, num_chromHMM_state, regression_machine, train_mode)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print "Process " + str(i) + " is finished!"



def main():
	num_mandatory_args = 8
	if len(sys.argv) < num_mandatory_args: 
		usage()
	train_segment_fn = sys.argv[1]
	helper.check_file_exist(train_segment_fn)
	all_ct_segment_folder = sys.argv[2] # where the segmentation data of all cell types are combined, and stored in files corresponding to different regions in the genome.
	if not os.path.isdir(all_ct_segment_folder):
		print "all_ct_segment_folder IS NOT VALID: " + all_ct_segment_folder
		usage()
	predict_outDir = sys.argv[3]
	helper.make_dir(predict_outDir)
	response_ct = sys.argv[4]
	try: 
		num_chromHMM_state = int(sys.argv[5])
		assert num_chromHMM_state > 0, "num_chromHMM_state needs to be positive"
		num_train_ct = int(sys.argv[6])
		assert num_train_ct > 0, "num_train_ct needs to be positive"
	except:
		print "num_chromHMM_state or num_train_ct is not valid"
		usage()
	train_mode = sys.argv[7]
	if len(sys.argv) != (num_train_ct + num_mandatory_args):
		print "num_train_ct is different from the number of arguments passed into the program"
		usage()
	print "Done getting command line arguments"
	train_cell_types = sys.argv[num_mandatory_args:] # the rest of the arguments are the cell types that we use to train the model
	# 1. Get the data of predictors and response for training
	Xtrain_segment_df, Y_df = get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_segment_fn, train_mode)
	print "Done getting one hot data"
	print Xtrain_segment_df.head()
	print 
	print Y_df.head()
	# 2. Get the regression machine
	regression_machine = train_model(Xtrain_segment_df, Y_df, num_chromHMM_state, train_mode)
	print "Done training"
	# 3. Based on the machine just created, process training data and then predict the segmentation at each position for the response_ct
	predict_segmentation (all_ct_segment_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state, train_mode)
	print "Done predicting whole genome"
def usage():
	print "train_multinomial_logistic_regression.py "
	print "train_segment_fn: where the state assignment and of training data are stored for all cell types"
	print "all_ct_segment_folder: where segmentation data of all cell types are stored, for the entire genome, so that we can get data for prediction out."
	print "predict_output_fn: where the beta values obtained from training the data will be stored"
	print "response_ct: the cell type that we are trying to predict from the training dataset. This data is the Y value in our model training"
	print "num_chromHMM_state: Number of chromHMM states that are shared across different cell types"
	print "num_train_ct: number of cell types that we will train"
	print "train_mode: \"Normal\" \"chromHMM_pos\", \"base_line\""
	print "train_cell_types: space-separated cell types that we will use to train the model. They will be used as X values"
	exit(1)

main()