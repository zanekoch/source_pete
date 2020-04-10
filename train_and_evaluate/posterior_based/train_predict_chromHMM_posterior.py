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
def get_state_number(state_annot):
	# 'E18' --> 18 (integers)
	return int(state_annot[1:])

def get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_segment_fn):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X) and response (Y). Predictor will be binarized with ct_state combinations as columns. Response will be just state labels 1 --> num_chromHMM_state. Rows of these two data frames correspond to different positions on the genome.
	train_segment_df = pd.read_csv(train_segment_fn, sep = '\t', header = 0)
	X_colnames_to_train = []
	for ct in train_cell_types:
		this_ct_colnames = map(lambda x: ct + "_S" + str(x+1), range(num_chromHMM_state)) # ex: 'E047_S16', 'E047_S17'
		X_colnames_to_train += this_ct_colnames
	Xtrain_segment_df = train_segment_df[X_colnames_to_train] # only get the data for the cell types that we will use to train the model 
	Y_df = train_segment_df[response_ct] # get the state train_segmentation for the response variable
	Y_df = Y_df.apply(get_state_number)
	return Xtrain_segment_df, Y_df # X is binarized, Y is just state label 1 --> num_chromHMM_state

def get_predictorX_segmentation_data(train_ct_pos_folder_list, train_cell_types, num_chromHMM_state, chrom_index):
	# given the segmentation data of multiple cell types, we want to extract information from the cell types that we will use as predictor (X). Predictor will be binarized with ct_state combinations as columns. Rows of this data frame correspond to different positions on the genome.
	pos_chrom_fn_list = map(lambda x: os.path.join(train_ct_pos_folder_list[x], train_cell_types[x] + '_18_core_K27ac_chr' + str(chrom_index) + '_posterior.txt.gz'), range(len(train_ct_pos_folder_list))) # list of:  'project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/E003/E003_18_core_K27ac_chr1_posterior.txt.gz', project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/E004/E004_18_core_K27ac_chr1_posterior.txt.gz, etc. 
	posterior_df = pd.DataFrame()
	posterior_df_list = [] # each entry corresponds to a ct --> df of that ct, in this specified chromosome and contains data of posterior probabilities. 
	for ct_index, ct in enumerate(train_cell_types):
		this_ct_pos_fn = pos_chrom_fn_list[ct_index] # the fn of containing posterior probabilities of the chromosome chrom_index in cell type ct
		this_ct_pos_df = pd.read_csv(this_ct_pos_fn, skiprows = [0], header = 0, sep = '\t') # posterior prob for this ct in this chromosome. Columns will be E1 --> E<num_chromHMM_state>
		this_ct_colnames = map(lambda x: ct + "_S" + x[1:], this_ct_pos_df.columns) # convert from 'E1' --> <ct>_S1 such as E003_S1
		this_ct_pos_df.columns = this_ct_colnames
		posterior_df_list.append(this_ct_pos_df)
	posterior_df = pd.concat(posterior_df_list, axis = 1)# concat by columns : columsn are put one next to the others. Same rows across different 
	return posterior_df

def train_multinomial_logistic_regression(X_df, Y_df, num_chromHMM_state):
	# give the Xtrain_segment_df and Y_df obtained from get_XY_segmentation_data --> train a logistic regression object
	regression_machine = LogisticRegression(random_state = 0, solver = 'lbfgs', multi_class = 'multinomial').fit(X_df, Y_df)
	return regression_machine 


def predict_regression_segmentation(predictor_df, regression_machine, num_chromHMM_state):
	response_df = regression_machine.predict_proba(predictor_df) # --> 2D: rows: positions (observations), columns: states (types) --> probability that each obs is of each type
	response_df = pd.DataFrame(response_df) # convert to a dataframe 
	# 3. Turn the results into readable format, then write to file. 
	response_df.columns = map(lambda x: "state_" + str(x + 1), range(num_chromHMM_state))
	return response_df

def predict_segmentation_one_genomic_window(train_ct_pos_folder_list, chrom_index, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine):
	# we will print out predictions for the response_ct one chromosome at a time
	# train_ct_pos_folder_list --: list of paths to folders that represent the chromHMM posterior for each cell types that we use to train to predict the segmentation of response_ct
	# based on the machine created through training, predict the segmentation corresponding to one specific window on the genome, specified in segment_fn. And print out, for each position, and for each chromHMM state, the probability that the region fall in to the state.
	# 1. Get the data of predictor cell types into the right format
	predictor_df = get_predictorX_segmentation_data(train_ct_pos_folder_list, train_cell_types, num_chromHMM_state, chrom_index) # columns: predictor cell type - state combinations --> posterior, rows: positions inside a  window on the genome. The column names in predictor_df are exactly the same as those used in training the model
	# 2. Do the prediction job. Different model has different prediction functions
	response_df = predict_regression_segmentation(predictor_df, regression_machine, num_chromHMM_state)
	response_df.to_csv(output_fn, header = True, index = False, compression = 'gzip', sep = '\t')
	print "Done producing file: " + output_fn


def one_job_run_predict_segmentation(train_ct_pos_folder_list, chrom_index_list, predict_outDir, train_cell_types, response_ct, num_chromHMM_state, regression_machine):
	'''
	segment_fn_list and output_fn_list: the orders of regions in these two lists are similar (look at function predict_segmentation).  
	Each element corresponds to a region on the genome
	'''
	for chrom_index in (chrom_index_list):
		output_fn = os.path.join(predict_outDir, 'chr' + str(chrom_index) + '_pred_out.txt.gz')
		predict_segmentation_one_genomic_window(train_ct_pos_folder_list, chrom_index, output_fn, train_cell_types, response_ct, num_chromHMM_state, regression_machine)
	return 


def predict_segmentation (all_ct_posterior_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state):
	# 1. Get list of posterior folder representing different cell types
	train_ct_pos_folder_list = map(lambda x: os.path.join(all_ct_posterior_folder, x), train_cell_types)

	# 2. partition the list of file names into groups, for later putting into jobs for multiple processes
	num_cores = 2
	chrom_list_partition = helper.partition_file_list(helper.CHROMOSOME_LIST, num_cores)
	print chrom_list_partition
	processes = [mp.Process(target = one_job_run_predict_segmentation, args = (train_ct_pos_folder_list, chrom_list_partition[i], predict_outDir, train_cell_types, response_ct, num_chromHMM_state, regression_machine)) for i in range(num_cores)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print "Process " + str(i) + " is finished!"



def main():
	num_mandatory_args = 7
	if len(sys.argv) < num_mandatory_args: 
		usage()
	train_segment_fn = sys.argv[1]
	helper.check_file_exist(train_segment_fn)
	all_ct_posterior_folder = sys.argv[2] # where the segmentation data of all cell types are combined, and stored in files corresponding to different regions in the genome.
	helper.check_dir_exist(all_ct_posterior_folder)
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
	if len(sys.argv) != (num_train_ct + num_mandatory_args):
		print "num_train_ct is different from the number of arguments passed into the program"
		usage()
	print "Done getting command line arguments"
	train_cell_types = sys.argv[num_mandatory_args:] # the rest of the arguments are the cell types that we use to train the model
	# 1. Get the data of predictors and response for training
	Xtrain_segment_df, Y_df = get_XY_segmentation_data (train_cell_types, response_ct, num_chromHMM_state, train_segment_fn) # Xtrain_segment_df: example colnames: 'E047_S16', 'E047_S17' --> posterior probabilities of each of the state in each cell type that are used to train
	# Y_df --> example colnames 'E047' --> state numbers 1 --> 18 of each position used to train data for the response cell type
	print "Done getting one hot data"
	print Xtrain_segment_df.head()
	print 
	print Y_df.head()
	# 2. Get the regression machine
	regression_machine = train_multinomial_logistic_regression(Xtrain_segment_df, Y_df, num_chromHMM_state)
	print "Done training"
	# 3. Based on the machine just created, process training data and then predict the segmentation at each position for the response_ct
	predict_segmentation (all_ct_posterior_folder, regression_machine, predict_outDir, train_cell_types, response_ct, num_chromHMM_state)
	print "Done predicting whole genome"

def usage():
	print "python train_predict_chromHMM_posterior.py "
	print "train_segment_fn: where the state assignment and of training data are stored for all cell types"
	print "all_ct_posterior_folder: where chromHMM posterior data of all cell types are stored, for the entire genome, so that we can get data for prediction out. Note, the structure of the data is such that each subfolder is a cell type. Within each subfolder, each file represent the posterior probabilities in one chromosome"
	print "predict_outDir: where the probabilities of each position being in a particular state will be stored for each cell type"
	print "response_ct: the cell type that we are trying to predict from the training dataset. This data is the Y value in our model training"
	print "num_chromHMM_state: Number of chromHMM states that are shared across different cell types"
	print "num_train_ct: number of cell types that we will train"
	print "train_cell_types: space-separated cell types that we will use to train the model. They will be used as X values"
	exit(1)

main()