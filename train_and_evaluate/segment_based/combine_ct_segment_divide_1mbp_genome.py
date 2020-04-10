import pandas as pd 
import numpy as np 
import os
import sys
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
import glob
import multiprocessing as mp
CHROMOSOME_LIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'] # we exclude chromosome Y because it's not available in all cell types
NUM_BP_PER_BIN = 200
NUM_BIN_PER_WINDOW = 50000 # each window is 10 000 000 base pairs
NUM_BP_PER_WINDOW = NUM_BP_PER_BIN * NUM_BIN_PER_WINDOW 

def open_segment_df(segment_fn):
	result_df = pd.read_csv(segment_fn, sep = '\t', header = None)
	result_df.columns = ['chr', 'start_bp', 'end_bp', 'state']
	return result_df


def get_ct_segment_df (org_ct_segment_folder):
	ct_dir_list = glob.glob(org_ct_segment_folder + "/E*/")
	ct_list = map(lambda x: x.split('/')[-2], ct_dir_list) # /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/E003/ --> E003
	ct_segment_fn_list = map(lambda x: os.path.join(x[0], x[1] + '_18_core_K27ac_segments.bed.gz'), zip(ct_dir_list, ct_list)) # /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/E003/E003_18_core_K27ac_segments.bed.gz
	ct_df_list = map(open_segment_df, ct_segment_fn_list) # a list of df representing the segmentation data for each of the cell type. Orders are preserved across all the list of cell types' entities.
	return ct_list, ct_df_list

def get_chromosome_length(org_ct_segment_folder, ct_list):
	first_ct_chrom_length_fn = os.path.join(org_ct_segment_folder, ct_list[0], ct_list[0] + '_chrom_length.bed')
	chrom_len_df = pd.read_csv(first_ct_chrom_length_fn, sep = '\t', header = None)
	chrom_len_df = chrom_len_df[[0,2]] # get the first and last column --> chromosome and length
	chrom_len_df.columns = ['chr', 'length']
	chrom_len_dict = pd.Series(chrom_len_df.length.values, index = chrom_len_df.chr).to_dict() # keys: chromosome, values: length of the chromosome
	return chrom_len_dict

def get_200bp_one_segment_df(ct, chromosome, this_chrom_df, start_bp_this_window, end_bp_this_window):
	'''
	this_chrom_df : segmentation data for this chromosome, sorted from coodinate 0 --> chromosome length
	'''
	start_window_row = this_chrom_df[(this_chrom_df['start_bp'] <= start_bp_this_window) & (this_chrom_df['end_bp'] >= start_bp_this_window)] # the row that contains the start of this window
	try:
		row_index_start_window = start_window_row.index[0] # find the index of line that marks the start of this 1mbp window
	except:
		print "start_bp_this_window: " + str(start_bp_this_window)
		print "end_bp_this_window: " + str(end_bp_this_window)
		print "ct: " + str(ct)
		print "chromosome: " + str(chromosome)
		num_segments = (end_bp_this_window - start_bp_this_window) / NUM_BP_PER_BIN
		return pd.Series(['E-1'] * num_segments)
	last_row_index = this_chrom_df.index[-1]
	results = []
	for row_index in range(row_index_start_window, last_row_index + 1):
		this_state = (this_chrom_df.loc[row_index])['state']
		if (row_index == row_index_start_window): # first segment of this window
			num_segments = ((this_chrom_df.loc[row_index])['end_bp'] - start_bp_this_window) / NUM_BP_PER_BIN
			if num_segments > NUM_BIN_PER_WINDOW: # this window only needs data from one row in the original segmentation, because this one row spans more than NUM_BP_PER_WINDOW bp
				results += [this_state] * NUM_BIN_PER_WINDOW
				break 
			else: # this row of the starting bp does not stretch beyond the NUM_BP_PER_WINDOW bp in one window
				results += [this_state] * num_segments
				continue
		if (this_chrom_df.loc[row_index])['end_bp']  >= end_bp_this_window: # last segment of this window
			num_segments = (end_bp_this_window - (this_chrom_df.loc[row_index])['start_bp']) / NUM_BP_PER_BIN
			results += [this_state] * num_segments
			break 
		num_segments = ((this_chrom_df.loc[row_index])['end_bp'] - (this_chrom_df.loc[row_index])['start_bp']) / NUM_BP_PER_BIN
		results += [this_state] * num_segments
	return pd.Series(results) 

def combine_segment_one_1Mb_window(ct_list, ct_df_list, chromosome, this_chromosome_length, start_bp_this_window, output_fn):
	'''
	chromosome: chr1, chr2, chr3, chr4, etc. 
	This function is tested and veried on two cell types E003 and E004
	'''
	end_bp_this_window = min(this_chromosome_length, start_bp_this_window + NUM_BP_PER_WINDOW)
	num_segment_this_window = (end_bp_this_window - start_bp_this_window) / NUM_BP_PER_BIN # number of 200bp segments in this window (5000 or something smaller)
	result_df = pd.DataFrame()
	# get the data of genomic positions first
	result_df['chrom'] = pd.Series([chromosome] * num_segment_this_window)
	result_df['start_bp_this_window'] = start_bp_this_window + NUM_BP_PER_BIN * np.array(range(num_segment_this_window))
	result_df['end_bp_this_window'] = result_df['start_bp_this_window'] + NUM_BP_PER_BIN
	for ct_index, ct in enumerate(ct_list):
		this_ct_segment_df = ct_df_list[ct_index]
		this_chrom_df = this_ct_segment_df[this_ct_segment_df['chr'] == (str(chromosome))]
		result_df[ct] = get_200bp_one_segment_df(ct, chromosome, this_chrom_df, start_bp_this_window, end_bp_this_window)
	result_df.to_csv(output_fn, sep = '\t', header = True, index = False, compression = 'gzip') 
	# --> chromosome, start_bp, end_bp, ct1, ct2, etc. --> segmentations data, broken down to 200bp segments in all the cell types that we investigate, across the window that we care about. 

def partition_genome_window(this_chromosome_length, num_cores):
	# this function is tested and verified
	num_window_this_chrom = int(np.ceil(float(this_chromosome_length) / NUM_BP_PER_WINDOW))
	num_window_per_core = int(num_window_this_chrom / num_cores)
	results = []
	for i in range(num_cores):
		start_window_index = int(i * num_window_per_core)
		end_window_index = int((i + 1) * num_window_per_core)
		if i == (num_cores - 1): 
			end_window_index = num_window_this_chrom
		this_core_start_window_bp_list = NUM_BP_PER_WINDOW * np.array(range(start_window_index, end_window_index))
		results.append(this_core_start_window_bp_list)
	return results

def call_combine_one_1Mb_window_from_one_core(ct_list, ct_df_list, output_folder, chromosome, this_chromosome_length, start_window_bp_list):
	output_fn_list = map(lambda x: os.path.join(output_folder, chromosome + "_" + str(x / NUM_BP_PER_WINDOW) + "_combined_segment.bed.gz"), start_window_bp_list)
	for w_i, w_start_bp in enumerate(start_window_bp_list):
		this_window_output_fn = output_fn_list[w_i]
		 (ct_list, ct_df_list, chromosome, this_chromosome_length, w_start_bp, this_window_output_fn)
		#print "Done creating file : " + this_window_output_fn


def combine_segment_in_parallel(ct_list, ct_df_list, chrom_len_dict, output_folder):
	num_cores = 4 # run 4 jobs in parallel
	for chromosome in chrom_len_dict: # loop through each chromosome at a time
		print "Running jobs for chromosome " + chromosome
		this_chromosome_length = chrom_len_dict[chromosome]
		partition_start_window_bp_list = partition_genome_window(this_chromosome_length, num_cores)
		processes = [mp.Process(target = call_combine_one_1Mb_window_from_one_core, args = (ct_list, ct_df_list, output_folder, chromosome, this_chromosome_length, partition_start_window_bp_list[i])) for i in range(num_cores)]
		print "Starting processes"
		for p in processes: 
			p.start()
		for i, p in enumerate(processes):
			p.join()
			print "Process " + str(i) + " is finished!"

def main():
	if len(sys.argv) != 3:
		usage()
	org_ct_segment_folder = sys.argv[1]
	if not os.path.isdir(org_ct_segment_folder): 
		print "org_ct_segment_folder IS NOT VALID: " + org_ct_segment_folder
		usage()
	output_folder = sys.argv[2]
	helper.make_dir(output_folder)
	print "Done getting command line arguments"
	ct_list, ct_df_list = get_ct_segment_df (org_ct_segment_folder)
	print "Done getting segment_df for all cell types"
	chrom_len_dict = get_chromosome_length(org_ct_segment_folder, ct_list)
	print "Done getting chromosome length"
	combine_segment_in_parallel(ct_list, ct_df_list, chrom_len_dict, output_folder)
	print ""
	print ""
	print "Done!"

def usage(): 
	print "python combine_ct_segment_divide_1percent_genome.py"
	print "org_ct_segment_folder: Example: /u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded"
	print "output_folder : where the files of different genomic 1mbp locations will be saved"
	exit(1)

main()