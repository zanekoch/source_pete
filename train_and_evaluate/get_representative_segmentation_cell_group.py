import pandas as pd 
import numpy as np 
import sys
import os
import glob
import helper
import multiprocessing as mp
NUM_CORES = 4
def read_avg_state_assignment(fn):
	df = pd.read_csv(fn, sep = '\t', header = 0, index_col = 0)
	return df

def calculate_avg_segmentation_one_genome_window(genome_window, ct_dir_list, avg_segmentaiton_dir):
	ct_segment_fn_list = map(lambda x: os.path.join(x, genome_window + "_avg_pred.txt.gz"), ct_dir_list) # get the list of files that represent results from this region in each of the cell type that we validated
	ct_segment_df_list = map(read_avg_state_assignment, ct_segment_fn_list) # read all the data from each of the cell type so that we can do the average of thses
	avg_df = pd.concat(ct_segment_df_list).groupby(level = 0).mean() # get the average of segmentation df across all cell types
	save_fn = os.path.join(avg_segmentaiton_dir, genome_window + "_avg_all_ct.txt.gz")
	avg_df.to_csv(save_fn, compression = 'gzip', header = True, index = True, sep = '\t')




def read_state_annot_fn(state_annotation_fn):
	df = pd.read_csv(state_annotation_fn, header = 0, sep = '\t')
	df['state'] = map(lambda x: 'E' + str(x), df['state'])
	return df

def get_rgb_format_right(rgb):
	# convert from (255, 245, 238) to 255,245,238
	# numbers = (rgb[1:-1]).split(',') # get rid of () 
	numbers = (rgb[:]).split(',')
	numbers = map(lambda x: str(int(x)), numbers)
	return ",".join(numbers)

def get_standard_ucsc_format_bed_file(segment_df, state_annot_df, output_fn, igv_track_name, genome_pos):
	segment_df = pd.merge(segment_df, state_annot_df
		, how = 'left', left_on = 'state', right_on = 'state', left_index = False, right_index = False)
	segment_df = segment_df [['chrom', 'start_bp', 'end_bp', 'state', 'itemRgb']]
	(nrow, ncol) = segment_df.shape
	segment_df['itemRgb'] = (segment_df['itemRgb']).apply(get_rgb_format_right)
	segment_df['score'] = ['1'] * nrow
	segment_df ['strand'] = ['.'] * nrow
	segment_df['thickStart'] = segment_df['start_bp']
	segment_df['thickEnd'] = segment_df['start_bp']
	segment_df = segment_df.rename(columns = {'start_bp' : 'chromStart', 'end_bp' : 'chromEnd', 'state' : 'name'})
	segment_df = segment_df[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']]
	outF = open(output_fn, 'a')
	header_comment = "track name=\"" + igv_track_name + "_" + genome_pos + "\" description=\"\" visibility=1 itemRgb=\"On\"\n"
	outF.write(header_comment) # write the comment first so that genome browser can read the file
	segment_df.to_csv(outF, sep = '\t', header = False, index = False, compression = 'gzip')
	outF.close()


def get_igv_format_segmentation(state_segment_fn, genome_pos, state_annot_df, out_dir, igv_track_name):
	print state_segment_fn
	segment_df = pd.read_csv(state_segment_fn, header = 0, index_col = 0, sep = '\t')
	max_state = segment_df.idxmax(axis = 1) # get the column max, but return the column names
	max_state = map(lambda x: 'E' + str(x.split('_')[1]), max_state) # chagne from state_18 to E18 for all elements in the max_state series
	chr_data = genome_pos.split('_')[0]
	offset_bp = int(genome_pos.split('_')[1])
	offset_bp = offset_bp * helper.NUM_BP_PER_WINDOW # each window contains data in terms of NUM_BP_PER_WINDOW
	result_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp', 'state'])
	current_start_bin_index = 0
	current_end_bin_index = 1 
	current_state = max_state[current_start_bin_index]
	for index in range(1, len(max_state)): # skip the first one because it's already reported into the current data
		if (current_state == max_state[index]):
			current_end_bin_index = index + 1 # add one more bin to the stread of bins that are segmented here
		else: # change of state --> report into the result_df
			start_bp_report = offset_bp +  current_start_bin_index * helper.NUM_BP_PER_BIN
			end_bp_report = offset_bp + current_end_bin_index * helper.NUM_BP_PER_BIN
			add_row = [chr_data, start_bp_report, end_bp_report, current_state]
			result_df.loc[result_df.shape[0]] = add_row
			current_state = max_state[index]
			current_start_bin_index = index
			current_end_bin_index = index + 1
	start_bp_report = offset_bp + current_start_bin_index * helper.NUM_BP_PER_BIN
	end_bp_report = offset_bp + current_end_bin_index * helper.NUM_BP_PER_BIN
	add_row = [chr_data, start_bp_report, end_bp_report, current_state]
	result_df.loc[result_df.shape[0]] = add_row
	# now we are done producing
	output_fn = os.path.join(out_dir, genome_pos + "_segment_igv_format.bed.gz")
	get_standard_ucsc_format_bed_file(result_df, state_annot_df, output_fn, igv_track_name, genome_pos)

def one_job_run_calculate_avg_state_assign_matrix(genome_pos_list_one_core, ct_dir_list, avg_segmentaiton_dir):
	map(lambda x: calculate_avg_segmentation_one_genome_window(x, ct_dir_list, avg_segmentaiton_dir), genome_pos_list_one_core) # create average state semgentation matrix for each genomic window

def get_average_state_assign_matrix(cg_dir, ct_list, num_chromHMM_state, out_dir):
	ct_dir_list = map(lambda x: os.path.join(cg_dir, "val_" + x, "average_predictions"), ct_list)  # list of path to the average prediction folder for each ct of this cg
	genome_pos_fn_list = glob.glob(ct_dir_list[0] + "/*_avg_pred.txt.gz")
	genome_pos_list = map(lambda x: (x.split('/')[-1]).split('_avg_pred.txt.gz')[0], genome_pos_fn_list) # genome_pos would be 'chrX_9', 'chrX_8', etc.
	partition_genome_pos_list = helper.partition_file_list(genome_pos_list, NUM_CORES) # [core_index] --> list of genomic position whose data will be processed in this process
	avg_segmentaiton_dir = os.path.join(out_dir, "avg_segmentation")
	helper.make_dir(avg_segmentaiton_dir)
	processes = [mp.Process(target = one_job_run_calculate_avg_state_assign_matrix, args = (partition_genome_pos_list[i], ct_dir_list, avg_segmentaiton_dir)) for i in range(NUM_CORES)]
	for p in processes:
		p.start()
	for i, p in enumerate(processes):
		p.join()
		print "Process " + str(i) + " is finished!"
	print "Done creating average state segmentation matrices over all cell types"

def create_igv_format_bed(out_dir, state_annot_df, draw_genome_pos_list, igv_track_name):
	draw_bed_dir = os.path.join(out_dir, 'igv_bed')
	helper.make_dir(draw_bed_dir)
	avg_segmentaiton_dir = os.path.join(out_dir, 'avg_segmentation')
	for genome_pos in draw_genome_pos_list:
		state_segment_fn = os.path.join(avg_segmentaiton_dir, genome_pos + "_avg_all_ct.txt.gz")
		get_igv_format_segmentation(state_segment_fn, genome_pos, state_annot_df, draw_bed_dir, igv_track_name)
		print "Done getting the igv bed file for " + genome_pos


def main():
	if len(sys.argv) != 7:
		usage()
	cg_dir = sys.argv[1]
	helper.check_dir_exist(cg_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	state_annotation_fn = sys.argv[3]
	helper.check_file_exist(state_annotation_fn)
	state_annot_df = read_state_annot_fn(state_annotation_fn)
	ct_list_fn = sys.argv[4]
	helper.check_file_exist(ct_list_fn)
	ct_list = helper.get_list_from_line_seperated_file(ct_list_fn)
	num_chromHMM_state = helper.get_command_line_integer(sys.argv[5])
	igv_track_name = sys.argv[6]
	print "Done getting command line arguments"
	get_average_state_assign_matrix(cg_dir, ct_list, num_chromHMM_state, out_dir)
	print "Done getting the representative state semgentation for the cellg group"
	draw_genome_pos_list = ['chr5_15']
	# create_igv_format_bed(out_dir, state_annot_df, draw_genome_pos_list, igv_track_name)
	print "Done!"



def usage():
	print "state_annotation_fn: fn of the annotations of the states such as colors and stuff like that"
	print "ct_list_fn: fn of list of ct that are in this cell group"
	print "num_chromHMM_state"
	print "igv_track_name: the name that we will put to the tracks that will later be produced for visualization on igv"
	exit(1)

main()