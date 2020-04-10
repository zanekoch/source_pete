import pandas as pd 
import numpy as np 
import sys
import os
import helper 

def get_cell_types_of_interest(ct_fn):
	ctF = open(ct_fn, 'r')
	ct_list = []
	for line in ctF: 
		ct = line.strip()
		ct_list.append(ct)
	ctF.close()
	return ct_list

def sample_genome_positions(cell_type_folder, ct_list, output_fn):
	NUM_BP_PER_BIN = 200
	first_ct = ct_list[0]
	chrom_length_fn = os.path.join(cell_type_folder, first_ct, first_ct + "_chrom_length.bed") # we pick only the first cell type in the list because I already confirmed that the length of all chromsomes in all the cell types' segmentation are similar
	ct_df = pd.read_table(chrom_length_fn, sep = '\t', header = None)
	ct_df = ct_df[[0,2]] # get rid of the second column because they are all 1's
	ct_df.columns = ['chrom', 'length'] 
	ct_df.drop(len(helper.CHROMOSOME_LIST), axis = 0) # drop the last row because it's chromosome Y
	# ct_df has the following columns: chrom, length, end_bp, start_bp, start_bin, end_bin --> data about the chromosome
	sample_df = pd.DataFrame(columns = ['chrom', 'start_bp', 'end_bp']) # all bp index are zero-based --> a dataframe of all the positions that we will sample
	for chrom_index in helper.CHROMOSOME_LIST:
		this_chrom_length = (ct_df[ct_df['chrom'] == 'chr' + chrom_index])['length']
		num_bin_this_chrom = this_chrom_length / NUM_BP_PER_BIN 
		sample_bins_this_chrom = np.random.choice(range(num_bin_this_chrom), size = int(num_bin_this_chrom / 10), replace = False)
		sample_bins_this_chrom.sort()
		sample_this_chrom_df = pd.DataFrame() 
		sample_this_chrom_df['chrom'] = ['chr' + chrom_index] * len(sample_bins_this_chrom)
		sample_this_chrom_df['start_bp'] = sample_bins_this_chrom * NUM_BP_PER_BIN
		sample_this_chrom_df['end_bp'] = sample_this_chrom_df['start_bp'] + NUM_BP_PER_BIN
		sample_df = sample_df.append(sample_this_chrom_df)
		print "Done with chromosome: " + chrom_index
	# save to file 
	sample_df.to_csv(output_fn, sep = '\t', header = None, index = False)
	return sample_df


def main():
	if len(sys.argv) != 4:
		usage()
	cell_type_folder=sys.argv[1]
	if not os.path.isdir(cell_type_folder):
		print "cell_type_folder DOES NOT EXIST"
		usage()
	ct_fn=sys.argv[2]
	helper.check_file_exist(ct_fn)
	ct_list = get_cell_types_of_interest(ct_fn) # list of cell types of interests example: ['E003', 'E004']
	output_fn = sys.argv[3]
	helper.create_folder_for_file(ct_fn)
	print "Done getting command line arguments"
	# select regions on the genome that we will sample from
	genome_sample_df = sample_genome_positions(cell_type_folder, ct_list, output_fn) # --> a dataframe of 3 columns: "chromosome", "start_bp", 'end_bp'

	
def usage():
	print "python sample_genome.py"
	print "cell_type_folder: Folder containing all the information about the cell types that we are interested"
	print "ct_fn: file where each line is the name of a cell type that we want to sample from"
	print "output_fn: where the data of sampled regions and state segementation will be stored for all the cell types that we chose"
	print "The result should give us around 1518140 bins"
	exit(1)
main()