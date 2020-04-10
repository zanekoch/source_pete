import pandas as pd
import numpy as np 
import sys
import os
import glob
sys.path.append("/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/")
import helper
def divide_one_chromosome_into_window(chrom_dir, chrom_index, out_dir):
	chrom_fn = os.path.join(0)
def main():
	if len(sys.argv) != 3: 
		usage()
	chrom_dir = sys.argv[1]
	helper.check_dir_exist(chrom_dir)
	out_dir = sys.argv[2]
	helper.make_dir(out_dir)
	print "Done getting command line arguments"

def usage():
	print "python divide_chromsome_data_into_window.py"
	print "chrom_dir: directory where the data are divided by chromosomes"
	print "out_dir: directory where the data that are divided by each 1MB or 10MB windows are stored" 
	print "this function is simply to split data from the chrom_dir into smaller windows"
	exit(1)
main()