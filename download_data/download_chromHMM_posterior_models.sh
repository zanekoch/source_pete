download_link_prefix="https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/POSTERIOR/"
s18_fn_phrase=18_core_K27ac
fn_suffix=posterior.txt.gz
ct_list_fn="/u/home/h/havu73/project-ernst/data/roadmap_epigenome/all_127_celltype"
s18_model_dir="/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_state_models"
s18_pos_dir=$s18_model_dir/posterior/
declare -a chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
write_list_of_cell_type=""
num_chromHMM_state=18
# while read ct 
# do
# 	# create a folder for the ct
# 	ct_outDir=$s18_pos_dir/$ct
# 	mkdir -p $ct_outDir
# 	for chrom_index in "${chrom_list[@]}"
# 	do
# 		link="${download_link_prefix}/${ct}_${s18_fn_phrase}_chr${chrom_index}_${fn_suffix}"
# 		command="wget -P $ct_outDir $link" # download the data to the right folder corresponding to the cell type of interest.
# 		$command
# 	done
# 	echo "Done downloading data for cell type: $ct"
# 	echo ""
#  	echo ""
#  	echo ""
# done < $ct_list_fn

