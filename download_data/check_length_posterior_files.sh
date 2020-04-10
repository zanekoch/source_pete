# this file will check that the length of all the posterior files that we downloaded are the same for different cell types, given that we are comparing files representing the same chromosome. 
# This step is necessary because it helps us ensure that when we predict results of posterior of states for different cell types, each line in the raw data (chromHMM posterior) corresponds to the same position on the genome
posterior_dir=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior
declare -a chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
ct_list=""
for f in $posterior_dir/E*/
do
	ct=$(echo $f | awk -F'/' '{print $(NF-1)}')
	ct_list="$ct $ct_list"
done

fn_middle='18_core_K27ac_chr'
fn_suffix='posterior.txt.gz'
for chrom_index in ${chrom_list[@]}
do
	echo "chrom_index: $chrom_index"
	for ct in $ct_list
	do
		fn=$posterior_dir/${ct}/${ct}_${fn_middle}${chrom_index}_${fn_suffix}
		num_line=$(zcat $fn | wc -l | awk '{print $1}')
		echo $ct $num_line
	done
	echo ""
	echo ""
	echo ""
done

chrom_len_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities/get_chrom_length_from_segmentation.sh

