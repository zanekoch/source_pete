#cg_dir=$1
cg_dir=/Users/vuthaiha/Desktop/window_hoff/diff_pete/roadmap/blood
###### CODE FILES ######
draw_hist_pos_code=/Users/vuthaiha/Desktop/window_hoff/source_pete/train_and_evaluate/draw_hist_pos_prob.R
########################
declare -a method_list=('baseline' 'chromHMM_posterior' 'multi_logistic')

small_bins_per_big_bin=50
num_chromHMM_state=18
for method in "${method_list[@]}"
do
	hist_dir=$cg_dir/${method}/representative_data/state_segment_histogram/
	hist_fn=${hist_dir}/prop_hist.txt.gz
	out_dir=${hist_dir}/bin_size_${small_bins_per_big_bin}/
	mkdir -p $out_dir
	command="Rscript $draw_hist_pos_code $hist_fn $out_dir $num_chromHMM_state $small_bins_per_big_bin"
	# echo $command
	echo "Done drawing for method $method"
done

hist_fn=/Users/vuthaiha/Desktop/window_hoff/data/roadmap_epigenome/18_core_K27ac_model_downloaded/histogram_posterior/blood/avgerage_posterior_histogram.txt.gz
blood_out_dir=/Users/vuthaiha/Desktop/window_hoff/diff_pete/roadmap/blood/hist_raw_segmentation/blood/bin_size_${small_bins_per_big_bin}
mkdir -p $blood_out_dir
command="Rscript $draw_hist_pos_code $hist_fn $blood_out_dir $num_chromHMM_state $small_bins_per_big_bin"

echo $command