ct_list_fn=$1
####
caclulate_hist_one_ct_code=/u/home/h/havu73/project-ernst/source_pete/process_raw_data/calculate_hist_pos_chromHMM_one_ct.py
####
posterior_dir=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior
hist_posterior_dir=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/histogram_posterior/
sh_fn=/u/home/h/havu73/project-ernst/source_pete/jobs/calculate_histogram_raw_chromHMM
rm -f $sh_fn # to write new file
mkdir -p $hist_posterior_dir
cg_list=""
while IFS= read -r ct; do
	cg_list="$cg_list $ct"
done < $ct_list_fn
echo $cg_list
for ct in $cg_list
do
	ct_posterior_dir=$posterior_dir/${ct}
	out_dir=${hist_posterior_dir}/${ct}
	mkdir -p $out_dir
	num_chromHMM_state=18
	prefix_pos_fn=${ct}_18_core_K27ac
	suffix_pos_fn=posterior.txt.gz
	command="python $caclulate_hist_one_ct_code $ct_posterior_dir $out_dir $num_chromHMM_state $ct $prefix_pos_fn $suffix_pos_fn"
	echo $command >> $sh_fn
done
chmod +x $sh_fn
