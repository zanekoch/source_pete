train_mode=$1 # ['baseline', 'multi_logistic', 'chromHMM_posterior']
train_mode_short=$2 # used to write into ucsc file: base, norm, cPos
#####################################
##### CODE FILES ####
cross_validation_segment_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/segment_based/run_cross_validation_pipeline_segment_based.py
cross_validation_posterior_code=/u/home/h/havu73/project-ernst//source_pete/train_and_evaluate/posterior_based/run_cross_validation_pipeline_posterior_based.py
calculate_pred_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/calculate_regression_performace.py
draw_confusion_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/draw_confusion_matrix_pred_results.R
calculate_roc_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/calculate_ROC_one_state.py
sum_roc_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/summarize_roc_data_for_cell_group.py
get_rep_state_segment_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/get_representative_segmentation_cell_group.py
calculate_histogram_state_segment_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/count_histogram_pos_prob_one_cg.py
#####################
shDir=/u/home/h/havu73/project-ernst/source_pete/jobs
sh_fn=${shDir}/sum_stats_cg_${train_mode}
rm -f $sh_fn # so that we can create a new one
roadmap_outDir=/u/home/h/havu73/project-ernst/diff_pete/roadmap/
num_chromHMM_state=18
num_score_bins=100
state18_annot_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/state_annotation.txt
# get the list of cell groups
cg_list=""
for f in $roadmap_outDir/*.list
do
	cg=$(echo $f | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}')
	cg_list="$cg_list $cg"
done

for cg in $cg_list
do
	cg_dir=${roadmap_outDir}/${cg}/${train_mode}
	out_dir=${cg_dir}/average_roc_all_ct/
	ct_list_fn=${roadmap_outDir}/${cg}.list
	# command="python $sum_roc_code $cg_dir $out_dir $num_chromHMM_state $num_score_bins $ct_list_fn"
	# command="python $get_rep_state_segment_code $cg_dir/ $cg_dir/representative_data $state18_annot_fn $ct_list_fn $num_chromHMM_state $train_mode_short"
	command="python $calculate_histogram_state_segment_code $cg_dir/representative_data/avg_segmentation/ $cg_dir/representative_data/state_segment_histogram/ $num_chromHMM_state"
	echo $command >> $sh_fn
	command="rm $cg_dir/representative_data/state_segment_histogram/hist_draft_*.txt.gz" 
	echo $command >> $sh_fn
	#echo "echo  " >> $sh_fn
done
chmod +x $sh_fn