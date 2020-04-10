	##### GET COMMANDLINE ARGUMENTS #####
cell_group=$1
train_mode=$2 # ['baseline', 'multi_logistic', 'chromHMM_posterior']
#####################################
##### CODE FILES ####
cross_validation_segment_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/segment_based/run_cross_validation_pipeline_segment_based.py
cross_validation_posterior_code=/u/home/h/havu73/project-ernst//source_pete/train_and_evaluate/posterior_based/run_cross_validation_pipeline_posterior_based.py
calculate_pred_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/calculate_regression_performace.py
draw_confusion_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/draw_confusion_matrix_pred_results.R
calculate_roc_code=/u/home/h/havu73/project-ernst/source_pete/train_and_evaluate/calculate_ROC_one_state.py
#####################
roadmap_outDir=/u/home/h/havu73/project-ernst/diff_pete/roadmap/
all_ct_posterior_folder=/u/home/h/havu73/project-ernst//data/roadmap_epigenome/18_core_K27ac_model_downloaded/posterior/
cg_outDir=$roadmap_outDir/${cell_group}/${train_mode}
mkdir -p $cg_outDir
cell_group_fn=${roadmap_outDir}/${cell_group}.list
train_sampled_data_fn=/u/home/h/havu73/project-ernst/diff_pete/roadmap/${cell_group}/all_replicates_segmentation.bed.gz
all_ct_segment_folder=/u/home/h/havu73/project-ernst/diff_pete/roadmap/all_ct_segments
num_chromHMM_state=18
shDir=/u/home/h/havu73/project-ernst/source_pete/jobs
draw_job_fn=${shDir}/draw_plots
#rm -f $draw_job_fn
shPrefix=cross_validation_${cell_group}
job_index=1
# first, we will get all the cell types related to ${cell_group}
ct_list=""
while read line;  
do
	ct_list="$ct_list $line"
done < $cell_group_fn

# now create the validation jobs for each cell type as the validation cell type
for ct in $ct_list
do
	sh_file_name=${shDir}/${shPrefix}_${job_index}.sh
	rm -f $sh_file_name # so that we can create new code file
	this_ct_outdir=${cg_outDir}/val_${ct}/average_predictions/
	this_ct_pred_performance_fn=${cg_outDir}/val_${ct}/prediction_performance
	if [ ! -f $this_ct_pred_performance_fn ]
	then
		echo "cell_group: $cell_group . Cell_type: $ct . validation of this cell type is not available"
		if [ $train_mode == "baseline" ] || [ $train_mode == "multi_logistic" ] # if baseline or multi_logistic using segmentation data, then use one piece of code
		then
			command="python $cross_validation_segment_code ${train_sampled_data_fn} ${cg_outDir} ${all_ct_segment_folder} ${num_chromHMM_state} $ct $train_mode $cell_group_fn"
			echo $command > $sh_file_name
		fi
		if [ $train_mode == "chromHMM_posterior" ] # if train using multi_logistic regression using chromHMM posterior data, use a different piece of code
		then 
			command="python $cross_validation_posterior_code ${train_sampled_data_fn} ${cg_outDir} ${all_ct_posterior_folder} ${num_chromHMM_state} $ct $cell_group_fn"
			echo $command > $sh_file_name
		fi
		command="python $calculate_pred_code $this_ct_outdir $all_ct_segment_folder $this_ct_pred_performance_fn $ct $num_chromHMM_state"
		echo $command >> $sh_file_name
		roc_outDir=${cg_outDir}/val_${ct}/roc/
		num_score_bins=100
		mkdir -p $roc_outDir
		command="python $calculate_roc_code $this_ct_outdir $all_ct_segment_folder $ct $roc_outDir $num_chromHMM_state $num_score_bins"
		echo $command >> $sh_file_name
		chmod +x $sh_file_name
		job_index=$(($job_index + 1))
	else
		echo ""
		echo ""
		echo ""
		state_annot_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/state_annotation.txt
		confusion_plot=${cg_outDir}/val_${ct}/confusion_matrix.png
		command="Rscript $draw_confusion_code $this_ct_pred_performance_fn $confusion_plot $state_annot_fn"
		echo $command >> $draw_job_fn
	fi 

	echo Done with cell type: $ct
done

# chmod +x $draw_job_fn


