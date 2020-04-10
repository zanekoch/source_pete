ct_list_fn=$1
out_dir=$2
##### CODE FILES #####
get_color_segment_code=/u/home/h/havu73/project-ernst/source/chromHMM_utilities/get_segmentation_file_with_color.py
######################
mkdir -p $out_dir
sh_fn=/u/home/h/havu73/project-ernst/source_pete/jobs/get_raw_segmentation_igv.sh
rm -f $sh_fn
raw_segmentation_folder=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/segmentation
region_to_extract_fn=/u/home/h/havu73/project-ernst/diff_pete/roadmap/raw_semgentation/chr5_15_region.bed
state_annot_fn=/u/home/h/havu73/project-ernst/data/roadmap_epigenome/18_core_K27ac_model_downloaded/state_annotation.txt
cg_list=""
while IFS= read -r ct; do
	cg_list="$cg_list $ct"
done < $ct_list_fn
echo $cg_list

for ct in $cg_list
do
	segmentation_fn=$raw_segmentation_folder/$ct/${ct}_18_core_K27ac_segments.bed.gz
	echo $segmentation_fn
	ct_outdir=$out_dir/$ct
	mkdir -p $ct_outdir
	ct_output_fn=$ct_outdir/${ct}_chr5_15_segment.bed
	command="bedtools intersect -a $segmentation_fn -b $region_to_extract_fn > $ct_output_fn"
	echo $command >> $sh_fn
	command="gzip $ct_output_fn"
	echo $command >> $sh_fn
	ct_output_fn=${ct_output_fn}.gz
	out_segment_color_fn=$ct_outdir/${ct}_chr5_15_igv.bed
	command="python $get_color_segment_code ${ct_output_fn} $state_annot_fn $out_segment_color_fn ${ct}"
	echo $command >> $sh_fn
	command="gzip $out_segment_color_fn"
	echo $command >> $sh_fn
	out_segment_color_fn=${out_segment_color_fn}.gz
done
chmod +x $sh_fn