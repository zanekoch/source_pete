library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')
read_avg_roc_df <- function(fn, train_mode){
	print(fn)
  df <- (read.table(fn, sep = '\t', stringsAsFactors = F))
  temp_col_names <- df[1,]
  df <- df[-1,] # get rid of the first row
  colnames(df) <- temp_col_names
  df <- df %>% mutate_each(as.double)
  df$tm <- rep(train_mode, nrow(df))
  return(df)
}

get_avg_roc_fn <- function(train_mode, cell_group_dir){
	return(file.path(cell_group_dir, train_mode, 'average_roc_all_ct/summary_all_ct_roc.txt.gz'))
}

draw_roc_plot_compare_train_mode <- function(plot_df, out_dir, state, cell_group){
	save_fn <- file.path(out_dir, paste0("avg_roc_state", state, ".png"))
	plot_title <- paste0(cell_group, "_avg_roc_state", state)
	p <- ggplot(data = plot_df, aes_string(y = paste0('avg_tpr_S', state), x = paste0('avg_fpr_S', state), color = "tm")) +
	geom_line() + 
	theme_bw() +
	ggtitle(plot_title)
	ggsave(save_fn)
}

main <- function(num_train_mode, cell_group_dir, out_dir, train_mode_list, num_chromHMM_state, cell_group){
	# 1. get the list of files that contain the data of average roc for different training modes
	print(train_mode_list)
	roc_fn_list <- sapply(train_mode_list, get_avg_roc_fn, cell_group_dir = cell_group_dir) 
	roc_fn_list <- as.character(roc_fn_list)
	print(roc_fn_list)
	# 2. get the list of dataframs from the list of file names that we just got
	roc_df_list <- list()
	for (tm_index in seq(num_train_mode)){
		this_tm_df <- read_avg_roc_df(roc_fn_list[tm_index], train_mode_list[tm_index])
		roc_df_list[[tm_index]] <- this_tm_df
	}
	# loop through each state
	for (state in 1:num_chromHMM_state){
		columns_to_select <- c(paste0(c('avg_fpr_S', 'avg_tpr_S', 'std_fpr_S', 'std_tpr_S'), state), 'tm')
		state_roc_df_list <- lapply(roc_df_list, FUN = function(x) (x %>% select(columns_to_select))) # for each df, select only the ones associated with the state we are trying to create a plot for
		plot_df <- bind_rows(state_roc_df_list)
		draw_roc_plot_compare_train_mode(plot_df, out_dir, state, cell_group)
		print(paste0("Done with drawing plots for state", state))
	}
}

args = commandArgs(trailingOnly=TRUE)
NUM_SURE_ARGS = 5
if (length(args) < NUM_SURE_ARGS)
{
	stop("wrong command line argument format", call.=FALSE)
}
num_train_mode <- as.integer(args[1]) # output of ChromHMM OverlapEnrichment, also input to this program
if(length(args) != (NUM_SURE_ARGS + num_train_mode)){
	stop("Number of arguments do not match the number of specified train_mode", call.=FALSE)
}
cell_group_dir <- args[2] # where the figures should be stored
if (! dir.exists(cell_group_dir)){
	stop(paste0("cell_group_dir DOES NOT EXIST: ", cell_group_dir), call.=FALSE)
}
out_dir <- args[3] # name of the enrichment context, so that it can be title of the figure
dir.create(out_dir, recursive = TRUE)
num_chromHMM_state <- as.integer(args[4])
cell_group <- args[5]
train_mode_list <- args[(NUM_SURE_ARGS + 1):length(args)]
main(num_train_mode, cell_group_dir, out_dir, train_mode_list, num_chromHMM_state, cell_group)