library("tidyverse")
library("tidyr")
library("dplyr")
library('ggplot2')

draw_hist_plot <- function(fn, num_chromHMM_state, small_bins_per_big_bin, y_col_to_draw, save_fn){
	df <- as.data.frame(read.table(fn, header = TRUE, sep = '\t'))
	colnames(df) <- c(paste0('E', 1:num_chromHMM_state), 'all_states')
	RMIN <- 0.0
	RMAX <- 1.0
	NBIN <- 1000
	# cat(NBIN)
	# cat("\n")
	# cat(small_bins_per_big_bin)
	# cat("\n")
	# cat(NBIN / small_bins_per_big_bin)
	# cat("\n")
	num_big_bins <- NBIN / small_bins_per_big_bin
	bin_size <- (RMAX - RMIN) / num_big_bins
	re_df <- data.frame(matrix(ncol = ncol(df), nrow = 0))
	colnames(re_df) <- colnames(df)
	for (i in (1:num_big_bins)){
		this_bin_df <- df[(1 + (i-1) * small_bins_per_big_bin): (i * small_bins_per_big_bin),]
		big_bin_data <- rowsum(this_bin_df, group = rep(1, nrow(this_bin_df)))
		re_df[i,] <- big_bin_data
	}
	re_df$x_values <- bin_size * (1:num_big_bins)
	p <- ggplot(data = re_df, aes_string(x = 'x_values', y = y_col_to_draw)) +
	geom_bar(stat = 'identity') +
	ylim(0,1) +
	theme_bw() +
	ggtitle(y_col_to_draw)
	ggsave(save_fn)
	return(p)
}


args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4)
{
	stop("wrong command line argument format", call.=FALSE)
}
args = commandArgs(trailingOnly=TRUE)
hist_fn <- args[1] # output of ChromHMM OverlapEnrichment, also input to this program
save_dir <- args[2] # where the figures should be stored
dir.create(save_dir, recursive = TRUE)
num_chromHMM_state <- args[3] # where data of states' mnenomics and colors are stored
small_bins_per_big_bin <- as.integer(args[4])
colnames_to_draw <- c(paste0('E', 1:num_chromHMM_state), 'all_states')
cat ("Drawing figures for all columns in the table")
for (cn in colnames_to_draw){
	save_fn <- file.path(save_dir, paste0('hist_', cn, '.png'))
	draw_hist_plot(hist_fn, num_chromHMM_state, small_bins_per_big_bin, cn, save_fn)
}
cat("Done!")