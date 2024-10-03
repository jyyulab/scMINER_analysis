rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)




#-------
datasets_id <- c('Zeisel','PBMC14K')
for (dataset in datasets_id ){
  # heatmap
 
  summary_file <- list.files(pattern = paste0(dataset,'.*_bin_size_.*.csv'))
  eval_df <- read.table(summary_file)
  eval_df <- eval_df[eval_df$bin_size !=1,]
  #eval_df <- eval_df[eval_df$louvain_res <=3, ]
  mat <- reshape2::acast(eval_df, bin_size ~ louvain_res, value.var = "ARI")
  rownames(mat) <- as.numeric(rownames(mat))
  colnames(mat) <- as.numeric(colnames(mat))
  sorted_mat <- mat[order(rownames(mat)), order(colnames(mat))] 
  
  pdf(paste0(dataset,'_scMINER_sweep_bin_size_Heatmap.pdf'))
  p <- Heatmap(mat,
          heatmap_legend_param = list(title_gp = gpar(fontsize = 16), # Legend title font size  
                                      labels_gp = gpar(fontsize = 16)), # Legend labels font size  
          row_names_gp = gpar(fontsize = 14), # Row names font size  
          column_names_gp = gpar(fontsize = 14), # Column names font size  
          row_title = "bin size", row_title_side = "left", row_title_gp = gpar(fontsize = 24), # Row title font size  
          column_title = "resolution", column_title_side = "bottom", column_title_gp = gpar(fontsize = 24),
          cluster_rows = FALSE,cluster_columns = FALSE,
          name='ARI')
  print(p)
  dev.off()
}
#--------------

