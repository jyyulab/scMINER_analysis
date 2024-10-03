rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)



#-------
datasets_id <- c('Yan','Buettner','Pollen','Kolod','Goolam','Chung','Usoskin','Klein')



for (dataset in datasets_id ){
  # heatmap
 
  summary_file <- list.files(pattern = paste0(dataset,'.*_binsize_summary.csv'))
  eval_df <- read.table(summary_file)
  #break

  mat <- reshape2::acast(eval_df, component ~ bin_size, value.var = "ARI")
  rownames(mat) <- as.numeric(rownames(mat))
  colnames(mat) <- as.numeric(colnames(mat))

  # Define the color mapping  
  col_fun = colorRamp2(c(0, 1), c("blue", "red"))  
  
  pdf(paste0(dataset,'_scMINER_sweep_bin_size_component_Heatmap.pdf'))
  
  if (dataset == 'Kolod'){
  p <- Heatmap(mat,
               col = col_fun, 
          heatmap_legend_param = list(title_gp = gpar(fontsize = 16), # Legend title font size  
                                      labels_gp = gpar(fontsize = 16)), # Legend labels font size  
          row_names_gp = gpar(fontsize = 14), # Row names font size  
          column_names_gp = gpar(fontsize = 14), # Column names font size  
          row_title = "component", row_title_side = "left", row_title_gp = gpar(fontsize = 24), # Row title font size  
          column_title = "bin size", column_title_side = "bottom", column_title_gp = gpar(fontsize = 24),
          cluster_rows = FALSE,cluster_columns = FALSE,
          name='ARI')
  } else{
    p <- Heatmap(mat,
                 #col = col_fun, 
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 16), # Legend title font size  
                                             labels_gp = gpar(fontsize = 16)), # Legend labels font size  
                 row_names_gp = gpar(fontsize = 14), # Row names font size  
                 column_names_gp = gpar(fontsize = 14), # Column names font size  
                 row_title = "dimension", row_title_side = "left", row_title_gp = gpar(fontsize = 24), # Row title font size  
                 column_title = "bin size", column_title_side = "bottom", column_title_gp = gpar(fontsize = 24),
                 cluster_rows = FALSE,cluster_columns = FALSE,
                 name='ARI')
  }
  print(p)
  dev.off()
}
#--------------

