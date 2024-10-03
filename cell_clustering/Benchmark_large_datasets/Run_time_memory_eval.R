rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsignif)
library(ggpubr)
library(circlize)
library(tidyverse)
library(tibble)

packages <- c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster')

ordered_dataset_ids <- c('Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein',
                         'Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K')

large_dataset_ids <- c('PBMC14K','HMC76K','Covid97K','Covid650K')

output_dir <- './'
# Run time (40 CPU on HPC)
# -----------

time_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(pattern = paste0(method,'.*.txt'))
  #print(metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  time_col <- as.numeric(gsub('\\D','',metric_df[['RunTime']]))
  time_col <- time_col/60 # seconds to minute
  print(time_col)
  #break
  time_df[1:length(time_col),method] <- time_col
}

# prepare for melt
time_df <- time_df[1:10,] # stops at PBMC14K
time_df$dataset_ID <- rownames(time_df)
time_df_long <- melt(time_df,id.vars='dataset_ID')
colnames(time_df_long) <- c('dataset_ID','method','time')
time_df_long$dataset_ID <- factor(time_df_long$dataset_ID, levels = unique(time_df_long$dataset_ID)) 
time_df_long$method <- factor(time_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))

# grouped barplot time
time_df_long[is.na(time_df_long)] <-0 
time_plot <- ggplot(time_df_long,aes(x= dataset_ID,y= time,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(time==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
  theme(text = element_text(size = 24),aspect.ratio = 3/5,
        plot.margin = margin(2,2,2,2,"cm"),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),angle = 45, hjust = 1),
        axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(linewidth = 1),
        axis.ticks.length = unit(-0.1,"cm")) +  
  
  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) +
  
  labs(x = "", y = "time(min)", fill = "")  
print(time_plot)
ggsave(paste0(output_dir,"scMINER_clustering_time.pdf"),plot=time_plot)




# Max memory
#----------
memory_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(pattern = paste0(method,'.*.txt'))
  #print(metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  memory_col <- as.numeric(gsub('\\D','',metric_df[['MaxMemory']]))
  memory_col <- memory_col/1024 # convert MB to GB
  print(memory_col)
  #break
  memory_df[1:length(memory_col),method] <- memory_col
}

# prepare for melt
memory_df <- memory_df[1:10,] # stops at PBMC14K
memory_df$dataset_ID <- rownames(memory_df)
memory_df_long <- melt(memory_df,id.vars='dataset_ID')
colnames(memory_df_long) <- c('dataset_ID','method','memory')
memory_df_long$dataset_ID <- factor(memory_df_long$dataset_ID, levels = unique(memory_df_long$dataset_ID)) 
memory_df_long$method <- factor(memory_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))

# grouped barplot memory
memory_df_long[is.na(memory_df_long)] <-0 
memory_plot <- ggplot(memory_df_long,aes(x= dataset_ID,y= memory,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(memory==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
  theme(text = element_text(size = 24),aspect.ratio = 3/5,
        plot.margin = margin(2,2,2,2,"cm"),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),angle = 45, hjust = 1),
        axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(linewidth = 1),
        axis.ticks.length = unit(-0.1,"cm")) +  
  
  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) +
  
  labs(x = "", y = "memory(GB)", fill = "")  
print(memory_plot)
ggsave(paste0(output_dir,"scMINER_clustering_memory.pdf"),plot=memory_plot)





