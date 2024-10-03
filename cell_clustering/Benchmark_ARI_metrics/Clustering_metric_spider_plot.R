# Radar/Spider plot for summary of ARI,AMI,NMI,ASW,AvgBIO
rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsignif)
library(ggpubr)
library(fmsb)
library(grDevices)

packages <- c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster')
ordered_dataset_ids <- c('Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein',
                         'Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K')
# mean_SD <- function(x) {  
#   var <- var(x)  
#   return(c(y = mean(x), ymin = mean(x) - sqrt(var), ymax = mean(x) + sqrt(var)))  
# }

# ECP ECA version
ARI_input_dir <- '.'
entropy_dir <- 'Cluster_entropy_figure/'

spider_mean_df <- data.frame(row.names = packages)
metrics <- c('ARI','AMI','NMI','AvgBIO')
for (method in packages){
  metric_file <- list.files(path=ARI_input_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  #print(metric_file)
  #metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  metric_df <- metric_df[1:10,metrics] # use PBMC14K and datasets before with real true label
  metric_mean <- colMeans(metric_df,na.rm=TRUE)
  spider_mean_df[method,1:length(metrics)] <- metric_mean
  # add ECP and ECA
  entropy_file <- list.files(path=entropy_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  entropy_df <- read.table(entropy_file,sep='\t',header=TRUE)
  ECP_mean <- mean(entropy_df$ECP)
  ECA_mean <- mean(entropy_df$ECA)
  spider_mean_df[method,length(metrics)+1] <- ECP_mean 
  spider_mean_df[method,length(metrics)+2] <- ECA_mean 
}
colnames(spider_mean_df) <- c(metrics,'ECP','ECA')




# add min 0 and max 1
spider_mean_df <- rbind(rep(1,ncol(spider_mean_df)),rep(0,ncol(spider_mean_df)),spider_mean_df)
spider_mean_df <- spider_mean_df[,c('ARI','AMI','ECA','AvgBIO','ECP','NMI')]
colors_method = c('#ff0000',"#AED581","#AED6F1","#8b5635","#fd9c32","#E1BEE7")
colors_method_opaque = grDevices::adjustcolor(colors_method , alpha.f = 0.1)  
# scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
#                              "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7"))

#
# pdf('ARI_metric_radar_plot_legendOff.pdf',width=7,height=7)
# par(cex=1.5)
# radarchart(spider_mean_df,
#            pcol=colors_method,
#            pfcol = colors_method_opaque,
#            plwd=1 , plty=1,
#            #custom the grid
#            cglcol="grey", cglty=1, axislabcol="black", caxislabels=c(0.5,1), cglwd=1.2,
#            #custom labels
#            vlcex=0.8)  
# dev.off()



#
pdf('ARI_metric_radar_plot_legendOn.pdf',width=7,height=7)
radarchart(spider_mean_df,
           pcol=colors_method,
           pfcol = colors_method_opaque,
           plwd=2 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", caxislabels=c(0.5,1), cglwd=1.2,
           #custom labels
           vlcex=0.8)  
legend("topright", legend = rownames(spider_mean_df[-c(1,2),]), 
       bty = "n", pch=20 , col=colors_method, text.col = "black", cex=1, pt.cex=3)
dev.off()

########################
#
# Use accuracy and purity instead of ECA and ECP
ARI_input_dir <- '.'
accuracy_purtiy_dir <- 'Cluster_Accuracy_Purity_figure/'

spider_mean_df <- data.frame(row.names = packages)
metrics <- c('ARI','AMI','NMI','AvgBIO')
for (method in packages){
  metric_file <- list.files(path=ARI_input_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  #print(metric_file)
  #metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  metric_df <- metric_df[1:10,metrics] # use PBMC14K and datasets before with real true label
  metric_mean <- colMeans(metric_df,na.rm=TRUE)
  spider_mean_df[method,1:length(metrics)] <- metric_mean
  # add accuracy and purity
  accuracy_purtiy_file <- list.files(path=accuracy_purtiy_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  accuracy_purtiy_df <- read.table(accuracy_purtiy_file,sep='\t',header=TRUE)
  accuray_mean <- mean(accuracy_purtiy_df$accuracy)
  purity_mean <- mean(accuracy_purtiy_df$purity)
  spider_mean_df[method,length(metrics)+1] <- accuray_mean
  spider_mean_df[method,length(metrics)+2] <- purity_mean
}
colnames(spider_mean_df) <- c(metrics,'Accuracy','Purity')

# Color scheme
spider_mean_df <- rbind(rep(1,ncol(spider_mean_df)),rep(0,ncol(spider_mean_df)),spider_mean_df)
spider_mean_df <- spider_mean_df[,c('ARI','AMI','Accuracy','AvgBIO','Purity','NMI')]
colors_method = c('#ff0000',"#AED581","#AED6F1","#8b5635","#fd9c32","#E1BEE7")
colors_method_opaque = grDevices::adjustcolor(colors_method , alpha.f = 0.1)  

# spider Plot
pdf('ARI_accuracy_purity_metric_radar_plot_legendOn.pdf',width=7,height=7)
radarchart(spider_mean_df,
           pcol=colors_method,
           pfcol = colors_method_opaque,
           plwd=2 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", caxislabels=c(0.5,1), cglwd=1.2,
           #custom labels
           vlcex=0.8)  
legend("topright", legend = rownames(spider_mean_df[-c(1,2),]), 
       bty = "n", pch=20 , col=colors_method, text.col = "black", cex=1, pt.cex=3)
dev.off()


############
# Use ASW with true label
# Use accuracy and purity instead of ECA and ECP
ARI_input_dir <- '.'
accuracy_purtiy_dir <- 'Cluster_Accuracy_Purity_figure/'

spider_mean_df <- data.frame(row.names = packages)
metrics <- c('ARI','AMI','NMI','ASW')
for (method in packages){
  metric_file <- list.files(path=ARI_input_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  #print(metric_file)
  #metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  metric_df <- metric_df[1:10,metrics] # use PBMC14K and datasets before with real true label
  metric_mean <- colMeans(metric_df,na.rm=TRUE)
  spider_mean_df[method,1:length(metrics)] <- metric_mean
  # add accuracy and purity
  accuracy_purtiy_file <- list.files(path=accuracy_purtiy_dir,pattern = paste0(method,'.*.csv'),full.names = TRUE)
  accuracy_purtiy_df <- read.table(accuracy_purtiy_file,sep='\t',header=TRUE)
  accuray_mean <- mean(accuracy_purtiy_df$accuracy)
  purity_mean <- mean(accuracy_purtiy_df$purity)
  spider_mean_df[method,length(metrics)+1] <- accuray_mean
  spider_mean_df[method,length(metrics)+2] <- purity_mean
}
colnames(spider_mean_df) <- c(metrics,'Accuracy','Purity')

# Color scheme
spider_mean_df <- rbind(rep(1,ncol(spider_mean_df)),rep(0,ncol(spider_mean_df)),spider_mean_df)
spider_mean_df <- spider_mean_df[,c('ARI','AMI','Accuracy','ASW','Purity','NMI')]
colors_method = c('#ff0000',"#AED581","#AED6F1","#8b5635","#fd9c32","#E1BEE7")
colors_method_opaque = grDevices::adjustcolor(colors_method , alpha.f = 0.1)  

# spider Plot
pdf('ARI_ASW_accuracy_purity_metric_radar_plot_legendOn.pdf',width=7,height=7)
radarchart(spider_mean_df,
           pcol=colors_method,
           pfcol = colors_method_opaque,
           plwd=2 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", caxislabels=c(0.5,1), cglwd=1.2,
           #custom labels
           vlcex=0.8)  
legend("topright", legend = rownames(spider_mean_df[-c(1,2),]), 
       bty = "n", pch=20 , col=colors_method, text.col = "black", cex=1, pt.cex=3)
dev.off()


