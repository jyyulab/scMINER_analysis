#------
rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(mclust)
library(ggsankey)
library(reshape2)
library(dplyr)
library(ggsignif)
library(ggpubr)
library(pdfCluster)

# set up datasets
# --------------
ordered_dataset_ids = c('Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein',
                        'Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K')
packages = c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster')

#---------------------
dir.create('Cluster_Accuracy_Purity_figure')

mean_SD <- function(x) {  
  var <- var(x)  
  return(c(y = mean(x), ymin = mean(x) - sqrt(var), ymax = mean(x) + sqrt(var)))  
}

# cluster purity
get_cluster_purity <- function(true_label,pred_label){
  if(length(true_label) != length(pred_label)){
    cat('Error, true label and pred label not equal length \n')
    return(NA)
  }
  
  df <- data.frame(True = true_label, Pred = pred_label)      
  
  pred_label_groups <- split(df, df$Pred)      
  
  majority_counts <- sapply(pred_label_groups, function(x) {      
    max(table(x$True))      
  })      
  
  pred_label_counts <- table(pred_label)    
  
  purity <- sum((majority_counts / pred_label_counts)*pred_label_counts / length(pred_label))   
  
  # Return purity     
  return(purity)   
}

# cluster accuracy
get_cluster_accuracy <- function(true_label,pred_label){
  if(length(true_label) != length(pred_label)){
    cat('Error, true label and pred label not equal length \n')
    return(NA)
  }
  df <- data.frame(True = true_label, Pred = pred_label)    
  
  true_label_groups <- split(df, df$True)    
  
  majority_counts <- sapply(true_label_groups, function(x) {    
    max(table(x$Pred))    
  })    
  
  true_label_counts <- table(true_label)  
  
  accuracy <- sum((majority_counts / true_label_counts)*true_label_counts / length(true_label))  
  
  # Return accuracy   
  return(accuracy)  
}
#----------


#------------
# start calculation of cluster Accuracy & Purity
#-----------------  
for (method in packages){
  #if (method != 'scMINER'){next}
  acc_purity_df = data.frame(matrix(ncol = 3, nrow = 0)) 
  for (dataset in ordered_dataset_ids){
    #if(dataset != 'PBMC14K'){next}
    #if (dataset %in% c('HMC76K','Covid97K','Covid650K')){next} # no "true" label 
    
    print(dataset)
    true_label_file <- list.files(path="../True_labels",pattern=dataset,full.names=TRUE)
    #print(true_label_file)
    true_label_df <- read.table(true_label_file,header=TRUE,sep='\t')
    
    pred_label_file <- list.files(path=paste0("Clustering_output/",method),
                                  pattern=paste0(dataset,'.*.txt'),full.names=TRUE)
    #print(pred_label_file)
    pred_label_df <- read.table(pred_label_file,header=TRUE,sep='\t') 
    #print(colnames(pred_label_df))
    if (method == 'scMINER'){
      if (dataset != 'Covid650K'){
        merge_df <- merge(pred_label_df,true_label_df,by.x ='ID',by.y = 'cell')
        colnames(merge_df)[colnames(merge_df) == 'label'] <-  "scMINER_label"
      }else{
        merge_df <- cbind(pred_label_df,true_label_df)
        colnames(merge_df)[colnames(merge_df) == 'label'] <-  "scMINER_label"
      }
    } else {
      merge_df <- merge(pred_label_df,true_label_df,by.x="cell_id",by.y = 'cell')
    }
    #print(colnames(merge_df))  
    #print(nrow(merge_df))
    #break
    
    # acc and purity 
    pred_label_key <- paste0(method,'_label')
    #print(pred_label_key)
    #print(unique(merge_df$true_label))
    #print(unique(merge_df[[pred_label_key]]))
    #break
    avg_accuracy <- round(get_cluster_accuracy(merge_df$true_label,merge_df[[pred_label_key]]),3) 
    avg_purity <- round(get_cluster_purity(merge_df$true_label,merge_df[[pred_label_key]]),3)


    n_cell <- nrow(merge_df)
    new_row <- c(dataset,n_cell,method,avg_accuracy,avg_purity)
    acc_purity_df <- rbind(acc_purity_df,new_row)
    
  }
  colnames(acc_purity_df) <- c('dataset_id','n_cell','method','accuracy','purity')
  filename = paste0('Cluster_Accuracy_Purity_figure/',method,'_avg_accuracy_purity.csv')
  write.table(acc_purity_df,file=filename,sep = '\t',row.names = FALSE)  
  #break
}



######################
# barplot of Accuracy & Purity
accuracy_df <- data.frame(row.names = ordered_dataset_ids)
purity_df <- data.frame(row.names = ordered_dataset_ids)
input_dir='./Cluster_Accuracy_Purity_figure/'
output_dir='./Cluster_Accuracy_Purity_figure/'
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  accuracy_col <- metric_df[['accuracy']]
  print(accuracy_col)
  accuracy_df[1:length(accuracy_col),method] <- accuracy_col

  # purity
  purity_col <- metric_df[['purity']]
  print(purity_col)
  purity_df[1:length(purity_col),method] <- purity_col
}



# Accuracy
# skip Covid650K temporarily
accuracy_df <- accuracy_df[1:10,] # stops at PBMC14K

# prepare for melt
accuracy_df$dataset_ID <- rownames(accuracy_df)
accuracy_df_long <- melt(accuracy_df,id.vars='dataset_ID')
colnames(accuracy_df_long) <- c('dataset_ID','method','accuracy')
accuracy_df_long$dataset_ID <- factor(accuracy_df_long$dataset_ID, levels = unique(accuracy_df_long$dataset_ID)) 
accuracy_df_long$method <- factor(accuracy_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))

# grouped barplot accuracy
accuracy_df_long[is.na(accuracy_df_long)] <-0 
accuracy_plot <- ggplot(accuracy_df_long,aes(x= dataset_ID,y= accuracy,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(accuracy==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "Accuracy", fill = "")  
print(accuracy_plot)
ggsave(paste0(output_dir,"scMINER_clustering_accuracy.pdf"),plot=accuracy_plot)

# mean accuracy + variance errorbar + wilcox.test of
#accuracy_df_long$ari <- accuracy_df_long$ari + runif(nrow(accuracy_df_long), min = -1e-10, max = 1e-10)  # break tied value
accuracy_summary_plot <- ggplot(accuracy_df_long, aes(x = method,y=accuracy,fill=method)) +
  
  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +  
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.9)) +
  
  geom_signif(
    comparisons = list(c("scMINER","Seurat"),c("scMINER","SC3s"),
                       c("scMINER","Scanpy"),c("scMINER","scVI"),
                       c("scMINER","scDeepCluster")),
    test = "wilcox.test",test.args = list(paired = TRUE, alternative = "two.sided"),
    map_signif_level = TRUE,y_position = c(1,1.1,1.2,1.3,1.4),textsize = 5
  ) +
  
  theme(aspect.ratio = 1,
        text = element_text(size = 24),
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
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),  
                     labels = seq(0, 1, by = 0.2)) +  
  
  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) + theme(legend.position ="none")+
  labs(x = "", y = "Accuracy", fill = "")
print(accuracy_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_accuracy_wilcoxon.pdf"),plot=accuracy_summary_plot)


# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_accuracy_Wilcoxon <- data.frame(  
  method = character(),   
  statistic = numeric(),   
  p.value = numeric(),   
  stringsAsFactors=FALSE  
) 

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(accuracy_df_long[accuracy_df_long$method == 'scMINER','accuracy'],
                             accuracy_df_long[accuracy_df_long$method == method,'accuracy'],
                             alternative="greater",paired=TRUE)
  scMINER_accuracy_Wilcoxon <- rbind(  
    scMINER_accuracy_Wilcoxon,   
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)  
  ) 
}

write.table(scMINER_accuracy_Wilcoxon,file=paste0(output_dir,"scMINER_accuracy_Wilcoxon.txt"),sep="\t",row.names = FALSE)

##############################
# grouped barplot purity
purity_df <- purity_df[1:10,] # stops at PBMC14K

# prepare for melt
purity_df$dataset_ID <- rownames(purity_df)
purity_df_long <- melt(purity_df,id.vars='dataset_ID')
colnames(purity_df_long) <- c('dataset_ID','method','purity')
purity_df_long$dataset_ID <- factor(purity_df_long$dataset_ID, levels = unique(purity_df_long$dataset_ID)) 
purity_df_long$method <- factor(purity_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))



purity_df_long[is.na(purity_df_long)] <-0 
purity_plot <- ggplot(purity_df_long,aes(x= dataset_ID,y= purity,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(purity==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "Purity", fill = "")  
print(purity_plot)
ggsave(paste0(output_dir,"scMINER_clustering_purity.pdf"),plot=purity_plot)

# mean purity + variance errorbar + wilcox.test of
#purity_df_long$ari <- purity_df_long$ari + runif(nrow(purity_df_long), min = -1e-10, max = 1e-10)  # break tied value
purity_summary_plot <- ggplot(purity_df_long, aes(x = method,y=purity,fill=method)) +
  
  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +  
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.9)) +
  
  geom_signif(
    comparisons = list(c("scMINER","Seurat"),c("scMINER","SC3s"),
                       c("scMINER","Scanpy"),c("scMINER","scVI"),
                       c("scMINER","scDeepCluster")),
    test = "wilcox.test",test.args = list(paired = TRUE, alternative = "two.sided"),
    map_signif_level = TRUE,y_position = c(1,1.1,1.2,1.3,1.4),textsize = 5
  ) +
  
  theme(aspect.ratio = 1,
        text = element_text(size = 24),
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
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),  
                     labels = seq(0, 1, by = 0.2)) +  
  
  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) + theme(legend.position ="none")+
  labs(x = "", y = "Purity", fill = "")
print(purity_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_purity_wilcoxon.pdf"),plot=purity_summary_plot)


# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_purity_Wilcoxon <- data.frame(  
  method = character(),   
  statistic = numeric(),   
  p.value = numeric(),   
  stringsAsFactors=FALSE  
) 

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(purity_df_long[purity_df_long$method == 'scMINER','purity'],
                             purity_df_long[purity_df_long$method == method,'purity'],
                             alternative="greater",paired=TRUE)
  scMINER_purity_Wilcoxon <- rbind(  
    scMINER_purity_Wilcoxon,   
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)  
  ) 
}

write.table(scMINER_purity_Wilcoxon,file=paste0(output_dir,"scMINER_purity_Wilcoxon.txt"),sep="\t",row.names = FALSE)



# plot Accuracy vs Purity for each dataset
#-------------------
df_all <- data.frame()
for (method in packages){
  
  accuracy_purity_file <- list.files(path = 'Cluster_Accuracy_Purity_figure/',
                             pattern=paste0(method,'_avg.*.csv'),full.names=TRUE)
  accuracy_purity_df <- read.table(accuracy_purity_file,header = TRUE)
  
  df_all <- rbind(df_all, accuracy_purity_df)
  #break
}
df_all$method <- as.factor(df_all$method)
df_all$method <- factor(df_all$method, levels = c("scMINER","Seurat","SC3s","Scanpy","scVI","scDeepCluster"))
levels(df_all$method)
# #------------


for (each in ordered_dataset_ids){
  #print(each)
  if (each %in% c('HMC76K','Covid97K','Covid650K')){next} # no "true" label 
  current_df = df_all[df_all$dataset_id == each,]

  p <- ggplot(current_df,aes(x=purity,y=accuracy,color=method,shape=method)) +
    geom_point(size = 10) +  
    scale_x_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) +  
    scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) + 
    labs(x = 'Purity', y= 'Accuracy',color='method',shape='method') +
    scale_color_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                                  "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) +
    ggtitle(paste0(each,' dataset')) + 
    theme(aspect.ratio = 1,
          axis.line = element_line(colour = "black",linewidth = 1),
          legend.key = element_blank(),
          plot.margin = margin(4,4,4,4,"cm"),
          text = element_text(size = 32),
          axis.text = element_text(size = 32),
          axis.text.x = element_text(size = 32,margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 32,margin = margin(t = 0, r = 10, b = 0, l = 0)),      
          panel.background = element_blank(),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.border = element_blank(),  
          axis.ticks.length = unit(-0.3,"cm")   
    )  
  print(p)
  #break
  plot_file <- paste0('Cluster_Accuracy_Purity_figure/',each,'_Accuracy_Purity.pdf')
  ggsave(plot_file,plot=p)
  #break
}
#---------------

# prepare manuscript figuer with 4 datset highlighted
highlight_df <- df_all[df_all$dataset_id %in% c('Buettner','Chung','Klein','PBMC14K'),]

ggplot(highlight_df, aes(x = purity, y = accuracy)) +
  geom_point(aes(color = method), size = 5) + # Scatter plot points with color
  scale_x_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) +  
  scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) + 
  facet_wrap(~ dataset_id,ncol = 4) + # Facet by dataset_id
  labs(
       x = "Purity",
       y = "Accuracy") +
  scale_color_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                                "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) +

  theme(aspect.ratio = 1,
        axis.line = element_line(colour = "black",linewidth = 1),
        legend.key = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 20,margin = margin(t = 0, r = 10, b = 0, l = 0)),      
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(),  
        axis.ticks.length = unit(-0.3,"cm")   
  )  
ggsave(filename = 'Cluster_Accuracy_Purity_figure/scMINER_accuracy_purity_Main4datasets.pdf',width = 16,height = 6)

# ED figure

highlight_df <- df_all[df_all$dataset_id %in% c('Yan','Goolam','Pollen','Usoskin','Kolod','Zeisel'),]
highlight_df$dataset_id <- factor(highlight_df$dataset_id,levels = c('Yan','Goolam','Pollen','Usoskin','Kolod','Zeisel'))

ggplot(highlight_df, aes(x = purity, y = accuracy)) +
  geom_point(aes(color = method), size = 5) + # Scatter plot points with color
  scale_x_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) +  
  scale_y_continuous(limits = c(0.4, 1.1), breaks = c(0.5, 0.75, 1)) + 
  facet_wrap(~ dataset_id,ncol = 3) + # Facet by dataset_id
  labs(
    x = "Purity",
    y = "Accuracy") +
  scale_color_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                                "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) +
  
  theme(aspect.ratio = 1,
        axis.line = element_line(colour = "black",linewidth = 1),
        legend.key = element_blank(),
        plot.margin = margin(2,2,2,2,"cm"),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 20,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 20,margin = margin(t = 0, r = 10, b = 0, l = 0)),      
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(),  
        axis.ticks.length = unit(-0.3,"cm")   
  )  
ggsave(filename = 'Cluster_Accuracy_Purity_figure/scMINER_accuracy_purity_ED6datasets.pdf',width = 16,height = 8)
