rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsignif)
library(ggpubr)

packages <- c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster')
ordered_dataset_ids <- c('Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein',
                       'Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K')

# mean, +/- 1 SD
mean_SD <- function(x) {  
  var <- var(x)  
  return(c(y = mean(x), ymin = mean(x) - sqrt(var), ymax = mean(x) + sqrt(var)))  
}

input_dir <- './'
output_dir <- 'ARI_figure/'
dir.create(output_dir)
# ARI
# -----------
ari_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  ari_col <- metric_df[['ARI']]
  print(ari_col)
  
  ari_df[1:length(ari_col),method] <- ari_col
}

# skip Covid650K temporarily
ari_df <- ari_df[1:10,] # stops at PBMC14K
# prepare for melt
ari_df$dataset_ID <- rownames(ari_df)
ari_df_long <- melt(ari_df,id.vars='dataset_ID')
colnames(ari_df_long) <- c('dataset_ID','method','ari')
ari_df_long$dataset_ID <- factor(ari_df_long$dataset_ID, levels = unique(ari_df_long$dataset_ID)) 
ari_df_long$method <- factor(ari_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))

# grouped barplot ari
ari_df_long[is.na(ari_df_long)] <-0 
ari_plot <- ggplot(ari_df_long,aes(x= dataset_ID,y= ari,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(ari==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "ARI", fill = "")  
print(ari_plot)
ggsave(paste0(output_dir,"scMINER_clustering_ARI.pdf"),plot=ari_plot)

# mean ari + variance errorbar + wilcox.test of
#ari_df_long$ari <- ari_df_long$ari + runif(nrow(ari_df_long), min = -1e-10, max = 1e-10)  # break tied value
ari_summary_plot <- ggplot(ari_df_long, aes(x = method,y=ari,fill=method)) +
  
  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +  
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.5)) +
  
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
  labs(x = "", y = "ARI", fill = "")
print(ari_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_ARI_wilcoxon.pdf"),plot=ari_summary_plot)


# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_ari_Wilcoxon <- data.frame(  
  method = character(),   
  statistic = numeric(),   
  p.value = numeric(),   
  stringsAsFactors=FALSE  
) 

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(ari_df_long[ari_df_long$method == 'scMINER','ari'],
                                  ari_df_long[ari_df_long$method == method,'ari'],
                                  alternative="two.sided",paired=TRUE)
  scMINER_ari_Wilcoxon <- rbind(  
    scMINER_ari_Wilcoxon,   
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)  
  ) 
}

write.table(scMINER_ari_Wilcoxon,file=paste0(output_dir,"scMINER_ARI_Wilcoxon.txt"),sep="\t",row.names = FALSE)

# ----------------------

# AMI
#---------------
# grouped bar plot ami
ami_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  ami_col <- metric_df[['AMI']]
  #print(ami_col)
  
  ami_df[1:length(ami_col),method] <- ami_col
}
# skip Covid650K temporarily
ami_df <- ami_df[1:10,]
# prepare for melt
ami_df$dataset_ID <- rownames(ami_df)
ami_df_long <- melt(ami_df,id.vars='dataset_ID')
colnames(ami_df_long) <- c('dataset_ID','method','ami')
ami_df_long$dataset_ID <- factor(ami_df_long$dataset_ID, levels = unique(ami_df_long$dataset_ID)) 
ami_df_long$method <- factor(ami_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))

# plot ami_long
ami_df_long[is.na(ami_df_long)] <-0 
ami_plot <- ggplot(ami_df_long,aes(x= dataset_ID,y= ami,fill= method)) + geom_bar(stat = "identity", position = "dodge",width=0.8) +  
  geom_text(aes(label = ifelse(ami==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "AMI", fill = "")  
print(ami_plot)
ggsave(paste0(output_dir,"scMINER_clustering_AMI.pdf"),plot=ami_plot)

# mean + variance errorbar + Wilcoxon test
#ami_df_long$ami <- ami_df_long$ami + runif(nrow(ami_df_long), min = -1e-10, max = 1e-10)  # break tied value

ami_summary_plot <- ggplot(ami_df_long, aes(x = method,y=ami,fill=method)) +
  
  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +  
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.9)) +
  
  geom_signif(
    comparisons = list(c("scMINER","Seurat"),c("scMINER","SC3s"),
                       c("scMINER","Scanpy"),c("scMINER","scVI"),
                       c("scMINER","scDeepCluster")),
    test = "wilcox.test",test.args = list(paired = TRUE, alternative = "two.sided"),
    #map_signif_level = TRUE,y_position = c(0.9,0.95,1,1.05,1.1)
    map_signif_level = TRUE,y_position = c(1,1.15,1.3,1.45,1.6),textsize = 5
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
  labs(x = "", y = "AMI", fill = "")

print(ami_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_AMI_wilcoxon.pdf"),plot=ami_summary_plot)

# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_ami_Wilcoxon <- data.frame(  
  method = character(),   
  statistic = numeric(),   
  p.value = numeric(),   
  stringsAsFactors=FALSE  
) 

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(ami_df_long[ami_df_long$method == 'scMINER','ami'],
                             ami_df_long[ami_df_long$method == method,'ami'],
                             alternative="two.sided",paired=TRUE)
  scMINER_ami_Wilcoxon <- rbind(  
    scMINER_ami_Wilcoxon,   
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)  
  ) 
}

write.table(scMINER_ami_Wilcoxon,file=paste0(output_dir,"scMINER_AMI_Wilcoxon.txt"),sep="\t",row.names = FALSE)


#----------

# NMI
#-------------
# grouped bar plot ami
nmi_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  nmi_col <- metric_df[['NMI']]
  #print(nmi_col)
  
  nmi_df[1:length(nmi_col),method] <- nmi_col
}
# skip Covid650K temporarily
nmi_df <- nmi_df[1:10,]
# prepare for melt
nmi_df$dataset_ID <- rownames(nmi_df)
nmi_df_long <- melt(nmi_df,id.vars='dataset_ID')
colnames(nmi_df_long) <- c('dataset_ID','method','nmi')
nmi_df_long$dataset_ID <- factor(nmi_df_long$dataset_ID, levels = unique(nmi_df_long$dataset_ID)) 
nmi_df_long$method <- factor(nmi_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))


#plot nmi
nmi_df_long[is.na(nmi_df_long)] <-0 
nmi_plot <- ggplot(nmi_df_long,aes(x= dataset_ID,y= nmi,fill= method)) + geom_bar(stat = "identity", position = "dodge",width=0.8) +  
  geom_text(aes(label = ifelse(nmi==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "NMI", fill = "")  
print(nmi_plot)
ggsave(paste0(output_dir,"scMINER_clustering_NMI.pdf"),plot=nmi_plot)

# mean + variance errorbar + Wilcoxon test
#nmi_df_long$nmi <- nmi_df_long$nmi + runif(nrow(nmi_df_long), min = -1e-10, max = 1e-10)  # break tied value

nmi_summary_plot <- ggplot(nmi_df_long, aes(x = method,y=nmi,fill=method)) +
  
  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +  
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.9)) +
  
  geom_signif(
    comparisons = list(c("scMINER","Seurat"),c("scMINER","SC3s"),
                       c("scMINER","Scanpy"),c("scMINER","scVI"),
                       c("scMINER","scDeepCluster")),
    test = "wilcox.test",test.args = list(paired = TRUE, alternative = "two.sided"),
    #map_signif_level = TRUE,y_position = c(0.9,0.95,1,1.05,1.1)
    map_signif_level = TRUE,y_position = c(1,1.1,1.2,1.3,1.4),textsize = 5,
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
  labs(x = "", y = "NMI", fill = "")

print(nmi_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_NMI_wilcoxon.pdf"),plot=nmi_summary_plot)

# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_nmi_Wilcoxon <- data.frame(  
  method = character(),   
  statistic = numeric(),   
  p.value = numeric(),   
  stringsAsFactors=FALSE  
) 

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(nmi_df_long[nmi_df_long$method == 'scMINER','nmi'],
                             nmi_df_long[nmi_df_long$method == method,'nmi'],
                             alternative="two.sided",paired=TRUE)
  scMINER_nmi_Wilcoxon <- rbind(  
    scMINER_nmi_Wilcoxon,   
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)  
  ) 
}

write.table(scMINER_nmi_Wilcoxon,file=paste0(output_dir,"scMINER_NMI_Wilcoxon.txt"),sep="\t",row.names = FALSE)

#----------

#---------
# ASW scaled to [0,1] skipped

# -----------
ASW_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  asw_col <- metric_df[['ASW']]
  print(asw_col)

  ASW_df[1:length(asw_col),method] <- asw_col
}

# skip Covid650K temporarily
ASW_df <- ASW_df[1:10,]
# prepare for melt
ASW_df$dataset_ID <- rownames(ASW_df)
ASW_df_long <- melt(ASW_df,id.vars='dataset_ID')
colnames(ASW_df_long) <- c('dataset_ID','method','ASW')
ASW_df_long$dataset_ID <- factor(ASW_df_long$dataset_ID, levels = unique(ASW_df_long$dataset_ID))
ASW_df_long$method <- factor(ASW_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))



# grouped barplot ari
ASW_df_long[is.na(ASW_df_long)] <-0
asw_plot <- ggplot(ASW_df_long,aes(x= dataset_ID,y= ASW,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) +
  geom_text(aes(label = ifelse(ASW==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) +
  theme(text = element_text(size = 14),aspect.ratio = 3/5,
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

  labs(x = "", y = "ASW", fill = "")
print(asw_plot)
ggsave(paste0(output_dir,"scMINER_clustering_ASW.pdf"),plot=asw_plot)

#ASW_df_long$ari <- ASW_df_long$ari + runif(nrow(ASW_df_long), min = -1e-10, max = 1e-10)  # break tied value

asw_summary_plot <- ggplot(ASW_df_long, aes(x = method,y=ASW,fill=method)) +

  stat_summary(fun = mean, geom = "bar", position = position_dodge(),width=0.5, na.rm = TRUE) +
  stat_summary(fun.data = mean_SD, geom = "errorbar", width = 0.1, position = position_dodge(.9)) +

  geom_signif(
    comparisons = list(c("scMINER","Seurat"),c("scMINER","SC3s"),
                       c("scMINER","Scanpy"),c("scMINER","scVI"),
                       c("scMINER","scDeepCluster")),
    test = "wilcox.test",test.args = list(paired = TRUE, alternative = "two.sided"),
    map_signif_level = TRUE,y_position = c(.9,1,1.1,1.2,1.3),textsize = 5
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
        axis.line = element_line(color = "black",size=1),
        axis.ticks.length = unit(-0.1,"cm")) +

  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = seq(0, 1, by = 0.2)) +

  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) + theme(legend.position ="none")+
  labs(x = "", y = "ASW", fill = "")
print(asw_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_ASW_wilcoxon.pdf"),plot=asw_summary_plot)


# sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_asw_Wilcoxon <- data.frame(
  method = character(),
  statistic = numeric(),
  p.value = numeric(),
  stringsAsFactors=FALSE
)

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(ASW_df_long[ASW_df_long$method == 'scMINER','ASW'],
                             ASW_df_long[ASW_df_long$method == method,'ASW'],
                             alternative="two.sided",paired=TRUE)
  scMINER_asw_Wilcoxon <- rbind(
    scMINER_ari_Wilcoxon,
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)
  )
}

write.table(scMINER_asw_Wilcoxon,file=paste0(output_dir,"scMINER_ASW_Wilcoxon.txt"),sep="\t",row.names = FALSE)

#--------------------
# AvgBIO

# -----------
AvgBIO_df <- data.frame(row.names = ordered_dataset_ids)
for (method in packages){
  metric_file <- list.files(path=input_dir,pattern = paste0(method,'.*.csv'))
  #print(metric_file)
  metric_file <- paste0(input_dir,metric_file)
  metric_df <- read.table(metric_file,sep = '\t',header = TRUE)
  #print(metric_df)
  ari_col <- metric_df[['AvgBIO']]
  print(ari_col)
  
  AvgBIO_df[1:length(ari_col),method] <- ari_col
}

# skip Covid650K temporarily
AvgBIO_df <- AvgBIO_df[1:10,]
# prepare for melt
AvgBIO_df$dataset_ID <- rownames(AvgBIO_df)
AvgBIO_df_long <- melt(AvgBIO_df,id.vars='dataset_ID')
colnames(AvgBIO_df_long) <- c('dataset_ID','method','AvgBIO')
AvgBIO_df_long$dataset_ID <- factor(AvgBIO_df_long$dataset_ID, levels = unique(AvgBIO_df_long$dataset_ID)) 
AvgBIO_df_long$method <- factor(AvgBIO_df_long$method,levels =c('scMINER','Seurat','SC3s','Scanpy','scVI','scDeepCluster'))


# grouped barplot ari
AvgBIO_df_long[is.na(AvgBIO_df_long)] <-0 
ari_plot <- ggplot(AvgBIO_df_long,aes(x= dataset_ID,y= AvgBIO,fill= method)) + geom_bar(stat = "identity", position = 'dodge',width=.8) + 
  geom_text(aes(label = ifelse(AvgBIO==0, "*", "")), na.rm = TRUE,position=position_dodge(width=0.7), vjust = 0.5,hjust=0.5) + 
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
  
  labs(x = "", y = "AvgBIO", fill = "")  
print(ari_plot)
ggsave(paste0(output_dir,"scMINER_clustering_AvgBIO.pdf"),plot=ari_plot)

#AvgBIO_df_long$ari <- AvgBIO_df_long$ari + runif(nrow(AvgBIO_df_long), min = -1e-10, max = 1e-10)  # break tied value

ari_summary_plot <- ggplot(AvgBIO_df_long, aes(x = method,y=AvgBIO,fill=method)) +
  
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
        axis.line = element_line(color = "black",size=1),
        axis.ticks.length = unit(-0.1,"cm")) +   
  
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),  
                     labels = seq(0, 1, by = 0.2)) +  
  
  scale_fill_manual(values = c("scMINER" = '#ff0000', "Seurat" = "#AED581","SC3s" = "#AED6F1","Scanpy" = "#8b5635",
                               "scVI" = "#fd9c32","scDeepCluster" = "#E1BEE7")) + theme(legend.position ="none")+
  labs(x = "", y = "AvgBIO", fill = "")
print(ari_summary_plot)
ggsave(paste0(output_dir,"scMINER_clustering_AvgBIO_wilcoxon.pdf"),plot=ari_summary_plot)


#sanity check of paired-Wilcoxon one-sided signed rank test
scMINER_ari_Wilcoxon <- data.frame(
  method = character(),
  statistic = numeric(),
  p.value = numeric(),
  stringsAsFactors=FALSE
)

for (method in c('Seurat','SC3s','Scanpy','scVI','scDeepCluster')){
  test_result <- wilcox.test(AvgBIO_df_long[AvgBIO_df_long$method == 'scMINER','AvgBIO'],
                                  AvgBIO_df_long[AvgBIO_df_long$method == method,'AvgBIO'],
                                  alternative="two.sided",paired=TRUE)
  scMINER_ari_Wilcoxon <- rbind(
    scMINER_ari_Wilcoxon,
    data.frame(method = method, statistic = test_result$statistic, p.value = test_result$p.value)
  )
}

write.table(scMINER_ari_Wilcoxon,file=paste0(output_dir,"scMINER_AvgBIO_Wilcoxon.txt"),sep="\t",row.names = FALSE)
#---------------


# Accuracy  
# For each true cell type, % of major prediction
get_cluster_accuracy <- function(true_label,pred_label){
  
  if(length(true_label) != length(pred_label)){
    cat('Error, true label and pred label not equal length \n')
    return(NA)
  }
  true_label = as.numeric(factor(true_label))
  pred_label = as.numeric(factor(pred_label))
  
  accuracy <- unlist(lapply(unique(pred_label),function(x){cal_entropy(true_label[pred_label==x])}))
  return(accuracy)
}


# Purity
# For each predicted cluster, % of true cell type
get_cluster_purity <- function(true_label,pred_label){

  if(length(true_label) != length(pred_label)){
    cat('Error, true label and pred label not equal length \n')
    return(NA)
  }
  
  true_label = as.numeric(factor(true_label))
  pred_label = as.numeric(factor(pred_label))
  purity <- unlist(lapply(unique(true_label),function(x){cal_entropy(pred_label[true_label==x])}))
  return(purity)
}






