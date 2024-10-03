rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(SingleCellExperiment)
library(Seurat)
library(scMINER)
library(Matrix)
library(dplyr)
library(anndata)
library(aricode)

filt_eset.log2 <- readRDS('log2cpm_filtered.8846_13605.rds')
cat(dim(exprs(filt_eset.log2))) # 8846 13605

# convert true label in eSet
unique(filt_eset.log2$true_label)
old_true_label <- c("CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T","CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic",
                    "CD14+ Monocyte","CD19+ B","CD56+ NK")
new_true_label <- c("CD4+ Treg","CD4+ Naive T","CD4+ TCM","CD8+Naive CTL","Monocyte","B cell","NK")
filt_eset.log2$true_label <- new_true_label[match(filt_eset.log2$true_label, old_true_label)]
unique(filt_eset.log2$true_label)

# set up clustering output folder
clustering_root <- './'
HVG_sweep <- c('HVG1000','HVG2000','HVG3000')
dataset <- 'PBMC14K'
packages = c('Seurat','SC3s','Scanpy','scVI','scDeepCluster') # scMINER uses all genes instead of picking HVG

#--------------
# Umap for different packages
ARI_df <- data.frame()

for (method in packages){
  # compare scMINER with true label with other package
  for (HVG in HVG_sweep){
    clustering_dir <- paste0(clustering_root,method,'/',HVG)
    cluster_file <- list.files(path=clustering_dir,
                               pattern=paste0(method,'.*.txt'),
                               full.names=TRUE)
    cluster_df <- read.table(cluster_file,header=TRUE,sep='\t',row.name=1)
    
    #break
    
    UMAP_X_key <- paste0('UMAP_X_',method)
    UMAP_Y_key <- paste0('UMAP_Y_',method)
    if ('X' %in% colnames(cluster_df)){
      X_cor <- cluster_df$X
      Y_cor <- cluster_df$Y
    } else if  ('UMAP_X' %in% colnames(cluster_df)) {
      X_cor <- cluster_df$UMAP_X
      Y_cor <- cluster_df$UMAP_Y    
    } else{
      cat('Not Umap coordinate column found.')
    }
    
    
    filt_eset.log2[[UMAP_X_key]] <- X_cor
    filt_eset.log2[[UMAP_Y_key]] <- Y_cor
    # add clustering label
    cluster_key <- paste0(method,'_label')
    if (method == 'scMINER'){
      filt_eset.log2[[cluster_key]] <- factor(cluster_df[['label']]) 
    } else if (method == 'Seurat'){
      filt_eset.log2[[cluster_key]] <- cluster_df[[cluster_key]]
      filt_eset.log2[[cluster_key]] <- factor(filt_eset.log2[[cluster_key]])  
    } else {
      filt_eset.log2[[cluster_key]] <- cluster_df[[cluster_key]] + 1
      filt_eset.log2[[cluster_key]] <- factor(filt_eset.log2[[cluster_key]])
    }
    # break  
    
    # ARI 
    ARI_df[method,HVG] <- round(ARI(filt_eset.log2$true_label,filt_eset.log2[[cluster_key]]),2)
    
    #next    
    
    # Umap plotting
    #--------------  
    MICAplot(input_eset = filt_eset.log2, X=UMAP_X_key,Y=UMAP_Y_key,show_label = FALSE,
             color_by = cluster_key, pct=0.5,label.size = 8) + 
      theme(aspect.ratio = 1,
            legend.text = element_text(size=10),
            text = element_text(size = 10),
            plot.margin = margin(1,1,1,1,"cm"),
            axis.text = element_text(size = 20),
            axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
            panel.background = element_blank(),  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            #panel.border = element_blank(),  
            panel.border = element_rect(color = "black", fill = NA, linewidth  = 1),  
            #axis.line = element_line(linewidth = 1),
            axis.line = element_blank(),
            axis.ticks.length = unit(-0.2,"cm")) 
    Umap_file <- paste0(clustering_dir,'/',dataset,'_',method,'_',HVG,'_Umap.pdf')
    ggsave(Umap_file)
    # # Treg projection
    # MICAplot(input_eset = filt_eset.log2, X=UMAP_X_key,Y=UMAP_Y_key,
    #          colors = c('#dfdfdf','#dfdfdf','#dfdfdf','#f8766d','#dfdfdf','#dfdfdf','#dfdfdf'),color_by = 'true_label', pct=1) + 
    #   theme(aspect.ratio = 1,
    #         text = element_text(size = 10),
    #         plot.margin = margin(1,1,1,1,"cm"),
    #         axis.text = element_text(size = 24),
    #         axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    #         axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
    #         panel.background = element_blank(),  
    #         panel.grid.major = element_blank(),  
    #         panel.grid.minor = element_blank(),  
    #         #panel.border = element_blank(),  
    #         panel.border = element_rect(color = "black", fill = NA, size = 1), 
    #         #axis.line = element_line(linewidth = 1),
    #         axis.line = element_blank(),
    #         axis.ticks.length = unit(-0.2,"cm")) 
    # Umap_file <- paste0(clustering_dir,'/',dataset,'_',method,'_',HVG,'_Treg_true_Umap.pdf')
    # ggsave(Umap_file)
    # Treg and TCM projection
    MICAplot(input_eset = filt_eset.log2, X=UMAP_X_key,Y=UMAP_Y_key,
             colors = c('#dfdfdf','#dfdfdf','#53b400','#f8766d','#dfdfdf','#dfdfdf','#dfdfdf'),color_by = 'true_label', pct=0.5,label.size = 8) + 
      theme(aspect.ratio = 1,
            legend.text = element_text(size=10),
            text = element_text(size = 10),
            plot.margin = margin(1,1,1,1,"cm"),
            axis.text = element_text(size = 20),
            axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
            panel.background = element_blank(),  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            #panel.border = element_blank(),  
            panel.border = element_rect(color = "black", fill = NA, size = 1), 
            #axis.line = element_line(linewidth = 1),
            axis.line = element_blank(),
            axis.ticks.length = unit(-0.2,"cm")) 
    Umap_file <- paste0(clustering_dir,'/',dataset,'_',method,'_',HVG,'_Treg_TCM_true_Umap.pdf')
    ggsave(Umap_file)
    #--------------
  }
}
#------------

# save ARI
write.table(ARI_df,'PBMC14K_HVG_ARI.txt',sep='\t',quote=FALSE)

