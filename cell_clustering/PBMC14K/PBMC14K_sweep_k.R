rm(list=ls())
current_script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_script_dir)

library(SingleCellExperiment)
library(Seurat)
library(scMINER)
library(Matrix)
library(dplyr)
library(anndata)
# cluster tree
library(clustree)
# load Sparse eset 
eset.log2 <- readRDS('log2cpm_filtered.8846_13605.rds')
cat(dim(exprs(eset.log2))) # 8846 13605
# convert true label in eSet
unique(eset.log2$true_label)
old_true_label <- c("CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T","CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic",
                    "CD14+ Monocyte","CD19+ B","CD56+ NK")
new_true_label <- c("CD4+ Treg","CD4+ Naive T","CD4+ TCM","CD8+Naive CTL","Monocyte","B cell","NK")
eset.log2$true_label <- new_true_label[match(eset.log2$true_label, old_true_label)]
unique(eset.log2$true_label)


# set up path
packages <- c('scMINER','Seurat','Scanpy')
dataset <- 'PBMC14K'

# add scMINER label and Seurat label to Sparse eSet
#------------
for (method in packages){
  cluster_files <- list.files(path=method,
                             pattern='*.txt',
                             full.names=TRUE)
 
  for (each_file in cluster_files){
    cluster_df <- read.table(each_file,header=TRUE,sep='\t',row.name=1)

    # check pred_k
    if (method == 'scMINER'){
      # rename df column
      colnames(cluster_df) <- c('UMAP_X','UMAP_Y','scMINER_label')
    }
    old_cluster_key <- paste0(method,'_label')
    pred_k <- length(unique(cluster_df[[old_cluster_key]]))
    # set new cluster label
    cluster_key <- paste0(old_cluster_key,'_k',pred_k)
    print(cluster_key)
  
    # add cluster label
    eset.log2[[cluster_key]] <- cluster_df[[old_cluster_key]]
    eset.log2[[cluster_key]] <- factor(eset.log2[[cluster_key]]) 
    
    
    if (method == 'Scanpy'){
      eset.log2[[cluster_key]] <- cluster_df[[old_cluster_key]] + 1
      eset.log2[[cluster_key]] <- factor(eset.log2[[cluster_key]]) 
    }
    
    # Umap coord
    X_cor <- cluster_df$UMAP_X
    Y_cor <- cluster_df$UMAP_Y    

    # set new Umap coord key
    UMAP_X_key <- paste0('UMAP_X_',method,'_k',pred_k)
    UMAP_Y_key <- paste0('UMAP_Y_',method,'_k',pred_k)

    # add Umap coord to eSet
    eset.log2[[UMAP_X_key]] <- X_cor
    eset.log2[[UMAP_Y_key]] <- Y_cor
    #break
  #skip plot
  #next
  # Umap 
  MICAplot(input_eset = eset.log2, X=UMAP_X_key,Y=UMAP_Y_key,
           color_by = cluster_key) + 
    theme(aspect.ratio = 1,
          text = element_text(size = 10),
          plot.margin = margin(1,1,1,1,"cm"),
          axis.text = element_text(size = 24),
          axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
          panel.background = element_blank(),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          #panel.border = element_blank(),  
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
          #axis.line = element_line(linewidth = 1),
          axis.line = element_blank(),
          axis.ticks.length = unit(-0.2,"cm")) 
  
  #break
  Umap_file <- paste0(method,'/',dataset,'_',method,'_k',pred_k,'_Umap.pdf')
  ggsave(Umap_file)
  
  
  }
  #break
}

#-------

# T subtype true label projection on Seurat clustering embedding
MICAplot(input_eset = eset.log2, X='UMAP_X_Seurat_k7',Y='UMAP_Y_Seurat_k7',show.cluster_label  = FALSE,
         colors = c('#dfdfdf','#39caa5','#f8766d','#53b400','#e79d38','#dfdfdf','#dfdfdf'),
         color_by = 'true_label') + 
  theme(aspect.ratio = 1,
        text = element_text(size = 10),
        plot.margin = margin(1,1,1,1,"cm"),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        #panel.border = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
        #axis.line = element_line(linewidth = 1),
        axis.line = element_blank(),
        axis.ticks.length = unit(-0.2,"cm")) 

ggsave('PBMC14K_Seurat_T_cell_true_label_projection.pdf')
# scMINER
# T subtype true label projection on Seurat clustering embedding
MICAplot(input_eset = eset.log2, X='UMAP_X_scMINER_k7',Y='UMAP_Y_scMINER_k7',show.cluster_label  = FALSE,
         colors = c('#dfdfdf','#39caa5','#f8766d','#53b400','#e79d38','#dfdfdf','#dfdfdf'),
         color_by = 'true_label') + 
  theme(aspect.ratio = 1,
        text = element_text(size = 10),
        plot.margin = margin(1,1,1,1,"cm"),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        #panel.border = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
        #axis.line = element_line(linewidth = 1),
        axis.line = element_blank(),
        axis.ticks.length = unit(-0.2,"cm")) 

ggsave('PBMC14K_scMINER_T_cell_true_label_projection.pdf')
######
# Scanpy
# T subtype true label projection on Seurat clustering embedding
MICAplot(input_eset = eset.log2, X='UMAP_X_Scanpy_k7',Y='UMAP_Y_Scanpy_k7',show.cluster_label = FALSE,
         colors = c('#dfdfdf','#39caa5','#f8766d','#53b400','#e79d38','#dfdfdf','#dfdfdf'),
         color_by = 'true_label') + 
  theme(aspect.ratio = 1,
        text = element_text(size = 10),
        plot.margin = margin(1,1,1,1,"cm"),
        axis.text = element_text(size = 24),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),           
        panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        #panel.border = element_blank(),  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
        #axis.line = element_line(linewidth = 1),
        axis.line = element_blank(),
        axis.ticks.length = unit(-0.2,"cm")) 

ggsave('PBMC14K_Scanpy_T_cell_true_label_projection.pdf')
#------


#

#clustree 

#-------
# Seurat
cell_meta_df <- eset.log2@phenoData@data
clustree(cell_meta_df,
         prefix = 'Seurat_label_k') 
ggsave('PBMC14K_Seurat_clustree.pdf')

# Scanpy
cell_meta_df <- eset.log2@phenoData@data
clustree(cell_meta_df,
         prefix = 'Scanpy_label_k') 
ggsave('PBMC14K_Scanpy_clustree.pdf')


# scMINER
clustree(cell_meta_df,
         prefix = 'scMINER_label_k')
ggsave('PBMC14K_scMINER_clustree.pdf')
# overlay clustree
clustree_overlay(cell_meta_df,
                 prefix = 'Seurat_label_k',
                 x_value='UMAP_X_Seurat_k7',
                 y_value='UMAP_Y_Seurat_k7')

#---------



