library(SingleCellExperiment)
library(Seurat)

library(Matrix)
library(ggplot2)

library(reticulate)
use_condaenv('scVI', required = TRUE)
library(argparse)
library(anndata)
library(sceasy)

# for evaluation
library(cluster)
library(aricode)
library(scales)


rm(list=ls())

# add command line argument
parser <- ArgumentParser()

parser$add_argument("-i", "--input_file", default=NULL,
    help=".h5ad format input file path")
parser$add_argument("-o", "--save_dir", default=NULL, 
    help="output folder")
parser$add_argument("-t", "--input_type", default=NULL, 
    help="type of input count matrix, can be raw,tpm,log1p")

# check output folder
args <- parser$parse_args()
if (!file.exists(args$save_dir)) {  
  dir.create(args$save_dir)  
} 
args$save_dir <- sub("/$","",args$save_dir)


# extract file name and dataset_id
#Seurat_file <- sub(".h5ad", ".rds", basename(args$input_file))
h5ad_file <- basename(args$input_file)
dataset_id <-  strsplit(h5ad_file,'_')[[1]][1]

#----------
# # sceasy convert h5ad to seurat
# seuratRDS <- paste0(args$save_dir,'/',dataset_id,'_Seurat.rds')
# sceasy::convertFormat(args$input_file, from="anndata", to="seurat",
#                       outFile=seuratRDS)
# 
# seuratObj <- readRDS(seuratRDS)
#------------


# function to convert adata to Seurat object
#----------------
adata2seurat <- function(adata){
  exp.mat <- t(adata$X)
  #exp.mat <- as.csc.matrix(exp.mat)
  
  meta.data <- adata$obs
  rownames(meta.data) <- colnames(exp.mat)
  seuratObj <- CreateSeuratObject(
    exp.mat,
    meta.data = adata$obs
  )
  return(seuratObj)
}
# read h5ad and convert
adata <- read_h5ad(args$input_file)
seuratObj <- adata2seurat(adata)

# normalize 
if (args$input_type == 'raw'){
  seuratObj <- NormalizeData(seuratObj) # default is normalization.method =  "LogNormalize", scale.factor = 10000
} else if (args$input_type == 'tpm'){
  print("TPM transformed")
  seuratObj <- SetAssayData(object=seuratObj, layer ="data", new.data = log1p(seuratObj@assays$RNA$counts))
} else if (args$input_type == 'log1p'){
  print("Input already log1p transformed.")
  seuratObj <- SetAssayData(object=seuratObj, layer ="data", new.data = seuratObj@assays$RNA$counts) 
} else {
  print("Error, wrong type of input")
  q("no")
}
# save cleaned SeuratObj to .Rdata file
seuratRdata_file <- paste0(args$save_dir,'/',dataset_id,'_Seurat.Rdata')
save(seuratObj,file=seuratRdata_file)

#----------------
# # function to convert adata to Seurat object
# adata2seurat <- function(adata){
#   exp.mat <- t(adata$X)
#   #exp.mat <- as.csc.matrix(exp.mat)
# 
#   meta.data <- adata$obs
#   rownames(meta.data) <- colnames(exp.mat)
#   seuratObj <- CreateSeuratObject(
#     exp.mat,
#     meta.data = adata$obs
#   )
#   return(seuratObj)
# }
# # read h5ad and convert
# adata <- read_h5ad(args$input_file)
# seuratObj <- adata2seurat(adata)
# 
# # normalize 
# if (args$input_type == 'raw'){
#     seuratObj <- NormalizeData(seuratObj) # default is normalization.method =  "LogNormalize", scale.factor = 10000
# } else if (args$input_type == 'tpm'){
#     print("TPM transformed")
#     seuratObj <- SetAssayData(object=seuratObj, layer ="data", new.data = log1p(seuratObj@assays$RNA$counts))
# } else if (args$input_type == 'log1p'){
#     print("Input already log1p transformed.")
#     seuratObj <- SetAssayData(object=seuratObj, layer ="data", new.data = seuratObj@assays$RNA$counts) 
# } else {
#     print("Error, wrong type of input")
#     q("no")
# }
# # save cleaned SeuratObj to .rds file
# save(seuratObj,file=paste0(args$save_dir,'/',dataset_id,'_Seurat.rds'))


#-------------
# select HVG from layer data
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000) # consistent with scVI, scDeepCluster

# Scale and PCA
seuratObj <- ScaleData(seuratObj) 
seuratObj <- RunPCA(seuratObj,features=VariableFeatures(object = seuratObj)) # PCA has to run ScaleData

# save ReducedDim to txt file
pca_matrix <- Embeddings(seuratObj, "pca")
write.table(pca_matrix,file=paste0(args$save_dir,'/pca_matrix.txt'),sep='\t',row.names=TRUE)
#pca_40 <- pca_matrix[,1:40]
#write.table(pca_40,file=paste0(args$save_dir,'/pca_40.txt'),sep='\t',row.names=TRUE)

# Umap
seuratObj <- RunUMAP(seuratObj,dims = 1:10)
umap_coordinates <- Embeddings(seuratObj, "umap")
# UMAP_X <- umap_coordinates [,1]
# UMAP_Y <- umap_coordinates [,2]

# test exit
#q("no")

# clustering and evaluation
seuratObj <- FindNeighbors(seuratObj,reduction="pca") 
seuratObj <- FindClusters(seuratObj) # default is Louvain, default resolution=0.8


# Silhouette analysis
# extract Seurat pred label
seurat_label <- as.integer(seuratObj$seurat_clusters)
n_cluster <- length(unique(seurat_label))

# plot Umap
Umap <- DimPlot(seuratObj, reduction = "umap")
Umap_file <- paste0(args$save_dir,'/Seurat_Umap_k',n_cluster,'.pdf')
ggsave(Umap_file,plot=Umap)

# wrtie clustering result to file 
df <- data.frame(
    cell_id = seuratObj$cell_id,
    UMAP_X = umap_coordinates [,1],
    UMAP_Y = umap_coordinates [,2],
    Seurat_label = seurat_label
)
label_file <- paste0(args$save_dir,'/Seurat_clustering_k',n_cluster,'.txt')
write.table(df,file=label_file,sep='\t',row.names=TRUE,quote = FALSE)
#break

# save clustered
clustered_file <- paste0(args$save_dir,'/','Seurat_clustered.Rdata')
save(seuratObj, file = clustered_file)








