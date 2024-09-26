library(Seurat)
library(scMINER)

dir <- "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao"

rna.3p <- Read10X(data.dir = paste0(dir, "/00_rawData/GSM5008737_RNA_3P"))
adt.3p <- Read10X(data.dir = paste0(dir, "/00_rawData/GSM5008738_ADT_3P"))

meta <- read.table(paste0(dir, "/00_rawData/GSE164378_sc.meta.data_3P.csv"), header = T, sep = ",", quote = "", stringsAsFactors = F)
row.names(meta) <- meta$X..
colnames(meta)[1] <- "cellid"
meta$celltype <- gsub(" ", "_", meta$celltype.l1); table(meta$celltype)

## Protein
rawCount.eset <- CreateSparseEset(data = as.matrix(adt.3p), meta.data = meta[colnames(adt.3p),], add.meta = F)
saveRDS(rawCount.eset, file = paste0(dir, "/01_sparseEset.ADT_rawCounts.rds"))

# normalization
norm = 1e6
exp.norm <- sweep(exprs(rawCount.eset), 2, norm/unname(Matrix::colSums(exprs(rawCount.eset))), '*')
exp.log2 <- log(exp.norm + 1, base = 2)

pd <- pData(rawCount.eset)[colnames(exp.log2),]

log2CPM.eset <- CreateSparseEset(data = exp.log2, meta.data = pd, add.meta = F)
saveRDS(log2CPM.eset, file = paste0(dir, "/01_sparseEset.ADT_log2CPM.rds"))


## RNA
rawCount.eset <- CreateSparseEset(data = as.matrix(rna.3p), meta.data = meta[colnames(rna.3p),], add.meta = F)
saveRDS(rawCount.eset, file = paste0(dir, "/01_sparseEset.RNA_rawCounts.rds"))

# normalization
norm = 1e6
exp.norm <- sweep(exprs(rawCount.eset), 2, norm/unname(Matrix::colSums(exprs(rawCount.eset))), '*')
exp.log2 <- log(exp.norm + 1, base = 2)

pd <- pData(rawCount.eset)[colnames(exp.log2),]

log2CPM.eset <- CreateSparseEset(data = exp.log2, meta.data = pd, add.meta = F)
saveRDS(log2CPM.eset, file = paste0(dir, "/01_sparseEset.RNA_log2CPM.rds"))


## generate SJARACNe inputs
log2CPM.eset$celltype <- factor(pData(log2CPM.eset)$celltype)
generateSJARACNeInput(input_eset = log2CPM.eset, funcType = "TF", ref = "hg", wd.src = paste0(dir, "/02_sjaracneAll"), group_name = "celltype")
generateSJARACNeInput(input_eset = log2CPM.eset, funcType = "SIG", ref = "hg", wd.src = paste0(dir, "/02_sjaracneAll"), group_name = "celltype")


## down-sample 1000
downsample_matrix <- function(mat, n = 1, keep_cluster_proportions = TRUE, metadata = NULL, cluster_col = "cluster") {
  cluster_info <- metadata
  if (keep_cluster_proportions == FALSE) {
    cluster_ids <- colnames(mat)
    if (n < 1) {
      n <- as.integer(ncol(mat) * n)
    }
    cluster_ids_new <- sample(cluster_ids, n)
  } else {
    if (is.vector(cluster_info)) {
      cluster_ids <- split(colnames(mat), cluster_info)
    } else if (is.data.frame(cluster_info) &
               !is.null(cluster_col)) {
      cluster_ids <- split(colnames(mat), cluster_info[[cluster_col]])
    } else if (is.factor(cluster_info)) {
      cluster_info <- as.character(cluster_info)
      cluster_ids <- split(colnames(mat), cluster_info)
    } else {
      stop("metadata not formatted correctly,
         supply either a  vector or a dataframe",
           call. = FALSE
      )
    }
    if (n < 1) {
      n2 <- vapply(cluster_ids, function(x) {
        as.integer(length(x) * n)
      }, FUN.VALUE = numeric(1))
      n <- n2
    }
    cluster_ids_new <-
      mapply(sample, cluster_ids, n, SIMPLIFY = FALSE)
  }
  return(mat[, unlist(cluster_ids_new)])
}

set.seed(111)
log2CPM.sel <- downsample_matrix(mat = exprs(log2CPM.eset),
                                 metadata = pData(log2CPM.eset),
                                 n = 1000,
                                 keep_cluster_proportions = T,
                                 cluster_col = "celltype")

log2CPM.eset.downsample <- log2CPM.eset[, colnames(log2CPM.sel)]

generateSJARACNeInput(input_eset = log2CPM.eset.downsample, funcType = "TF", ref = "hg", wd.src = paste0(dir, "/02_sjaracneDownSample1000"), group_name = "celltype")
generateSJARACNeInput(input_eset = log2CPM.eset.downsample, funcType = "SIG", ref = "hg", wd.src = paste0(dir, "/02_sjaracneDownSample1000"), group_name = "celltype")

