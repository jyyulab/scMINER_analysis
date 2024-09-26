library(Seurat)
library(scMINER)

ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}
Read10X_new <- function (data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- Matrix::readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                sep = "\t", row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                              yes = gene.loc, no = features.loc), header = FALSE, 
                                stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "CsparseMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}

obj <- Read10X_new(data.dir = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/00_rawData/outs", strip.suffix = T)
rna <- CreateSeuratObject(counts = obj$`Gene Expression`, project = "PBMC", min.cells = 0, min.features = 0)
adt <- CreateSeuratObject(counts = obj$`Antibody Capture`, project = "PBMC", min.cells = 0, min.features = 0)

rna.raw <- rna@assays$RNA@counts; dim(rna.raw)
adt.raw <- adt@assays$RNA@counts; dim(adt.raw)

pd <- read.csv(file = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/00_rawData/GSE213282_metadata.csv", sep = ",", header = TRUE); head(pd)
colnames(pd)[1] <- "cellID"; head(pd) 
row.names(pd) <- pd$cellID
row.names(pd) <- gsub("\\.1","", row.names(pd))
row.names(pd) <- gsub("HTO\\d_","", row.names(pd))
length(row.names(pd))

table(colnames(rna.raw) %in% row.names(pd))
rna.sel <- as.matrix(rna.raw[,colnames(rna.raw) %in% row.names(pd)]); dim(rna.sel)
rna.sel <- rna.sel[rowSums(rna.sel) > 0,]; dim(rna.sel)

table(colnames(adt.raw) %in% row.names(pd))
adt.sel <- as.matrix(adt.raw[,colnames(adt.raw) %in% row.names(pd)]); dim(adt.sel)
adt.sel <- adt.sel[rowSums(adt.sel) > 0,]; dim(adt.sel)

###### RNA outputs
rawCount.eset <- CreateSparseEset(data = rna.sel, meta.data = pd[colnames(rna.sel),], add.meta = F)
saveRDS(rawCount.eset, file = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.RNA_rawCounts.rds")

# normalization
norm = 1e6
exp.norm <- sweep(exprs(rawCount.eset), 2, norm/unname(Matrix::colSums(exprs(rawCount.eset))), '*')
exp.log2 <- log(exp.norm + 1, base = 2)

pd <- pData(rawCount.eset)[colnames(exp.log2),]

log2CPM.eset <- CreateSparseEset(data = exp.log2, meta.data = pd, add.meta = F)
saveRDS(log2CPM.eset, file = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.RNA_log2CPM.rds")

## generate SJARACNe inputs
log2CPM.eset$celltype <- factor(pData(log2CPM.eset)$type)
generateSJARACNeInput(input_eset = log2CPM.eset, funcType = "TF", ref = "hg", wd.src = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne", group_name = "celltype")
generateSJARACNeInput(input_eset = log2CPM.eset, funcType = "SIG", ref = "hg", wd.src = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne", group_name = "celltype")



###### ADT outputs
rawCount.eset <- CreateSparseEset(data = adt.sel, meta.data = pd[colnames(adt.sel),], add.meta = F)
saveRDS(rawCount.eset, file = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.ADT_rawCounts.rds")

# normalization
norm = 1e6
exp.norm <- sweep(exprs(rawCount.eset), 2, norm/unname(Matrix::colSums(exprs(rawCount.eset))), '*')
exp.log2 <- log(exp.norm + 1, base = 2)

pd <- pData(rawCount.eset)[colnames(exp.log2),]

log2CPM.eset <- CreateSparseEset(data = exp.log2, meta.data = pd, add.meta = F)
saveRDS(log2CPM.eset, file = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.ADT_log2CPM.rds")




