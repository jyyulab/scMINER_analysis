library(Seurat)
library(scMINER)

dir <- "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation"

################### prepare the esets of geneSymbols ################### 
## 1. Protein
obj.pro <- readRDS("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.ADT_log2CPM.rds")
obj.pro <- obj.pro[,pData(obj.pro)$type %in% c("B","CD4T","CD8T","CM","NK")]; dim(obj.pro)
d.pro <- data.frame(proteinID = row.names(obj.pro), exprs(obj.pro)); dim(d.pro)

p2g <- read.table(paste0(dir, "/00_proteinID2geneSymbol.txt"), header = T, sep = "\t", quote = "", stringsAsFactors = F); head(p2g)
d.pro <- dplyr::left_join(p2g, d.pro, by = join_by("proteinID" == "proteinID")); dim(d.pro)
d.pro.sel <- d.pro[d.pro$ifConvertable %in% c("TRUE"),]; dim(d.pro.sel)
d.pro.sel <- d.pro.sel[,-c(1,3)]; d.pro.sel[1:10,1:10]

d.pro.agg <- aggregate(.~geneSymbol, d.pro.sel, mean, na.rm = T); d.pro.agg[1:10,1:10]
row.names(d.pro.agg) <- d.pro.agg$geneSymbol
d.pro.agg <- d.pro.agg[,-1]

pd.pro <- pData(obj.pro)[colnames(d.pro.agg),]; head(pd.pro)
table(pd.pro$type)

fd.pro <- data.frame(row.names = row.names(d.pro.agg), geneSymbol = row.names(d.pro.agg))

pro.eset <- CreateSparseEset(data = as.matrix(d.pro.agg), meta.data = pd.pro, feature.data = fd.pro, add.meta = F)
save(pro.eset, file = paste0(dir, "/01_sparseEset.pro.eset"))
generateMICAinput(pro.eset, filepath = paste0(dir, "/02_MICA_pro/01_micaInput.txt"))

## 2. Expression
obj.rna <- readRDS("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/01_sparseEset.RNA_log2CPM.rds")
obj.rna <- obj.rna[,pData(obj.rna)$type %in% c("B","CD4T","CD8T","CM","NK")]; dim(obj.rna)
exp.eset <- obj.rna[rowSums(exprs(obj.rna))>0,]; dim(exp.eset)
save(exp.eset, file = paste0(dir, "/01_sparseEset.exp.eset"))


## 3. Activity
# 3.1 mean
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne/B_13807_13807_304/activityMatrix_mean.eset")
d.act_B <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne/CD4T_17399_17399_2878/activityMatrix_mean.eset")
d.act_CD4T <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne/CD8T_14501_14501_539/activityMatrix_mean.eset")
d.act_CD8T <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne/CM_12484_12484_127/activityMatrix_mean.eset")
d.act_CM <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/02_sjaracne/NK_15039_15039_754/activityMatrix_mean.eset")
d.act_NK <- exprs(actMatrix.eset)

d.act_mean <- merge(d.act_B, d.act_CD4T, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)
d.act_mean <- merge(d.act_mean, d.act_CD8T, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)
d.act_mean <- merge(d.act_mean, d.act_CM, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)
d.act_mean <- merge(d.act_mean, d.act_NK, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)

d.act_mean2 <- data.frame(geneSymbol = gsub("_TF$|_SIG$", "", row.names(d.act_mean)), d.act_mean); d.act_mean2[1:10,1:10]
d.act <- aggregate(.~geneSymbol, d.act_mean2, mean, na.action = na.pass)

row.names(d.act) <- d.act$geneSymbol; d.act[1:10,1:10]
d.act <- d.act[,-1]

minAct <- min(d.act, na.rm = T); minAct
d.act[is.na(d.act)] <- minAct

fd.act <- data.frame(row.names = row.names(d.act), GeneSymbol = row.names(d.act)); head(fd.act)
pd.act <- pd.pro[colnames(d.act),]

act.eset <- CreateSparseEset(data = as.matrix(d.act), meta.data = pd.act, feature.data = fd.act, add.meta = F)
save(act.eset, file = paste0(dir, "/01_sparseEset.act.eset"))
