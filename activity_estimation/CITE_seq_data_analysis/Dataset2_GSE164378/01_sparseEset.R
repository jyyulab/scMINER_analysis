library(Seurat)
library(scMINER)

dir <- "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation"

################### prepare the esets of geneSymbols ################### 
## 1. Protein
obj.pro <- readRDS("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/01_sparseEset.ADT_log2CPM.rds")
obj.pro <- obj.pro[,pData(obj.pro)$celltype %in% c("B","CD4_T","CD8_T","DC","Mono","NK")]; dim(obj.pro)
d.pro <- data.frame(proteinID = row.names(obj.pro), exprs(obj.pro)); dim(d.pro)

p2g <- read.table(paste0(dir, "/00_proteinID2geneSymbol.txt"), header = T, sep = "\t", quote = "", stringsAsFactors = F); head(p2g)
d.pro <- dplyr::left_join(p2g, d.pro, by = join_by("proteinID" == "proteinID"))
d.pro.sel <- d.pro[d.pro$ifTranferable %in% c("YES"),]; dim(d.pro.sel)
d.pro.sel <- d.pro.sel[,-c(1,3)]; d.pro.sel[1:10,1:10]

d.pro.agg <- aggregate(.~geneSymbol, d.pro.sel, mean, na.rm = T); d.pro.agg[1:10,1:10]
row.names(d.pro.agg) <- d.pro.agg$geneSymbol
d.pro.agg <- d.pro.agg[,-1]

pd.pro <- pData(obj.pro)[colnames(d.pro.agg),]; head(pd.pro)
table(pd.pro$celltype)

fd.pro <- data.frame(row.names = row.names(d.pro.agg), geneSymbol = row.names(d.pro.agg))

pro.eset <- CreateSparseEset(data = as.matrix(d.pro.agg), meta.data = pd.pro, feature.data = fd.pro, add.meta = F)
save(pro.eset, file = paste0(dir, "/01_sparseEset.pro.eset"))
generateMICAinput(pro.eset, filepath = paste0(dir, "/02_MICA_pro/01_micaInput.txt"))

## 2. Expression
obj.rna <- readRDS("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/01_sparseEset.RNA_log2CPM.rds")
obj.rna <- obj.rna[,pData(obj.rna)$celltype %in% c("B","CD4_T","CD8_T","DC","Mono","NK")]; dim(obj.rna)
exp.eset <- obj.rna[rowSums(exprs(obj.rna))>0,]; dim(obj.rna.sel)
save(exp.eset, file = paste0(dir, "/01_sparseEset.exp.eset"))


## 3. Activity
# 3.1 mean
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/B_17336_17336_1000/activityMatrix_mean.eset")
d.act_B <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/CD4_T_16458_16458_1000/activityMatrix_mean.eset")
d.act_CD4T <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/CD8_T_16575_16575_1000/activityMatrix_mean.eset")
d.act_CD8T <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/DC_18578_18578_1000/activityMatrix_mean.eset")
d.act_DC <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/Mono_17554_17554_1000/activityMatrix_mean.eset")
d.act_CM <- exprs(actMatrix.eset)
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/02_sjaracneDownSample1000/NK_16509_16509_1000/activityMatrix_mean.eset")
d.act_NK <- exprs(actMatrix.eset)

d.act_mean <- merge(d.act_B, d.act_CD4T, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)
d.act_mean <- merge(d.act_mean, d.act_CD8T, by = "row.names", all = T)
row.names(d.act_mean) <- d.act_mean$Row.names; d.act_mean <- d.act_mean[,-1]; dim(d.act_mean)
d.act_mean <- merge(d.act_mean, d.act_DC, by = "row.names", all = T)
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
save(act.eset, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.act.eset")
