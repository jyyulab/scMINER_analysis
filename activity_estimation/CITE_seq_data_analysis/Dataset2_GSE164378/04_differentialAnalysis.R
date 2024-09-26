library(scMINER)
library(stringr)
library(NetBID2)

load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.exp.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.act.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.pro.eset")

genes.pct <- rowMeans(exprs(exp.eset) > 0)*100 
genes.sel <- names(genes.pct)[genes.pct >= 1]; length(unique(genes.sel))
exp.eset <- exp.eset[genes.sel,]

min = min(exprs(act.eset)); min
genes.pct <- rowMeans(exprs(act.eset) > min)*100 
genes.sel <- names(genes.pct)[genes.pct >= 1]; length(unique(genes.sel))
act.eset <- act.eset[genes.sel,]

genes.pct <- rowMeans(exprs(pro.eset) > 0)*100 
genes.sel <- names(genes.pct)[genes.pct >= 1]; length(unique(genes.sel))
pro.eset <- pro.eset[genes.sel,]


########### 1. differential analysis
treatment <- c("B", "CD4_T,CD8_T,DC,Mono,NK", "B", "CD4_T", "B", "CD8_T", "B", "DC",  "B", "Mono", "B", "NK",
               "CD4_T", "B,CD8_T,DC,Mono,NK", "CD4_T", "B", "CD4_T", "CD8_T", "CD4_T", "DC", "CD4_T", "Mono", "CD4_T", "NK",
               "CD8_T", "B,CD4_T,DC,Mono,NK", "CD8_T", "B", "CD8_T", "CD4_T", "CD8_T", "DC", "CD8_T", "Mono", "CD8_T", "NK",
               "DC", "B,CD4_T,CD8_T,Mono,NK", "DC", "B", "DC", "CD4_T", "DC", "CD8_T", "DC", "Mono", "DC", "NK",
               "Mono", "B,CD4_T,CD8_T,DC,NK", "Mono", "B", "Mono", "CD4_T", "Mono", "CD8_T", "Mono", "DC", "Mono", "NK",
               "NK", "B,CD4_T,CD8_T,DC,Mono", "NK", "B", "NK", "CD4_T", "NK", "CD8_T", "NK", "DC", "NK", "Mono")
compare <- matrix(treatment, ncol = 2, byrow = T)
colnames(compare) <- c("Foreground","Background"); head(compare)

DE.limma <- list(); DA.limma <- list(); DP.limma <- list()
for (i in 1:nrow(compare)) {
  cat("Comparison", i, ':', compare[i,1], 'VS', compare[i,2], '\n')
  if (str_count(compare[i,1], ",") > 0) {c_fg <- unlist(strsplit(compare[i,1], ","))} else {c_fg <- compare[i,1]}
  if (str_count(compare[i,2], ",") > 0) {c_bg <- unlist(strsplit(compare[i,2], ","))} else {c_bg <- compare[i,2]}
  
  comp_name <- paste0(compare[i,1], ".VS.", compare[i,2])
  phe_info <- pData(exp.eset)
  foreground <- row.names(phe_info)[which(phe_info$celltype %in% c_fg)]
  background <- row.names(phe_info)[which(phe_info$celltype %in% c_bg)]
  
  # DE by limma
  DE_gene_limma <- getDE.limma.2G(eset = exp.eset, G1 = foreground, G0 = background, G1_name = compare[i,1], G0_name = compare[i,2],
                                  random_effect = NULL, # a vector of character, vector or factor specifying a blocking variable
                                  verbose = F)
  
  # DA by limma
  DA_driver_limma <- getDE.limma.2G(eset = act.eset, G1 = foreground, G0 = background, G1_name = compare[i,1], G0_name = compare[i,2],
                                    random_effect = NULL, # a vector of character, vector or factor specifying a blocking variable
                                    verbose = F)
  
  # DP by limma
  DP_protein_limma <- getDE.limma.2G(eset = pro.eset, G1 = foreground, G0 = background, G1_name = compare[i,1], G0_name = compare[i,2],
                                     random_effect = NULL, # a vector of character, vector or factor specifying a blocking variable
                                     verbose = F)
  
  DE.limma[[comp_name]] <- DE_gene_limma; DA.limma[[comp_name]] <- DA_driver_limma; DP.limma[[comp_name]] <- DP_protein_limma;
}

all_comp <- names(DE.limma) # Get all comparison names
masterTable.DE_limma <- data.frame(); masterTable.DA_limma <- data.frame(); masterTable.DP_limma <- data.frame()
for (i in 1:length(all_comp)) {
  
  tmp.DE_limma <- DE.limma[[i]][, c(9,10,2,5,6,8)]
  colnames(tmp.DE_limma) <- paste0(colnames(tmp.DE_limma), ".", all_comp[i], ".DE")
  masterTable.DE_limma <- merge(masterTable.DE_limma, tmp.DE_limma, by = "row.names", all = T)
  row.names(masterTable.DE_limma) <- masterTable.DE_limma$Row.names; masterTable.DE_limma <- masterTable.DE_limma[, -1]
  
  tmp.DA_limma <- DA.limma[[i]][, c(9,10,2,5,6,8)]
  colnames(tmp.DA_limma) <- paste0(colnames(tmp.DA_limma), ".", all_comp[i], ".DA")
  masterTable.DA_limma <- merge(masterTable.DA_limma, tmp.DA_limma, by = "row.names", all = T)
  row.names(masterTable.DA_limma) <- masterTable.DA_limma$Row.names; masterTable.DA_limma <- masterTable.DA_limma[, -1]
  
  tmp.DP_limma <- DP.limma[[i]][, c(9,10,2,5,6,8)]
  colnames(tmp.DP_limma) <- paste0(colnames(tmp.DP_limma), ".", all_comp[i], ".DP")
  masterTable.DP_limma <- merge(masterTable.DP_limma, tmp.DP_limma, by = "row.names", all = T)
  row.names(masterTable.DP_limma) <- masterTable.DP_limma$Row.names; masterTable.DP_limma <- masterTable.DP_limma[, -1]
}

write.table(masterTable.DE_limma, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.DE_limma.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(masterTable.DA_limma, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.DA_limma.txt", row.names = T, col.names = T, quote = F, sep = "\t")
write.table(masterTable.DP_limma, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.DP_limma.txt", row.names = T, col.names = T, quote = F, sep = "\t")


########### 2. percentage
## 2.1 expression
exp.b <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("B")]; dim(exp.b)
exp_pc.b <- data.frame(geneSymbol = row.names(exp.b), cellType = "B", pct = rowSums(exp.b > 0)/ncol(exp.b)); head(exp_pc.b)

exp.cd4t <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("CD4_T")]; dim(exp.cd4t)
exp_pc.cd4t <- data.frame(geneSymbol = row.names(exp.cd4t), cellType = "CD4_T", pct = rowSums(exp.cd4t > 0)/ncol(exp.cd4t)); head(exp_pc.cd4t)

exp.cd8t <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("CD8_T")]; dim(exp.cd8t)
exp_pc.cd8t <- data.frame(geneSymbol = row.names(exp.cd8t), cellType = "CD8_T", pct = rowSums(exp.cd8t > 0)/ncol(exp.cd8t)); head(exp_pc.cd8t)

exp.dc <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("DC")]; dim(exp.dc)
exp_pc.dc <- data.frame(geneSymbol = row.names(exp.dc), cellType = "DC", pct = rowSums(exp.dc > 0)/ncol(exp.dc)); head(exp_pc.dc)

exp.cm <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("Mono")]; dim(exp.cm)
exp_pc.cm <- data.frame(geneSymbol = row.names(exp.cm), cellType = "Mono", pct = rowSums(exp.cm > 0)/ncol(exp.cm)); head(exp_pc.cm)

exp.nk <- exprs(exp.eset)[, pData(exp.eset)$celltype %in% c("NK")]; dim(exp.nk)
exp_pc.nk <- data.frame(geneSymbol = row.names(exp.nk), cellType = "NK", pct = rowSums(exp.nk > 0)/ncol(exp.nk)); head(exp_pc.nk)

exp_pc <- rbind(exp_pc.b, exp_pc.cd4t, exp_pc.cd8t, exp_pc.dc, exp_pc.cm, exp_pc.nk); dim(exp_pc)
write.table(exp_pc, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.pct_Exp.txt", row.names = F, col.names = T, quote = F, sep = "\t")

## 2.2 activity
act.b <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("B")]; dim(act.b)
act_pc.b <- data.frame(geneSymbol = row.names(act.b), cellType = "B", pct = rowSums(act.b > min(exprs(act.eset)))/ncol(act.b)); head(act_pc.b)

act.cd4t <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("CD4_T")]; dim(act.cd4t)
act_pc.cd4t <- data.frame(geneSymbol = row.names(act.cd4t), cellType = "CD4_T", pct = rowSums(act.cd4t > min(exprs(act.eset)))/ncol(act.cd4t)); head(act_pc.cd4t)

act.cd8t <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("CD8_T")]; dim(act.cd8t)
act_pc.cd8t <- data.frame(geneSymbol = row.names(act.cd8t), cellType = "CD8_T", pct = rowSums(act.cd8t > min(exprs(act.eset)))/ncol(act.cd8t)); head(act_pc.cd8t)

act.dc <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("DC")]; dim(act.dc)
act_pc.dc <- data.frame(geneSymbol = row.names(act.dc), cellType = "DC", pct = rowSums(act.dc > min(exprs(act.eset)))/ncol(act.dc)); head(act_pc.dc)

act.cm <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("Mono")]; dim(act.cm)
act_pc.cm <- data.frame(geneSymbol = row.names(act.cm), cellType = "Mono", pct = rowSums(act.cm > min(exprs(act.eset)))/ncol(act.cm)); head(act_pc.cm)

act.nk <- exprs(act.eset)[, pData(act.eset)$celltype %in% c("NK")]; dim(act.nk)
act_pc.nk <- data.frame(geneSymbol = row.names(act.nk), cellType = "NK", pct = rowSums(act.nk > min(exprs(act.eset)))/ncol(act.nk)); head(act_pc.nk)

act_pc <- rbind(act_pc.b, act_pc.cd4t, act_pc.cd8t, act_pc.dc, act_pc.cm, act_pc.nk); dim(act_pc)
write.table(act_pc, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.pct_Act.txt", row.names = F, col.names = T, quote = F, sep = "\t")

## 2.3 protein
pro.b <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("B")]; dim(pro.b)
pro_pc.b <- data.frame(geneSymbol = row.names(pro.b), cellType = "B", pct = rowSums(pro.b > 0)/ncol(pro.b)); head(pro_pc.b)

pro.cd4t <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("CD4_T")]; dim(pro.cd4t)
pro_pc.cd4t <- data.frame(geneSymbol = row.names(pro.cd4t), cellType = "CD4_T", pct = rowSums(pro.cd4t > 0)/ncol(pro.cd4t)); head(pro_pc.cd4t)

pro.cd8t <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("CD8_T")]; dim(pro.cd8t)
pro_pc.cd8t <- data.frame(geneSymbol = row.names(pro.cd8t), cellType = "CD8_T", pct = rowSums(pro.cd8t > 0)/ncol(pro.cd8t)); head(pro_pc.cd8t)

pro.dc <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("DC")]; dim(pro.dc)
pro_pc.dc <- data.frame(geneSymbol = row.names(pro.dc), cellType = "DC", pct = rowSums(pro.dc > 0)/ncol(pro.dc)); head(pro_pc.dc)

pro.cm <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("Mono")]; dim(pro.cm)
pro_pc.cm <- data.frame(geneSymbol = row.names(pro.cm), cellType = "Mono", pct = rowSums(pro.cm > 0)/ncol(pro.cm)); head(pro_pc.cm)

pro.nk <- exprs(pro.eset)[, pData(pro.eset)$celltype %in% c("NK")]; dim(pro.nk)
pro_pc.nk <- data.frame(geneSymbol = row.names(pro.nk), cellType = "NK", pct = rowSums(pro.nk > 0)/ncol(pro.nk)); head(pro_pc.nk)

pro_pc <- rbind(pro_pc.b, pro_pc.cd4t, pro_pc.cd8t, pro_pc.dc, pro_pc.cm, pro_pc.nk); dim(pro_pc)
write.table(pro_pc, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.pct_Pro.txt", row.names = F, col.names = T, quote = F, sep = "\t")


########### 3. mastertable including differential analysis and percentage
## prepare ms_diff
de <- masterTable.DE_limma; da <- masterTable.DA_limma; dp <- masterTable.DP_limma
row.names(dp)[!row.names(dp) %in% row.names(de)]
row.names(dp)[!row.names(dp) %in% row.names(da)]

de <- de[row.names(de) %in% row.names(dp),]; dim(de)
da <- da[row.names(da) %in% row.names(dp),]; dim(da)

ms_diff.de.b <- data.frame(geneSymbol = row.names(de), cellType = "B", modality = "Expression", de[,2], de[,c(3,9,15,21,27,33)], de[,c(4,10,16,22,28,34)], de[,c(6,12,18,24,30,36)])
ms_diff.de.cd4t <- data.frame(geneSymbol = row.names(de), cellType = "CD4_T", modality = "Expression", de[,38], de[,c(39,45,51,57,63,69)], de[,c(40,46,52,58,64,70)], de[,c(42,48,54,60,66,72)])
ms_diff.de.cd8t <- data.frame(geneSymbol = row.names(de), cellType = "CD8_T", modality = "Expression", de[,74], de[,c(75,81,87,93,99,105)], de[,c(76,82,88,94,100,106)], de[,c(78,84,90,96,102,108)])
ms_diff.de.dc <- data.frame(geneSymbol = row.names(de), cellType = "DC", modality = "Expression", de[,110], de[,c(111,117,123,129,135,141)], de[,c(112,118,124,130,136,142)], de[,c(114,120,126,132,138,144)])
ms_diff.de.cm <- data.frame(geneSymbol = row.names(de), cellType = "Mono", modality = "Expression", de[,146], de[,c(147,153,159,165,171,177)], de[,c(148,154,160,166,172,178)], de[,c(150,156,162,168,174,180)])
ms_diff.de.nk <- data.frame(geneSymbol = row.names(de), cellType = "NK", modality = "Expression", de[,182], de[,c(183,189,195,201,207,213)], de[,c(184,190,196,202,208,214)], de[,c(186,192,198,204,210,216)])

ms_diff.da.b <- data.frame(geneSymbol = row.names(da), cellType = "B", modality = "Activity", da[,2], da[,c(3,9,15,21,27,33)], da[,c(4,10,16,22,28,34)], da[,c(6,12,18,24,30,36)])
ms_diff.da.cd4t <- data.frame(geneSymbol = row.names(da), cellType = "CD4_T", modality = "Activity", da[,38], da[,c(39,45,51,57,63,69)], da[,c(40,46,52,58,64,70)], da[,c(42,48,54,60,66,72)])
ms_diff.da.cd8t <- data.frame(geneSymbol = row.names(da), cellType = "CD8_T", modality = "Activity", da[,74], da[,c(75,81,87,93,99,105)], da[,c(76,82,88,94,100,106)], da[,c(78,84,90,96,102,108)])
ms_diff.da.dc <- data.frame(geneSymbol = row.names(da), cellType = "DC", modality = "Activity", da[,110], da[,c(111,117,123,129,135,141)], da[,c(112,118,124,130,136,142)], da[,c(114,120,126,132,138,144)])
ms_diff.da.cm <- data.frame(geneSymbol = row.names(da), cellType = "Mono", modality = "Activity", da[,146], da[,c(147,153,159,165,171,177)], da[,c(148,154,160,166,172,178)], da[,c(150,156,162,168,174,180)])
ms_diff.da.nk <- data.frame(geneSymbol = row.names(da), cellType = "NK", modality = "Activity", da[,182], da[,c(183,189,195,201,207,213)], da[,c(184,190,196,202,208,214)], da[,c(186,192,198,204,210,216)])

ms_diff.dp.b <- data.frame(geneSymbol = row.names(dp), cellType = "B", modality = "Protein", dp[,2], dp[,c(3,9,15,21,27,33)], dp[,c(4,10,16,22,28,34)], dp[,c(6,12,18,24,30,36)])
ms_diff.dp.cd4t <- data.frame(geneSymbol = row.names(dp), cellType = "CD4_T", modality = "Protein", dp[,38], dp[,c(39,45,51,57,63,69)], dp[,c(40,46,52,58,64,70)], dp[,c(42,48,54,60,66,72)])
ms_diff.dp.cd8t <- data.frame(geneSymbol = row.names(dp), cellType = "CD8_T", modality = "Protein", dp[,74], dp[,c(75,81,87,93,99,105)], dp[,c(76,82,88,94,100,106)], dp[,c(78,84,90,96,102,108)])
ms_diff.dp.dc <- data.frame(geneSymbol = row.names(dp), cellType = "DC", modality = "Protein", dp[,110], dp[,c(111,117,123,129,135,141)], dp[,c(112,118,124,130,136,142)], dp[,c(114,120,126,132,138,144)])
ms_diff.dp.cm <- data.frame(geneSymbol = row.names(dp), cellType = "Mono", modality = "Protein", dp[,146], dp[,c(147,153,159,165,171,177)], dp[,c(148,154,160,166,172,178)], dp[,c(150,156,162,168,174,180)])
ms_diff.dp.nk <- data.frame(geneSymbol = row.names(dp), cellType = "NK", modality = "Protein", dp[,182], dp[,c(183,189,195,201,207,213)], dp[,c(184,190,196,202,208,214)], dp[,c(186,192,198,204,210,216)])


headers <- c("geneSymbol", "cellType", "modality", "AverageValue", "log2FC_vsRest", "log2FC_vsCellType1", "log2FC_vsCellType2", "log2FC_vsCellType3", "log2FC_vsCellType4", "log2FC_vsCellType5",
             "P_vsRest", "P_vsCellType1", "P_vsCellType2", "P_vsCellType3", "P_vsCellType4", "P_vsCellType5", "Z_vsRest", "Z_vsCellType1", "Z_vsCellType2", "Z_vsCellType3", "Z_vsCellType4", "Z_vsCellType5")

colnames(ms_diff.de.b) <- headers; colnames(ms_diff.de.cd4t) <- headers; colnames(ms_diff.de.cd8t) <- headers; colnames(ms_diff.de.dc) <- headers; colnames(ms_diff.de.cm) <- headers; colnames(ms_diff.de.nk) <- headers
colnames(ms_diff.da.b) <- headers; colnames(ms_diff.da.cd4t) <- headers; colnames(ms_diff.da.cd8t) <- headers; colnames(ms_diff.da.dc) <- headers; colnames(ms_diff.da.cm) <- headers; colnames(ms_diff.da.nk) <- headers
colnames(ms_diff.dp.b) <- headers; colnames(ms_diff.dp.cd4t) <- headers; colnames(ms_diff.dp.cd8t) <- headers; colnames(ms_diff.dp.dc) <- headers; colnames(ms_diff.dp.cm) <- headers; colnames(ms_diff.dp.nk) <- headers

ms_diff <- rbind(ms_diff.de.b, ms_diff.de.cd4t, ms_diff.de.cd8t, ms_diff.de.dc, ms_diff.de.cm, ms_diff.de.nk, ms_diff.da.b, ms_diff.da.cd4t, ms_diff.da.cd8t, ms_diff.da.dc, ms_diff.da.cm, ms_diff.da.nk, ms_diff.dp.b, ms_diff.dp.cd4t, ms_diff.dp.cd8t, ms_diff.dp.dc, ms_diff.dp.cm, ms_diff.dp.nk)
dim(ms_diff)

## prepare ms_pct
exp_pc$modality <- "Expression"; act_pc$modality <- "Activity"; pro_pc$modality <- "Protein"
ms_pct <- rbind(exp_pc, act_pc, pro_pc); dim(ms_pct)

## merge ms_diff and ms_pct
ms_diff$tag <- paste0(ms_diff$geneSymbol, "_", ms_diff$cellType, "_", ms_diff$modality)
ms_pct$tag <- paste0(ms_pct$geneSymbol, "_", ms_pct$cellType, "_", ms_pct$modality)
ms <- merge(ms_pct, ms_diff[,-c(1:3)], by = "tag"); dim(ms)
write.table(ms, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.final.txt", row.names = F, col.names = T, quote = F, sep = "\t")

overlap <- row.names(dp)[row.names(dp) %in% row.names(de)]
overlap <- overlap[overlap %in% row.names(da)]
ms.sel <- ms[ms$geneSymbol %in% overlap,]; dim(ms.sel)
write.table(ms.sel, file = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.final_overlapOnly.txt", row.names = F, col.names = T, quote = F, sep = "\t")



