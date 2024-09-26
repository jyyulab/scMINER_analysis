library(NetBID2)

ms <- read.table("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.final_overlapOnly.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F)

## 1. selected markers FC > 1.5: bubble plot
marker.b <- c("CD19", "FCGR2A", "CD40", "MS4A1", "TNFRSF13C", "CD22", "CCR6", "NT5E", "CXCR5", "BTLA")
marker.cd4t <- c("CD4", "CD3D", "CD3E", "CD3G", "DPP4", "ICOS", "CD28", "SELL")
marker.cd8t <- c("CD8A", "CD8B", "KLRK1")
marker.cm <- c("CD36", "CD33", "ITGAX", "CD14", "CD86", "ENTPD1", "PECAM1", "ITGAM", "CD1D", "LILRB1")
marker.nk <- c("NCAM1", "FCGR3A", "CD244", "KLRD1", "NCR1", "IL2RB", "KIR3DL1", "KLRG1", "KLRB1", "NCR3")

########################
marker.b <- c("CD19", "MS4A1", "CD22", "TNFRSF13C", "CCR6", "CD40", "CXCR5", "FCGR2A", "BTLA", "CD24")
marker.cd4t <- c("CD4", "DPP4", "ICOS", "SELL", "CD28", "CD3D", "CD3E", "CD3G")
marker.cd8t <- c("CD8A", "CD8B", "KLRK1")
marker.cm <- c("CD36", "CD33", "CD14", "CD86", "ITGAX", "CD163", "SLC3A2", "ITGB1", "GGT1", "PECAM1")
marker.nk <- c("NCAM1", "FCGR3A", "KLRD1", "NCR1", "IL2RB", "KIR3DL1", "KLRG1", "NCR3", "KLRB1", "CD244")
########################

########################
marker.b <- c("CD19", "FCGR2A", "CD40", "MS4A1", "TNFRSF13C", "CD22", "CCR6", "NT5E", "CXCR5", "BTLA", "CD24", "FCER2", "CD79B")
marker.cd4t <- c("CD4", "CD3D", "CD3E", "CD3G", "DPP4", "ICOS", "CD28", "SELL")
marker.cd8t <- c("CD8A", "CD8B", "KLRK1")
marker.cm <- c("CD36", "CD33", "ITGAX", "CD14", "CD86", "ENTPD1", "PECAM1", "ITGAM", "CD1D", "LILRB1", "CD163", "SLC3A2", "ITGB1", "ICAM1", "GGT1", "C5AR1", "ITGAE", "CD44", "ITGAL", "HAVCR2", "TNFSF14", "ITGB2", "IL7R", "TNFSF9", "FAS", "CD58", "ICOSLG", "LAMP1")
marker.nk <- c("NCAM1", "FCGR3A", "CD244", "KLRD1", "NCR1", "IL2RB", "KIR3DL1", "KLRG1", "KLRB1", "NCR3")
########################

######################## by seurat
marker.b <- c("CD19", "CD22", "CR2", "MS4A1", "FCGR2A", "CD40", "ENTPD1", "CR1", "CCR6", "TNFRSF13C")
marker.cd4t <- c("CD4", "CD27", "CD3D", "CD3E", "CD3G", "DPP4", "IL2RA", "ICOS", "CD5", "IL6R")
marker.cd8t <- c("CD8A", "CD8B", "KLRK1", "PECAM1", "NT5E", "CD27")
marker.cm <- c("CD36", "CD33", "CD86", "CD14", "FCGR1A",  "CLEC12A", "ITGAX", "SIGLEC7", "ENTPD1", "SELP")
marker.nk <- c("SIGLEC7", "NCAM1", "ITGAX", "KLRD1", "CD244", "FCGR3A", "KIR3DL1", "NCR1", "KLRB1", "KLRK1")
########################

ms.sel <- ms[ms$geneSymbol %in% c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk),]; dim(ms.sel)
ms.sel$modalityNew <- factor(ms.sel$modality, levels = c("Protein", "Activity", "Expression"))

sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="red", space ="Lab")

## z score
p1 <- ggplot(ms.sel, aes(x = geneSymbol, y = cellType, colour = Z_vsRest, size = pct)) + geom_point() + facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p1
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/02_bubblePlot.selectedMarker_FC1.5.Zscore_byMin.pdf",
       p1, width = 12, height = 6, units = "in", useDingbats = F)

## fold change
p2 <- ggplot(ms.sel, aes(x = geneSymbol, y = cellType, colour = log2FC_vsRest, size = pct)) + geom_point() + facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p2
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/02_bubblePlot.selectedMarker_FC1.5.log2FC_byMin.pdf",
       p2, width = 12, height = 6, units = "in", useDingbats = F)


## expression average
markers <- data.frame(celltype = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk),
                      markers = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk), weight = 1)

load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.exp.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.act.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.pro.eset")


draw.marker.bbp<-function(ref = NULL,input_eset, feature='geneSymbol',group_name="ClusterRes"){
  
  #exp<-apply(exprs(eset),2,std)
  #filter reference marker sets
  
  if (!feature%in%colnames(fData(input_eset))) stop('Please check your feature!')
  colnames(ref)<-c("celltype","markers","weight")
  ref<-dplyr::filter(ref,markers%in%fData(input_eset)[,feature])
  indx<-which(fData(input_eset)[,feature]%in%ref$markers)
  if(length(indx)==0) stop("No genes from the reference list could be found in data!","\n")
  
  exp<-as.matrix(exprs(input_eset))[indx,]
  rownames(exp)<-fData(input_eset)[,feature][indx]
  
  celltypes<-unique(ref$celltype)
  
  ac<-matrix(NA,nrow=ncol(exp),ncol=length(celltypes),dimnames = list(colnames(exp),celltypes))
  for(i in 1:length(celltypes)){
    cat(i,"\n")
    ref.sel<-dplyr::filter(ref,celltype==celltypes[i])
    n <- length(unique(ref.sel$markers))
    
    if(n>1){
      mat<-t(exp[ref.sel$markers,])%*%as.vector(ref.sel$weight)
      ac[,i]<-mat[,1]/n
    }else if (n==1){
      ac[,i]<-exp[ref.sel$markers,]
    }
  }
  
  ac_norm<-apply(ac,2,scale) #column normalization
  
  n_mtx<-(ac>0.5)
  df_n<-data.frame(label=pData(input_eset)[,group_name],n_mtx)
  df_n<-aggregate(.~label,df_n,mean)
  library(reshape2)
  df_n_melt<-melt(df_n,id.vars = "label")
  
  df<-data.frame(label=pData(input_eset)[,group_name],ac_norm);
  df<-df[,colSums(is.na(df))<nrow(df)];#remove NA columns
  df<-aggregate(.~label,df,mean)
  input<-t(apply(df[,-1],1,scale))#row normalization
  input<-as.data.frame(cbind(df[,1],input))
  rownames(input)<-rownames(df)
  colnames(input)<-colnames(df)
  df_melt<-melt(input,id.vars = "label")
  
  
  if(all(df_melt[,c(1,2)]==df_n_melt[,c(1,2)])){
    
    d<-cbind(df_melt,df_n_melt[,3])
    colnames(d)<-c("Cluster", "CellType", "MarkerScore","ExpressionPercentage")
    d$Cluster<-as.factor(d$Cluster)
    d$MarkerScore <- as.numeric(d$MarkerScore)
    d$ExpressionPercentage <- as.numeric(d$ExpressionPercentage)
  }
  return(d)
}


d.pro <- draw.marker.bbp(ref = markers, input_eset = pro.eset, feature = "geneSymbol", group_name = "type"); head(d.pro)
d.act <- draw.marker.bbp(ref = markers, input_eset = act.eset, feature = "GeneSymbol", group_name = "type"); head(d.act)
d.exp <- draw.marker.bbp(ref = markers, input_eset = exp.eset, feature = "geneSymbol", group_name = "type"); head(d.exp)

d.pro$modality <- "Protein"
d.act$modality <- "Activity"
d.exp$modality <- "Expression"

d <- rbind(d.pro, d.act, d.exp); dim(d)
d$modalityNew <- factor(d$modality, levels = c("Protein", "Activity", "Expression"))

sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="red", space ="Lab")
p3 <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p3
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/02_bubblePlot.selectedMarker_FC1.5.quantification_byMin.pdf",
       p3, width = 12, height = 6, units = "in", useDingbats = F)
