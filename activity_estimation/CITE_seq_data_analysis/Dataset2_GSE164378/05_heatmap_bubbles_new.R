library(NetBID2)

ms <- read.table("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/01_masterTable_scMINER.final_overlapOnly.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F)

## 1. selected markers FC > 1.5: bubble plot
marker.b <- c("MS4A1", "CR2", "CD72", "CD19", "CD22", "TNFRSF13C",  "CD40", "BTLA", "CXCR5", "NT5E")
marker.cd4t <- c("CD4", "IL7R", "CD28")
marker.cd8t <- c("CD8A", "CD8B")
marker.dc <- c("IL3RA", "TFRC","NRP1","IL6R","THBD","CCR5","FLT3","ITGA4")
marker.cm <- c("CD36", "PVR", "CD86", "FCGR1A", "TREM1", "ITGAM", "ITGAX", "CD14", "ANPEP", "CD93")
marker.nk <- c("NCR1", "FCGR3A", "KLRB1", "NCR3", "NCAM1", "IL2RB", "CD38", "SPN", "TIGIT")

## expression average
markers <- data.frame(celltype = c(marker.b, marker.cd4t, marker.cd8t, marker.dc, marker.cm, marker.nk),
                      markers = c(marker.b, marker.cd4t, marker.cd8t, marker.dc, marker.cm, marker.nk), weight = 1)

load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.exp.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.act.eset")
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/01_sparseEset.pro.eset")


draw.marker.bbp<-function(ref = NULL,input_eset,
                          feature='geneSymbol',group_name="ClusterRes"){
  
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


d.pro <- draw.marker.bbp(ref = markers, input_eset = pro.eset, feature = "geneSymbol", group_name = "celltype"); head(d.pro)
d.act <- draw.marker.bbp(ref = markers, input_eset = act.eset, feature = "GeneSymbol", group_name = "celltype"); head(d.act)
d.exp <- draw.marker.bbp(ref = markers, input_eset = exp.eset, feature = "geneSymbol", group_name = "celltype"); head(d.exp)

d.pro$modality <- "Protein"
d.act$modality <- "Activity"
d.exp$modality <- "Expression"

d <- rbind(d.pro, d.act, d.exp); dim(d)
d$modalityNew <- factor(d$modality, levels = c("Protein", "Activity", "Expression"))

sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="red", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.red.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="blue", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.blue.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)


sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#fdc70c", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.yellow.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)


sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="green", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.green.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)


sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#fff33b", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.orange.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)


sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#1b8a5a", space ="Lab")
p <- ggplot(d, aes(x = CellType, y = Cluster, colour = MarkerScore, size = ExpressionPercentage)) + geom_point() + 
  facet_wrap(vars(modalityNew), nrow = 3) + sc +
  labs(x = "", y = "") +
  theme(legend.position = "right", legend.box = "vertical", legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", size = 1, fill = NA),
        axis.title = element_text(size = 8, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 6, face = "bold", hjust = 1, color = "black", angle = 45),
        axis.text.y = element_text(size = 6, face = "bold", hjust = 1, color = "black")) +
  scale_x_discrete(limits = c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))
p
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/04_Hao/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_new.darkblue.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)


