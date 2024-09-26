library(Seurat)

################################################# 1. read inputs
### protein
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.pro.eset")
d <- exprs(pro.eset); pd <- pData(pro.eset)

obj.pro <- CreateSeuratObject(counts=d)
Idents(obj.pro) <- pd$type; table(Idents(obj.pro))

obj.pro@assays$RNA@data <- obj.pro@assays$RNA@counts

#### expression
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.exp.eset")
d <- exprs(exp.eset); pd <- pData(exp.eset)

obj.exp <- CreateSeuratObject(counts=d)
Idents(obj.exp) <- pd$type; table(Idents(obj.exp))

obj.exp@assays$RNA@data <- obj.exp@assays$RNA@counts

#### activity
load("/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/01_sparseEset.act.eset")
d <- exprs(act.eset); pd <- pData(act.eset)

obj.act <- CreateSeuratObject(counts=d)
Idents(obj.act) <- pd$type; table(Idents(obj.act))

obj.act@assays$RNA@data <- obj.act@assays$RNA@counts


################################################# 2. specify markers
marker.b <- c("CD19", "FCGR2A", "CD40", "MS4A1", "TNFRSF13C", "CD22", "CCR6", "NT5E", "CXCR5", "BTLA")
marker.cd4t <- c("CD4", "CD3D", "CD3E", "CD3G", "DPP4", "ICOS", "CD28", "SELL")
marker.cd8t <- c("CD8A", "CD8B", "KLRK1")
marker.cm <- c("CD36", "CD33", "ITGAX", "CD14", "CD86", "ENTPD1", "PECAM1", "ITGAM", "CD1D", "LILRB1")
marker.nk <- c("NCAM1", "FCGR3A", "CD244", "KLRD1", "NCR1", "IL2RB", "KIR3DL1", "KLRG1", "KLRB1", "NCR3")

features <- unique(c(marker.b, marker.cd4t, marker.cd8t, marker.cm, marker.nk))

################################################# 3. bubble plot
DotPlot.new <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                         scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  return(data.plot)
}

d.pro <- DotPlot.new(obj.pro, features = features)
d.act <- DotPlot.new(obj.act, features = features)
d.exp <- DotPlot.new(obj.exp, features = features)

d.pro$modality <- "Protein"
d.act$modality <- "Activity"
d.exp$modality <- "Expression"

d <- rbind(d.pro, d.act, d.exp); dim(d)
d$modalityNew <- factor(d$modality, levels = c("Protein", "Activity", "Expression"))

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="red", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.red.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="blue", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.blue.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#fdc70c", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.yellow.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="green", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.green.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#fff33b", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.orange.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

##
sc <- scale_color_gradient2(midpoint=0, low="grey", mid="white", high="#1b8a5a", space ="Lab")
p <- ggplot(d, aes(x = features.plot, y = id, colour = avg.exp.scaled, size = pct.exp)) + geom_point() + 
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
ggsave(filename = "/research_jude/rgs01_jude/groups/yu3grp/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/CITEseq/02_PBMC/03_evaluation/03_diffAnalysis/05_heatmap_bubbles_scaledBySeurat.darkblue.pdf",
       p, width = 12, height = 6, units = "in", useDingbats = F)

