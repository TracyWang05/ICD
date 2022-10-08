#### whole score ####
rm(list=ls())
options(stringsAsFactors = F)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(Seurat)
geneset <- getGmt(file.path("data/catabolism.gmt")) 
data=readRDS("Rdata/SingleCell/data.merge.pro.rds")
DefaultAssay(data) <- "RNA"
expr <- GetAssayData(data, slot = "counts")
expr=as.matrix(expr)
score <- gsva(expr, gset.idx.list = geneset, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
score=as.data.frame(t(score))
score=score[,1]
data@meta.data$score <-score

pdf("figure/SingleCell/whole-score.pdf")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "score")
dev.off()

#### tumor score ####
rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggplot2)
library(viridis)
library(ggpointdensity)
set.seed(101)
library(future)
library(data.table)
plan("multisession", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2)
setwd("yuanwang/pan-cancer/ICD")

data=readRDS("Rdata/SingleCell/data.merge.pro.rds")
DefaultAssay(data) <- "RNA"
cell.Type <- c("epithelial cell")
sce = data[,data@meta.data$cellType_low == cell.Type]

# set.resolutions <- seq(0.5, 2, by = 0.1)
set.resolutions=0.4
T.pro <- SCTransform(sce, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T.pro <- RunPCA(T.pro, npcs = 100, verbose = FALSE)
# ElbowPlot(object = T.pro, ndims = 100)
T.pro  <- FindNeighbors(object = T.pro , dims = 1:50, verbose = FALSE)
T.pro  <- FindClusters(object = T.pro , resolution = set.resolutions, verbose = FALSE) 
# clustree(T.pro) # 0.4
T.pro  <- RunUMAP(T.pro , dims = 1:50)

cell.type.markers <- read.table(file = "data/SingleCell/iCMS.csv", header = T, stringsAsFactors = F, sep = ",")
marker.gene=c(cell.type.markers$iCMS2_Up[1:308],cell.type.markers$iCMS2_Down[1:279],
              cell.type.markers$iCMS3_Up[1:74],cell.type.markers$iCMS3_Down[1:52])
names(marker.gene)=c(rep('iCMS2_Up',308),
                     rep('iCMS2_Down',279),
                     rep('iCMS3_Up',74),
                     rep('iCMS3_Down',52))
marker.gene=as.data.frame(marker.gene)
marker.gene$type=c(rep('iCMS2_Up',308),
                   rep('iCMS2_Down',279),
                   rep('iCMS3_Up',74),
                   rep('iCMS3_Down',52))
colnames(marker.gene)=c('gene','type')
exp.matrix <- GetAssayData(T.pro, slot = "data")
index <- match(marker.gene$gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- marker.gene[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, T.pro$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) 
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$type, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)
annotation.colors <- rainbow(length(unique(cell.type.markers$type)))
names(annotation.colors) <- unique(cell.type.markers$type)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$type, 
                                               levels = unique(cell.type.markers$type))
                                 # col = list(Type = annotation.colors)
)
# pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 3), c("white", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$type, levels = unique(cell.type.markers$type))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, cluster_columns = F, cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))

a <- T.pro
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 6)) + theme(axis.text.y = element_text(size = 6))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + NoLegend()+ theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()

#####annotated the cell type
#RNA low resolution
cluster.label <- T.pro@meta.data$seurat_clusters
cluster.label <- gsub("^0$", "unknown", cluster.label)
cluster.label <- gsub("^1$", "unknown", cluster.label)
cluster.label <- gsub("^2$", "unknown", cluster.label)
cluster.label <- gsub("^3$", "unknown", cluster.label)
cluster.label <- gsub("^4$", "iCMS2", cluster.label)
cluster.label <- gsub("^5$", "iCMS2", cluster.label)
cluster.label <- gsub("^6$", "iCMS3", cluster.label)
cluster.label <- gsub("^7$", "iCMS2", cluster.label)
cluster.label <- gsub("^8$", "iCMS2", cluster.label)
cluster.label <- gsub("^9$", "iCMS3", cluster.label)
cluster.label <- gsub("^10$", "iCMS3", cluster.label)
cluster.label <- gsub("^11$", "iCMS3", cluster.label)
cluster.label <- gsub("^12$", "iCMS2", cluster.label)
cluster.label <- gsub("^13$", "iCMS2", cluster.label)
cluster.label <- gsub("^14$", "unknown", cluster.label)
cluster.label <- gsub("^15$", "iCMS2", cluster.label)
T.pro <- AddMetaData(T.pro, cluster.label, col.name = "iCMS")
# Remove cells that cannot be annotated or the cell number less than 100 cells
CMS <- subset(T.pro, subset = seurat_clusters %in% c(4:13,15))
CMS$seurat_clusters <- factor(CMS$seurat_clusters, levels = c(4:13,15))
saveRDS(CMS, "Rdata/SingleCell/CMS.rds")

# pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cellType.pro.pdf")
DimPlot(object = CMS, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")
DimPlot(object = CMS, reduction = 'umap',label = TRUE, group.by = "iCMS")
dev.off()

library(msigdbr)
library(GSVA)
library(GSEABase)
library(Seurat)
library(escape)
geneset <- getGmt(file.path("data/catabolism.gmt")) 
data=readRDS("Rdata/SingleCell/CMS.rds")
DefaultAssay(data) <- "RNA"
expr <- GetAssayData(data, slot = "counts")
expr=as.matrix(expr)
# score <- gsva(expr, gset.idx.list = geneset, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
ES <- enrichIt(obj = expr, gene.sets = geneset)
score=ES[,1]
data@meta.data$score <-score

df=cbind(Embeddings(object = data[['umap']]),FetchData(data,'iCMS'))
pdf("figure/SingleCell/CMS-score.pdf",width=8,height=6)
DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")
DimPlot(object = data, reduction = 'umap',label = TRUE, group.by = "iCMS")
ggplot(data = df,mapping = aes(x=UMAP_1,y=UMAP_2))+
  geom_pointdensity()+
  scale_color_viridis()+theme_bw()+
  theme(panel.grid=element_blank())
dev.off()

#### immune score ####
rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggplot2)
library(viridis)
library(ggpointdensity)
set.seed(101)
library(future)
library(data.table)
plan("multisession", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2)
setwd("yuanwang/pan-cancer/ICD")

data=readRDS("Rdata/SingleCell/data.merge.pro.rds")
DefaultAssay(data) <- "RNA"
cell.Type <- c("immune cell")
sce = data[,data@meta.data$cellType_low == cell.Type]

# set.resolutions <- seq(0.5, 2, by = 0.1)
set.resolutions=0.4
T.pro <- SCTransform(sce, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T.pro <- RunPCA(T.pro, npcs = 100, verbose = FALSE)
# ElbowPlot(object = T.pro, ndims = 100)
T.pro  <- FindNeighbors(object = T.pro , dims = 1:50, verbose = FALSE)
T.pro  <- FindClusters(object = T.pro , resolution = set.resolutions, verbose = FALSE) 
# clustree(T.pro) # 0.4
T.pro  <- RunUMAP(T.pro , dims = 1:50)

cell.type.markers <- read.table(file = "data/SingleCell/immune cell marker.csv", header = F, stringsAsFactors = F, sep = ",")
marker.gene=cell.type.markers
colnames(marker.gene)=c('gene','type')
exp.matrix <- GetAssayData(T.pro, slot = "data")
index <- match(marker.gene$gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- marker.gene[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, T.pro$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) 
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$type, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)
annotation.colors <- rainbow(length(unique(cell.type.markers$type)))
names(annotation.colors) <- unique(cell.type.markers$type)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$type, 
                                               levels = unique(cell.type.markers$type))
                                 # col = list(Type = annotation.colors)
)
# pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 3), c("white", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$type, levels = unique(cell.type.markers$type))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, cluster_columns = F, cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))

a <- T.pro
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 6)) + theme(axis.text.y = element_text(size = 6))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + NoLegend()+ theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()

#####annotated the cell type
#RNA low resolution
cluster.label <- T.pro@meta.data$seurat_clusters
cluster.label <- gsub("^0$", "Th1", cluster.label)
cluster.label <- gsub("^1$", "Th1", cluster.label)
cluster.label <- gsub("^2$", "pDC", cluster.label)
cluster.label <- gsub("^3$", "unknown", cluster.label)
cluster.label <- gsub("^4$", "immature B cells", cluster.label)
cluster.label <- gsub("^5$", "Treg", cluster.label)
cluster.label <- gsub("^6$", "Th1", cluster.label)
cluster.label <- gsub("^7$", "NKT", cluster.label)
cluster.label <- gsub("^8$", "activated B cells", cluster.label)
cluster.label <- gsub("^9$", "activated B cells", cluster.label)
cluster.label <- gsub("^10$", "mast cells", cluster.label)
cluster.label <- gsub("^11$", "activated B cells", cluster.label)
cluster.label <- gsub("^12$", "pDC", cluster.label)
cluster.label <- gsub("^13$", "NK56 dim", cluster.label)
cluster.label <- gsub("^14$", "Th17", cluster.label)
cluster.label <- gsub("^15$", "unknown", cluster.label)
cluster.label <- gsub("^16$", "Th17", cluster.label)
cluster.label <- gsub("^17$", "neutrophils", cluster.label)
cluster.label <- gsub("^18$", "unknown", cluster.label)
cluster.label <- gsub("^19$", "unknown", cluster.label)
T.pro <- AddMetaData(T.pro, cluster.label, col.name = "immune")
# Remove cells that cannot be annotated or the cell number less than 100 cells
Immune <- subset(T.pro, subset = seurat_clusters %in% c(0:2,4:14,16,17))
Immune$seurat_clusters <- factor(Immune$seurat_clusters, levels = c(0:2,4:14,16,17))
saveRDS(Immune, "Rdata/SingleCell/immune.rds")

pdf("figure/SingleCell/immuneCellCluster.pdf")
DimPlot(object = Immune, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")
DimPlot(object = Immune, reduction = 'umap',label = TRUE, group.by = "immune")
dev.off()

library(msigdbr)
library(GSVA)
library(GSEABase)
library(Seurat)
library(escape)
geneset <- getGmt(file.path("data/catabolism.gmt")) 
data=readRDS("Rdata/SingleCell/immune.rds")
Tcell <- subset(data, subset = seurat_clusters %in% c(0:1,5:7,14,16))
Tcell$seurat_clusters <- factor(Tcell$seurat_clusters, levels = c(0:1,5:7,14,16))
DefaultAssay(Tcell) <- "RNA"
expr <- GetAssayData(Tcell, slot = "counts")
expr=as.matrix(expr)
# score <- gsva(expr, gset.idx.list = geneset, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
ES <- enrichIt(obj = expr, gene.sets = geneset)
score=ES[,1]
Tcell@meta.data$score <-score

df=cbind(Embeddings(object = Tcell[['umap']]),FetchData(Tcell,'immune'))
pdf("figure/SingleCell/T-score.pdf",width=8,height=6)
DimPlot(object = Tcell, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")
DimPlot(object = Tcell, reduction = 'umap',label = TRUE, group.by = "immune")
ggplot(data = df,mapping = aes(x=UMAP_1,y=UMAP_2))+
  geom_pointdensity()+
  scale_color_viridis()+theme_bw()+
  theme(panel.grid=element_blank())
dev.off()
