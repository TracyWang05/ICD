library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
set.seed(101)
library(openxlsx)
library(future)
plan("multisession", workers = 10) 
options(future.globals.maxSize = 50000 * 1024^2) # set 50G RAM
setwd('yuanwang/pan-cancer/ICD')
source("scripts/SingleCell/Combined.P.FC.R")
source("scripts/SingleCell/ratio.plot.R")

#### Harmony corrected result
data.merge <- readRDS("Rdata/SingleCell/data.merge.harmony.PC40.SCT.feature3000.rds")
DefaultAssay(data.merge) <- "RNA"

pdf("figure/SingleCell/2.Cluster/cluster.pdf")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "seurat_clusters")+NoLegend()
dev.off()

#### Plot the ratio of each cell type and sample situation
pdf("figure/SingleCell/2.Cluster/cluster.number.ratio.pdf", height = 4, width = 7)
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters")
dev.off()

expMatrix <- GetAssayData(data.merge, slot = "scale.data")
highVariableGenes <- VariableFeatures(data.merge)
expMatrix.high <- expMatrix[highVariableGenes,]
meanExpCluster <- apply(expMatrix.high, 1, function(x){
  mean.value <- tapply(x, data.merge$seurat_clusters, mean)
  return(mean.value)
})

corrMatrix <- (1- cor(t(meanExpCluster), method="pearson"))/2
library(ape)
## dd <- dist(M)
hc <- hclust(as.dist(corrMatrix),method="complete")
pdf("figure/SingleCell/2.Cluster/clustersimilarity.HVG.pdf")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)
dev.off()

#### Differential expression
Idents(data.merge) <- data.merge$seurat_clusters
cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "seurat_clusters", test.use = "MAST", latent.vars = "orig.ident")
cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
saveRDS(cluster.sig.markers, file = "Rdata/SingleCell/2.Cluster/AnnotateCellType/cluster.sig.markers.rds")
write.table(cluster.sig.markers, file = "data/SingleCell/2.Cluster/AnnotateCellType/cluster.sig.markers.txt", quote = F, sep = "\t", row.names = F)
cluster.sig.markers.top <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(cluster.sig.markers.top, file = "data/SingleCell/2.Cluster/AnnotateCellType/cluster.sig.markers.top30.txt", col.names = T, row.names = F, sep = "\t", quote = F)
top.genes <- cluster.sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cluster.topgenes.pdf")
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

####Method 1: The expression of the classic marker
cell.type.markers <- read.table(file = "data/SingleCell/simplified cell marker.csv", header = T, stringsAsFactors = F, sep = ",")
exp.matrix <- GetAssayData(data.merge, slot = "data")
index <- match(cell.type.markers$Gene, rownames(exp.matrix))
gene.matrix <- exp.matrix[na.omit(index),]
cell.type.markers <- cell.type.markers[which(!is.na(index)),]

# Calculate the average expression of the cell type of each cluster
cluster.score <- apply(gene.matrix, 1, function(x){
  a <- tapply(x, data.merge$seurat_clusters, mean)
  return(a)
})
cluster.score.normailzed <- decostand(cluster.score, "range", 2) 
cellType.cluster.score <- apply(cluster.score, 1, function(x){
  a <- tapply(x, cell.type.markers$Cell.Type, mean)
  return(a)
})
cellType.cluster.score.normailzed <- decostand(cellType.cluster.score, "range", 1)
annotation.colors <- rainbow(length(unique(cell.type.markers$Cell.Type)))
names(annotation.colors) <- unique(cell.type.markers$CellType)
row.annotations <- rowAnnotation(Type = factor(cell.type.markers$Cell.Type, 
                                               levels = unique(cell.type.markers$Cell.Type))
                                 # col = list(Type = annotation.colors)
                                 )
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cluster.signature.expression.pdf")
col_fun1 <- colorRamp2(c(0, 3), c("white", "#ff5a36"))
col_fun2 <- colorRamp2(c(0, 0.5, 1), c("#1e90ff", "white", "#ff5a36"))

row_split <- factor(cell.type.markers$Cell.Type, levels = unique(cell.type.markers$Cell.Type))
Heatmap(t(cluster.score.normailzed), col = col_fun2, row_split = row_split, width = unit(10, "cm"), height = unit(12, "cm"), cluster_columns = F, cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))
Heatmap(cellType.cluster.score, name = "Expression", col = col_fun1, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = T, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
Heatmap(cellType.cluster.score.normailzed, col = col_fun2, width = unit(8, "cm"), height = unit(8, "cm"), cluster_columns = T , cluster_rows = F, show_column_names = T, show_row_names = T, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6), 
        heatmap_legend_param = list(
          title = "Expression", at = c(0, 1), 
          labels = c("min", "max")))

a <- data.merge
a$seurat_clusters <- factor(a$seurat_clusters, levels = rownames(cluster.score.normailzed))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + theme(axis.text.x = element_text(size = 6)) + theme(axis.text.y = element_text(size = 6))
DotPlot(a, cols = c("#1e90ff", "#ff5a36"), features = cell.type.markers$Gene, group.by = "seurat_clusters", dot.scale = 4) + RotatedAxis() + NoLegend()+ theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 8))
dev.off()

#####annotated the cell type
#RNA low resolution
cluster.label <- data.merge@meta.data$seurat_clusters
cluster.label <- gsub("^0$", "immune cell", cluster.label)
cluster.label <- gsub("^1$", "immune cell", cluster.label)
cluster.label <- gsub("^2$", "immune cell", cluster.label)
cluster.label <- gsub("^3$", "epithelial cell", cluster.label)
cluster.label <- gsub("^4$", "immune cell", cluster.label)
cluster.label <- gsub("^5$", "immune cell", cluster.label)
cluster.label <- gsub("^6$", "immune cell", cluster.label)
cluster.label <- gsub("^7$", "epithelial cell", cluster.label)
cluster.label <- gsub("^8$", "epithelial cell", cluster.label)
cluster.label <- gsub("^9$", "epithelial cell", cluster.label)
cluster.label <- gsub("^10$", "unknown", cluster.label)
cluster.label <- gsub("^11$", "immune cell", cluster.label)
cluster.label <- gsub("^12$", "epithelial cell", cluster.label)
cluster.label <- gsub("^13$", "epithelial cell", cluster.label)
cluster.label <- gsub("^14$", "stromal cell", cluster.label)
cluster.label <- gsub("^15$", "immune cell", cluster.label)
cluster.label <- gsub("^16$", "epithelial cell", cluster.label)
cluster.label <- gsub("^17$", "immune cell", cluster.label)
cluster.label <- gsub("^18$", "immune cell", cluster.label)
cluster.label <- gsub("^19$", "stromal cell", cluster.label)
cluster.label <- gsub("^20$", "immune cell", cluster.label)
cluster.label <- gsub("^21$", "unknown", cluster.label)
cluster.label <- gsub("^22$", "immune cell", cluster.label)
cluster.label <- gsub("^23$", "unknown", cluster.label)
data.merge <- AddMetaData(data.merge, cluster.label, col.name = "cellType_low")
# Remove cells that cannot be annotated or the cell number less than 100 cells
data.merge.pro <- subset(data.merge, subset = seurat_clusters %in% c(0:9, 11:20,22))
data.merge.pro$seurat_clusters <- factor(data.merge.pro$seurat_clusters, levels = c(0:9, 11:20,22))
saveRDS(data.merge.pro, "Rdata/SingleCell/data.merge.pro.rds")

pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cellType.pro.pdf")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")+NoLegend()
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
dev.off()

##Plot--- celltype marker plot
plot.cellType <- c("immune cell", "epithelial cell", "stromal cell")
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cellType.ratio.pdf")
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "cellType_low", angle = 60)
dev.off()

#### Cell type specific gene
# data.merge=readRDS('Rdata/SingleCell/data.merge.pro.rds')
Idents(data.merge) <- data.merge$cellType_low
idents <- as.character(levels(data.merge))
cellType.all.markers <- FindAllMarkers(data.merge, 
                                       group.by = "cellType_low", 
                                       logfc.threshold = 0, 
                                       min.pct = 0.1, 
                                       test.use = "MAST", 
                                       latent.vars = "orig.ident")
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.all.markers$cluster == x)
  DEGs <- cellType.all.markers[index,]
  DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
  DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
  DEGs <- rbind(DEGs.up, DEGs.down)
  return(DEGs)
})
write.xlsx(saveFormat, file = "data/SingleCell/2.Cluster/AnnotateCellType/celltype.all.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.all.markers, file = "Rdata/SingleCell/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")

cellType.all.markers=readRDS("Rdata/SingleCell/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
#require logfc.threshold = 0.25 & p_val_adj < 0.05
cellType.sig.DEGs <- cellType.all.markers %>% filter(abs(avg_log2FC) >1 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
saveFormat <- lapply(idents, function(x){
  index <- which(cellType.sig.DEGs$cluster == x)
  DEGs <- cellType.sig.DEGs[index,]
  DEGs <- DEGs %>% arrange(desc(avg_log2FC))
  return(DEGs)
})
write.xlsx(saveFormat, file = "data/SingleCell/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellType.sig.DEGs, file = "Rdata/SingleCell/2.Cluster/AnnotateCellType/cellType.sig.pos.DEGs.rds")
top.genes <- cellType.sig.DEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/cellType.topgenes.pdf")
DoHeatmap(data.merge, features = unique(top.genes$gene), size = 2) + NoLegend()
dev.off()

Tumor.info <- data.merge@meta.data[which(data.merge@meta.data$cellType=="epithelial cell"),]
cell.number <- as.data.frame(table(Tumor.info$orig.ident))
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/Tumor.patient.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = rainbow(length(unique(Tumor.info$orig.ident))),
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

cellType.ratio <- as.data.frame(table(data.merge$cellType_low))
cellType.ratio$Type <- rep("immune cells", nrow(cellType.ratio))
idx <- which(cellType.ratio$Var1 %in% "stromal cell")
cellType.ratio$Type[idx] <- "stromal cell"
idx <- which(cellType.ratio$Var1 %in% "epitheial cell")
cellType.ratio$Type[idx] <- "epithelial cell"
group.ratio <- tapply(cellType.ratio$Freq, cellType.ratio$Type, sum)
group.ratio <- data.frame(Type = names(group.ratio), Ratio = group.ratio)
labs <- paste0(group.ratio$Type, " (", round(group.ratio$Ratio/sum(group.ratio$Ratio), 4)*100, "%)")
pdf("figure/SingleCell/2.Cluster/AnnotateCellType/lineage.ratio.pdf")
p <- ggpie(group.ratio, "Ratio", label = labs,
           fill = "Type", color = "white", lab.pos = "in",
           palette = "npgs")
print(p)
dev.off()

#### functional enrichment ####
cellType.all.DEGs <- readRDS("Rdata/SingleCell/2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
cellType.sig.DEGs <- cellType.all.DEGs %>% filter(abs(avg_log2FC) >1 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
Tumor.DEGs <- cellType.sig.DEGs[cellType.sig.DEGs$cluster=="epithelial cell",]
geneList <- Tumor.DEGs$avg_log2FC
names(geneList) <- Tumor.DEGs$gene # 587 genes

a <- data.frame(geneName = Tumor.DEGs$gene, Rank = Tumor.DEGs$avg_log2FC)
write.csv(a, file = "data/SingleCell/2.Cluster/AnnotateCellType/GSEA/Tumor.DEGs.csv")

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

gene=a$geneName
tmp <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
entrez=tmp$ENTREZID
x <-enrichPathway(gene=entrez,pvalueCutoff=0.05, readable=T)
x1 <- pairwise_termsim(x)
pdf("figure/SingleCell/DEG-enrichment.pdf",width = 13,height = 13)
emapplot(x1)
dev.off()
