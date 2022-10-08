# load package
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
library(future)
library(data.table)
plan("multisession", workers = 10) 
options(future.globals.maxSize = 100000 * 1024^2)
setwd("yuanwang/pan-cancer/ICD")

#### 12 samples
T1.data <- fread('data/SingleCell/GSE166555/GSM5075660_p007t.tsv.gz')
T1.data=tibble::column_to_rownames(T1.data,'gene')

T2.data <- fread('data/SingleCell/GSE166555/GSM5075662_p008t.tsv.gz')
T2.data=tibble::column_to_rownames(T2.data,'gene')

T3.data <- fread('data/SingleCell/GSE166555/GSM5075665_p009t1.tsv.gz')
T3.data=tibble::column_to_rownames(T3.data,'gene')

T4.data <- fread('data/SingleCell/GSE166555/GSM5075668_p012t.tsv.gz')
T4.data=tibble::column_to_rownames(T4.data,'gene')

T5.data <- fread('data/SingleCell/GSE166555/GSM5075670_p013t.tsv.gz')
T5.data=tibble::column_to_rownames(T5.data,'gene')

T6.data <- fread('data/SingleCell/GSE166555/GSM5075672_p014t.tsv.gz')
T6.data=tibble::column_to_rownames(T6.data,'gene')

T7.data <- fread('data/SingleCell/GSE166555/GSM5075674_p016t.tsv.gz')
T7.data=tibble::column_to_rownames(T7.data,'gene')

T8.data <- fread('data/SingleCell/GSE166555/GSM5075676_p017t.tsv.gz')
T8.data=tibble::column_to_rownames(T8.data,'gene')

T9.data <- fread('data/SingleCell/GSE166555/GSM5075678_p020t.tsv.gz')
T9.data=tibble::column_to_rownames(T9.data,'gene')

T10.data <- fread('data/SingleCell/GSE166555/GSM5075680_p021t.tsv.gz')
T10.data=tibble::column_to_rownames(T10.data,'gene')

T11.data <- fread('data/SingleCell/GSE166555/GSM5075682_p025t.tsv.gz')
T11.data=tibble::column_to_rownames(T11.data,'gene')

T12.data <- fread('data/SingleCell/GSE166555/GSM5075683_p026t.tsv.gz')
T12.data=tibble::column_to_rownames(T12.data,'gene')

#### Preliminary filtration
#min.cell >3 & min.features >200                        
T1 <- CreateSeuratObject(counts = T1.data,
                         project = "T1",
                         min.cells = 3,
                         min.features = 200)
T2 <- CreateSeuratObject(counts = T2.data,
                         project = "T2",
                         min.cells = 3,
                         min.features = 200) 
T3 <- CreateSeuratObject(counts = T3.data,
                         project = "T3",
                         min.cells = 3,
                         min.features = 200)
T4 <- CreateSeuratObject(counts = T4.data,
                         project = "T4",
                         min.cells = 3,
                         min.features = 200) 
T5 <- CreateSeuratObject(counts = T5.data,
                         project = "T5",
                         min.cells = 3,
                         min.features = 200)
T6 <- CreateSeuratObject(counts = T6.data,
                         project = "T6",
                         min.cells = 3,
                         min.features = 200)
T7 <- CreateSeuratObject(counts = T7.data,
                         project = "T7",
                         min.cells = 3,
                         min.features = 200)
T8 <- CreateSeuratObject(counts = T8.data,
                         project = "T8",
                         min.cells = 3,
                         min.features = 200)
T9 <- CreateSeuratObject(counts = T9.data,
                         project = "T9",
                         min.cells = 3,
                         min.features = 200)
T10 <- CreateSeuratObject(counts = T10.data,
                         project = "T10",
                         min.cells = 3,
                         min.features = 200)
T11 <- CreateSeuratObject(counts = T11.data,
                         project = "T11",
                         min.cells = 3,
                         min.features = 200)
T12 <- CreateSeuratObject(counts = T12.data,
                         project = "T12",
                         min.cells = 3,
                         min.features = 200)

#### Mitochondrial gene ratio
T1[["percent.mt"]] <- PercentageFeatureSet(T1, pattern = "^MT-")
T2[["percent.mt"]] <- PercentageFeatureSet(T2, pattern = "^MT-")
T3[["percent.mt"]] <- PercentageFeatureSet(T3, pattern = "^MT-")
T4[["percent.mt"]] <- PercentageFeatureSet(T4, pattern = "^MT-")
T5[["percent.mt"]] <- PercentageFeatureSet(T5, pattern = "^MT-")
T6[["percent.mt"]] <- PercentageFeatureSet(T6, pattern = "^MT-")
T7[["percent.mt"]] <- PercentageFeatureSet(T7, pattern = "^MT-")
T8[["percent.mt"]] <- PercentageFeatureSet(T8, pattern = "^MT-")
T9[["percent.mt"]] <- PercentageFeatureSet(T9, pattern = "^MT-")
T10[["percent.mt"]] <- PercentageFeatureSet(T10, pattern = "^MT-")
T11[["percent.mt"]] <- PercentageFeatureSet(T11, pattern = "^MT-")
T12[["percent.mt"]] <- PercentageFeatureSet(T12, pattern = "^MT-")

######################### Detect the resolution parameters of each sample cluster. After the parameters are determined, you can block them without executing [test]
# set.resolutions <- seq(0.5, 2, by = 0.1)
# pdf(file = "figure/SingleCell/1.QualityControl/PCA-test.pdf")

set.resolutions <- 0.8
T1.pro <- subset(T1, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T1.pro <- SCTransform(T1.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T1.pro <- RunPCA(T1.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T1.pro, ndims = 100)
T1.pro  <- FindNeighbors(object = T1.pro , dims = 1:50, verbose = FALSE)
T1.pro  <- FindClusters(object = T1.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T1.pro)
T1.pro  <- RunUMAP(T1.pro , dims = 1:50)
# T1.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T1.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

# set.resolutions <- 1.9
T2.pro <- subset(T2, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000)
T2.pro <- SCTransform(T2.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T2.pro <- RunPCA(T2.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T2.pro, ndims = 100)
T2.pro  <- FindNeighbors(object = T2.pro , dims = 1:50, verbose = FALSE)
T2.pro  <- FindClusters(object = T2.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T2.pro)
T2.pro  <- RunUMAP(T2.pro , dims = 1:50)
# T2.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T2.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.8
T3.pro <- subset(T3, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T3.pro <- SCTransform(T3.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T3.pro <- RunPCA(T3.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T3.pro, ndims = 100)
T3.pro  <- FindNeighbors(object = T3.pro , dims = 1:50, verbose = FALSE)
T3.pro  <- FindClusters(object = T3.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T3.pro)
T3.pro  <- RunUMAP(T3.pro , dims = 1:50)
# T3.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T3.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.5
T4.pro <- subset(T4, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T4.pro <- SCTransform(T4.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T4.pro <- RunPCA(T4.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T4.pro, ndims = 100)
T4.pro  <- FindNeighbors(object = T4.pro , dims = 1:50, verbose = FALSE)
T4.pro  <- FindClusters(object = T4.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T4.pro)
T4.pro  <- RunUMAP(T4.pro , dims = 1:50)
# T4.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T4.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.6
T5.pro <- subset(T5, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T5.pro <- SCTransform(T5.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T5.pro <- RunPCA(T5.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T5.pro, ndims = 100)
T5.pro  <- FindNeighbors(object = T5.pro , dims = 1:50, verbose = FALSE)
T5.pro  <- FindClusters(object = T5.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T5.pro)
T5.pro  <- RunUMAP(T5.pro , dims = 1:50)
# T5.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T5.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.5
T6.pro <- subset(T6, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000)
T6.pro <- SCTransform(T6.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T6.pro <- RunPCA(T6.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T6.pro, ndims = 100)
T6.pro  <- FindNeighbors(object = T6.pro , dims = 1:50, verbose = FALSE)
T6.pro  <- FindClusters(object = T6.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T6.pro)
T6.pro  <- RunUMAP(T6.pro , dims = 1:50)
# T6.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T6.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.6
T7.pro <- subset(T7, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T7.pro <- SCTransform(T7.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T7.pro <- RunPCA(T7.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T7.pro, ndims = 100)
T7.pro  <- FindNeighbors(object = T7.pro , dims = 1:50, verbose = FALSE)
T7.pro  <- FindClusters(object = T7.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T7.pro)
T7.pro  <- RunUMAP(T7.pro , dims = 1:50)
# T7.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T7.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.6
T8.pro <- subset(T8, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T8.pro <- SCTransform(T8.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T8.pro <- RunPCA(T8.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T8.pro, ndims = 100)
T8.pro  <- FindNeighbors(object = T8.pro , dims = 1:50, verbose = FALSE)
T8.pro  <- FindClusters(object = T8.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T8.pro)
T8.pro  <- RunUMAP(T8.pro , dims = 1:50)
# T8.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T8.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.9
T9.pro <- subset(T9, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T9.pro <- SCTransform(T9.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T9.pro <- RunPCA(T9.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T9.pro, ndims = 100)
T9.pro  <- FindNeighbors(object = T9.pro , dims = 1:50, verbose = FALSE)
T9.pro  <- FindClusters(object = T9.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T9.pro)
T9.pro  <- RunUMAP(T9.pro , dims = 1:50)
# T9.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T9.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.8
T10.pro <- subset(T10, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000)
T10.pro <- SCTransform(T10.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T10.pro <- RunPCA(T10.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T10.pro, ndims = 100)
T10.pro  <- FindNeighbors(object = T10.pro , dims = 1:50, verbose = FALSE)
T10.pro  <- FindClusters(object = T10.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T10.pro)
T10.pro  <- RunUMAP(T10.pro , dims = 1:50)
# T10.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T10.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 2
T11.pro <- subset(T11, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T11.pro <- SCTransform(T11.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T11.pro <- RunPCA(T11.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T11.pro, ndims = 100)
T11.pro  <- FindNeighbors(object = T11.pro , dims = 1:50, verbose = FALSE)
T11.pro  <- FindClusters(object = T11.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T11.pro)
T11.pro  <- RunUMAP(T11.pro , dims = 1:50)
# T11.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T11.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

set.resolutions <- 0.7
T12.pro <- subset(T12, subset = nFeature_RNA > 200 & percent.mt < 10 & nCount_RNA > 1000 & nFeature_RNA < 6000) 
T12.pro <- SCTransform(T12.pro, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
T12.pro <- RunPCA(T12.pro, npcs = 100, verbose = FALSE)
ElbowPlot(object = T12.pro, ndims = 100)
T12.pro  <- FindNeighbors(object = T12.pro , dims = 1:50, verbose = FALSE)
T12.pro  <- FindClusters(object = T12.pro , resolution = set.resolutions, verbose = FALSE) 
clustree(T12.pro)
T12.pro  <- RunUMAP(T12.pro , dims = 1:50)
# T12.res <- sapply(set.resolutions, function(x){
#   p <- DimPlot(object = T12.pro, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
#   print(p)
# })

dev.off()

#### remove doublet
library(DoubletFinder) # Require cleanup of low-quality cells in advance
source(file = "scripts/SingleCell/doubletDetect.R")
pdf("figure/SingleCell/1.QualityControl/doublet.pdf")

T1.pro1 <- doubletDetect(Seurat.object = T1.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.8", sct = T)
T2.pro1 <- doubletDetect(Seurat.object = T2.pro, PCs = 1:50, doublet.rate = 0.106, annotation = "SCT_snn_res.1.9", sct = T)
T3.pro1 <- doubletDetect(Seurat.object = T3.pro, PCs = 1:50, doublet.rate = 0.091, annotation = "SCT_snn_res.0.8", sct = T)
T4.pro1 <- doubletDetect(Seurat.object = T4.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.5", sct = T)
T5.pro1 <- doubletDetect(Seurat.object = T5.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.6", sct = T)
T6.pro1 <- doubletDetect(Seurat.object = T6.pro, PCs = 1:50, doublet.rate = 0.106, annotation = "SCT_snn_res.0.5", sct = T)
T7.pro1 <- doubletDetect(Seurat.object = T7.pro, PCs = 1:50, doublet.rate = 0.091, annotation = "SCT_snn_res.0.6", sct = T)
T8.pro1 <- doubletDetect(Seurat.object = T8.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.6", sct = T)
T9.pro1 <- doubletDetect(Seurat.object = T9.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.9", sct = T)
T10.pro1 <- doubletDetect(Seurat.object = T10.pro, PCs = 1:50, doublet.rate = 0.106, annotation = "SCT_snn_res.0.8", sct = T)
T11.pro1 <- doubletDetect(Seurat.object = T11.pro, PCs = 1:50, doublet.rate = 0.091, annotation = "SCT_snn_res.2", sct = T)
T12.pro1 <- doubletDetect(Seurat.object = T12.pro, PCs = 1:50, doublet.rate = 0.061, annotation = "SCT_snn_res.0.7", sct = T)
dev.off()
saveRDS(T1.pro1, file = "Rdata/SingleCell/T1.pro1.rds")
saveRDS(T2.pro1, file = "Rdata/SingleCell/T2.pro1.rds")
saveRDS(T3.pro1, file = "Rdata/SingleCell/T3.pro1.rds")
saveRDS(T4.pro1, file = "Rdata/SingleCell/T4.pro1.rds")
saveRDS(T5.pro1, file = "Rdata/SingleCell/T5.pro1.rds")
saveRDS(T6.pro1, file = "Rdata/SingleCell/T6.pro1.rds")
saveRDS(T7.pro1, file = "Rdata/SingleCell/T7.pro1.rds")
saveRDS(T8.pro1, file = "Rdata/SingleCell/T8.pro1.rds")
saveRDS(T9.pro1, file = "Rdata/SingleCell/T9.pro1.rds")
saveRDS(T10.pro1, file = "Rdata/SingleCell/T10.pro1.rds")
saveRDS(T11.pro1, file = "Rdata/SingleCell/T11.pro1.rds")
saveRDS(T12.pro1, file = "Rdata/SingleCell/T12.pro1.rds")

pdf("figure/SingleCell/1.QualityControl/doublet.cell.pdf")
DimPlot(object = T1.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T2.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T3.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T4.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T5.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T6.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T7.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T8.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T9.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T10.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T11.pro1, reduction = 'umap', group.by = "Doublet")
DimPlot(object = T12.pro1, reduction = 'umap', group.by = "Doublet")
dev.off()

T1.pro2 <- subset(T1.pro1, subset = Doublet == "Singlet")
T2.pro2 <- subset(T2.pro1, subset = Doublet == "Singlet")
T3.pro2 <- subset(T3.pro1, subset = Doublet == "Singlet")
T4.pro2 <- subset(T4.pro1, subset = Doublet == "Singlet")
T5.pro2 <- subset(T5.pro1, subset = Doublet == "Singlet")
T6.pro2 <- subset(T6.pro1, subset = Doublet == "Singlet")
T7.pro2 <- subset(T7.pro1, subset = Doublet == "Singlet")
T8.pro2 <- subset(T8.pro1, subset = Doublet == "Singlet")
T9.pro2 <- subset(T9.pro1, subset = Doublet == "Singlet")
T10.pro2 <- subset(T10.pro1, subset = Doublet == "Singlet")
T11.pro2 <- subset(T11.pro1, subset = Doublet == "Singlet")
T12.pro2 <- subset(T12.pro1, subset = Doublet == "Singlet")
saveRDS(T1.pro2, file = "Rdata/SingleCell/T1.pro2.rds")
saveRDS(T2.pro2, file = "Rdata/SingleCell/T2.pro2.rds")
saveRDS(T3.pro2, file = "Rdata/SingleCell/T3.pro2.rds")
saveRDS(T4.pro2, file = "Rdata/SingleCell/T4.pro2.rds")
saveRDS(T5.pro2, file = "Rdata/SingleCell/T5.pro2.rds")
saveRDS(T6.pro2, file = "Rdata/SingleCell/T6.pro2.rds")
saveRDS(T7.pro2, file = "Rdata/SingleCell/T7.pro2.rds")
saveRDS(T8.pro2, file = "Rdata/SingleCell/T8.pro2.rds")
saveRDS(T9.pro2, file = "Rdata/SingleCell/T9.pro2.rds")
saveRDS(T10.pro2, file = "Rdata/SingleCell/T10.pro2.rds")
saveRDS(T11.pro2, file = "Rdata/SingleCell/T11.pro2.rds")
saveRDS(T12.pro2, file = "Rdata/SingleCell/T12.pro2.rds")

############################################## merge data and correct the batch effect
DefaultAssay(T1.pro2) <- "RNA"
DefaultAssay(T2.pro2) <- "RNA"
DefaultAssay(T3.pro2) <- "RNA"
DefaultAssay(T4.pro2) <- "RNA"
DefaultAssay(T5.pro2) <- "RNA"
DefaultAssay(T6.pro2) <- "RNA"
DefaultAssay(T7.pro2) <- "RNA"
DefaultAssay(T8.pro2) <- "RNA"
DefaultAssay(T9.pro2) <- "RNA"
DefaultAssay(T10.pro2) <- "RNA"
DefaultAssay(T11.pro2) <- "RNA"
DefaultAssay(T12.pro2) <- "RNA"

source(file = "scripts/SingleCell/variableFeatureSelection.R")
CRC.list <- list(T1 = T1.pro2, T2 = T2.pro2, T3 = T3.pro2, T4 = T4.pro2,
                   T5 = T5.pro2, T6 = T6.pro2, T7 = T7.pro2, T8 = T8.pro2,
                   T9 = T9.pro2, T10 = T10.pro2, T11 = T11.pro2, T12 = T12.pro2)
CRC.list.Stardard <- variableFeatureSelection(seurat.lists = CRC.list, method = "Stardard", nfeatures = 3000)
saveRDS(CRC.list.Stardard, file = "Rdata/SingleCell/CRC.list.Stardard.3000.rds")

CRC.list.SCT <- variableFeatureSelection(seurat.lists = CRC.list, method = "SCT", nfeatures = 3000)
saveRDS(CRC.list.SCT, file = "Rdata/SingleCell/CRC.list.SCT.3000.rds")

#### 
# assay=SCT
data.merge <- merge(CRC.list.SCT[[1]], y = CRC.list.SCT[2:length(CRC.list.SCT)], project = "CRC")
DefaultAssay(data.merge) <- "SCT"
seurat.features.SCT <- SelectIntegrationFeatures(object.list = CRC.list.SCT, nfeatures = 3000)
VariableFeatures(data.merge) <- seurat.features.SCT
# Remove previous clustering results
index <- match(paste0("SCT_snn_res.", seq(0.5, 2, by=0.1)), colnames(data.merge@meta.data))
data.merge@meta.data <- data.merge@meta.data[,-index]
# assay=RNA
seurat.features.RNA <- SelectIntegrationFeatures(object.list = CRC.list.Stardard, nfeatures = 3000)
DefaultAssay(data.merge) <- "RNA"
VariableFeatures(data.merge) <- seurat.features.RNA
data.merge <- NormalizeData(data.merge, verbose = FALSE)
data.merge <- ScaleData(data.merge, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mt"), features = rownames(data.merge@assays$RNA@data))
DefaultAssay(data.merge) <- "SCT"
saveRDS(data.merge, file = "Rdata/SingleCell/data.merge.rds")

##############################################2.Evaluation of patient bias
DefaultAssay(data.merge) <- "SCT"
pdf("figure/SingleCell/1.QualityControl/filtered.statistics.pdf")
VlnPlot(object = data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", 
        # cols = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))], 
        pt.size = 0)
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
FeatureScatter(object = data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
dev.off()
# Draw the distribution of the number of samples
cell.number <- as.data.frame(table(data.merge$orig.ident))
pdf("figure/SingleCell/1.QualityControl/highQuality.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", 
          # palette = Palettes$group_pal[1:length(unique(data.merge@meta.data$orig.ident))],
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()

# set.resolutions <- seq(0.2, 1.2, by = 0.1)
set.resolutions <- 0.7
#### First observe whether the clustering effect will depend on the sample
a <- data.merge
a <- RunPCA(a, npcs = 100, verbose = T)
pdf("figure/SingleCell/1.QualityControl/merge.observe.batch.pdf")
ElbowPlot(object = a, ndims = 100)
a <- FindNeighbors(a, dims = 1:50, verbose = T)
a <- FindClusters(object = a, resolution = set.resolutions, verbose = T) 
clustree(a)
a <- RunUMAP(a, dims = 1:50)
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = "orig.ident")
DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = "Phase")
merge.res <- sapply(set.resolutions, function(x){
  p <- DimPlot(object = a, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
  print(p)
})
dev.off()

##############################################3.correct batch effect
source(file = "scripts/SingleCell/scRNA.Integrate.multipleSample.R")
pdf("figure/SingleCell/2.Cluster/SCT.Harmony.Integration.PC40.feature3000.pdf")
data.merge.harmony.PC40.SCT <- Harmony.integration.reduceDimension(seurat.object = data.merge, assay = "SCT",set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 40, nfeatures = 3000, npcs = 50)
dev.off()
saveRDS(data.merge.harmony.PC40.SCT, file = "Rdata/SingleCell/data.merge.harmony.PC40.SCT.feature3000.rds")
