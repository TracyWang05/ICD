#### tumor ####
rm(list=ls())
options(stringsAsFactors = F)

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
ES <- enrichIt(obj = expr, gene.sets = geneset)
score=ES[,1]
data@meta.data$score <-score
data@meta.data$ICDgroup=ifelse(data@meta.data$score < median(score),'low-ICD','high-ICD')
diff_tumor <- FindMarkers(data,group.by = "ICDgroup",
                               ident.1 ="high-ICD",
                               ident.2="low-ICD") 
sig_diff_tumor <- subset(diff_tumor, p_val_adj<0.05&abs(avg_log2FC)>1) # PI3
saveRDS(data, "Rdata/SingleCell/tumor-ICD.rds")

#### T score ####
rm(list=ls())
options(stringsAsFactors = F)

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
ES <- enrichIt(obj = expr, gene.sets = geneset)
score=ES[,1]
Tcell@meta.data$score <-score
Tcell@meta.data$ICDgroup=ifelse(Tcell@meta.data$score < median(score),'low-ICD','high-ICD')
diff_T <- FindMarkers(Tcell,group.by = "ICDgroup",
                          ident.1 ="high-ICD",
                          ident.2="low-ICD") 
sig_diff_T <- subset(diff_T, p_val_adj<0.05&abs(avg_log2FC)>1) # TNF
saveRDS(Tcell, "Rdata/SingleCell/T-ICD.rds")
