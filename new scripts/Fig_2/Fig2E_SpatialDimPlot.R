#Script to generate UMAP projection for spatial RNAseq data
#load in relavant dataset and libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
C1 <- Load10X_Spatial(data.dir = "/Volumes/20190816.singlecellRNAseq/20200916.visium_analysis/20200910.Novagenerun/space_C1/outs")

#Data preprocessing
plot1 <- VlnPlot(C1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(C1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
C1 <- SCTransform(C1, assay = "Spatial", verbose = FALSE)

#Dim reduction
C1 <- RunPCA(C1, assay = "SCT", verbose = FALSE)
C1 <- FindNeighbors(C1, reduction = "pca", dims = 1:30)
C1 <- FindClusters(C1, verbose = FALSE)
C1 <- RunUMAP(C1, reduction = "pca", dims = 1:30)

#generate Spatial UMAP
SpatialDimPlot(C1)

#save
saveRDS(C1, file = 'C1.rds')

