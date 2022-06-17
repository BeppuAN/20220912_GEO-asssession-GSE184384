#we will be perform dataset integration of all epithelial datasets
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)

#read in singlecell data containing hashed samples
#Ammend batch and condition info
flugenes.n357.singlet <- readRDS("~/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/flugenes.n357.singlet.rds")
flugenes.n357.singlet@meta.data
flugenes.n357.singlet$tech <- "batch 1"

flugenes.91114.singlet <- readRDS("~/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/flugenes.91114.singlet.rds")
flugenes.91114.singlet@meta.data
flugenes.91114.singlet$tech <- "batch 2"

fourmpi.trachea.singlet <- readRDS("~/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/fourmpi.trachea.singlet.rds")  
fourmpi.trachea.singlet@meta.data
fourmpi.trachea.singlet$tech <- "batch 3"

#read in singlecell data containing non-hashed samples
#Ammend batch and condition info
redo17dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/redo-17dpi-noflu/filtered_gene_bc_matrices/mm10")
redo17dpiflugenes <- CreateSeuratObject(counts = redo17dpiflugenes)
redo17dpiflugenes <- NormalizeData(object = redo17dpiflugenes)
redo17dpiflugenes <- FindVariableFeatures(object = redo17dpiflugenes)
redo17dpiflugenes <- ScaleData(object = redo17dpiflugenes)
redo17dpiflugenes <- RunPCA(object = redo17dpiflugenes)
redo17dpiflugenes <- FindNeighbors(object = redo17dpiflugenes)
redo17dpiflugenes <- FindClusters(object = redo17dpiflugenes)
redo17dpiflugenes <- RunTSNE(object = redo17dpiflugenes)
DimPlot(object = redo17dpiflugenes, reduction = "tsne")
saveRDS(object = redo17dpiflugenes, file = "redo17dpiflugenes.rds")
redo17dpiflugenes@meta.data
redo17dpiflugenes$tech <- "batch 4"
redo17dpiflugenes$timepoint <- "17 dpi"

redo21dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/redo-21dpi-flugenes/filtered_gene_bc_matrices/mm10")
redo21dpiflugenes <- CreateSeuratObject(counts = redo21dpiflugenes)
redo21dpiflugenes <- NormalizeData(object = redo21dpiflugenes)
redo21dpiflugenes <- FindVariableFeatures(object = redo21dpiflugenes)
redo21dpiflugenes <- ScaleData(object = redo21dpiflugenes)
redo21dpiflugenes <- RunPCA(object = redo21dpiflugenes)
redo21dpiflugenes <- FindNeighbors(object = redo21dpiflugenes)
redo21dpiflugenes <- FindClusters(object = redo21dpiflugenes)
redo21dpiflugenes <- RunTSNE(object = redo21dpiflugenes)
DimPlot(object = redo21dpiflugenes, reduction = "tsne")
saveRDS(object = redo21dpiflugenes, file = "redo21dpiflugenes.rds")
redo21dpiflugenes@meta.data
redo21dpiflugenes$tech <- "batch 5"
redo21dpiflugenes$timepoint <- "21 dpi"

redo8mpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/20190823.epi-8mpi-flugenes/filtered_feature_bc_matrix/mm10")
redo8mpiflugenes <- CreateSeuratObject(counts = redo8mpiflugenes)
redo8mpiflugenes <- NormalizeData(object = redo8mpiflugenes)
redo8mpiflugenes <- FindVariableFeatures(object = redo8mpiflugenes)
redo8mpiflugenes <- ScaleData(object = redo8mpiflugenes)
redo8mpiflugenes <- RunPCA(object = redo8mpiflugenes)
redo8mpiflugenes <- FindNeighbors(object = redo8mpiflugenes)
redo8mpiflugenes <- FindClusters(object = redo8mpiflugenes)
redo8mpiflugenes <- RunTSNE(object = redo8mpiflugenes)
DimPlot(object = redo8mpiflugenes, reduction = "tsne")
saveRDS(object = redo8mpiflugenes, file = "redo8mpiflugenes.rds")
redo8mpiflugenes@meta.data
redo8mpiflugenes$tech <- "batch 6"
redo8mpiflugenes$timepoint <- "8 mpi"

redotracheaflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/20190823.trachea-flugenes/filtered_feature_bc_matrix/mm10")
redotracheaflugenes <- CreateSeuratObject(counts = redotracheaflugenes)
redotracheaflugenes <- NormalizeData(object = redotracheaflugenes)
redotracheaflugenes <- FindVariableFeatures(object = redotracheaflugenes)
redotracheaflugenes <- ScaleData(object = redotracheaflugenes)
redotracheaflugenes <- RunPCA(object = redotracheaflugenes)
redotracheaflugenes <- FindNeighbors(object = redotracheaflugenes)
redotracheaflugenes <- FindClusters(object = redotracheaflugenes)
redotracheaflugenes <- RunTSNE(object = redotracheaflugenes)
DimPlot(object = redotracheaflugenes, reduction = "tsne")
saveRDS(object = redotracheaflugenes, file = "redotracheaflugenes.rds")
redotracheaflugenes@meta.data
redotracheaflugenes$tech <- "batch 6"
redotracheaflugenes$timepoint <- "trachea"

oldNdpiflugenes <-Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-nepi-flugenes/filtered_gene_bc_matrices/mm10")
oldNdpiflugenes <- CreateSeuratObject(counts = oldNdpiflugenes)
oldNdpiflugenes <- NormalizeData(object = oldNdpiflugenes)
oldNdpiflugenes <- FindVariableFeatures(object = oldNdpiflugenes)
oldNdpiflugenes <- ScaleData(object = oldNdpiflugenes)
oldNdpiflugenes <- RunPCA(object = oldNdpiflugenes)
oldNdpiflugenes <- FindNeighbors(object = oldNdpiflugenes)
oldNdpiflugenes <- FindClusters(object = oldNdpiflugenes)
oldNdpiflugenes <- RunTSNE(object = oldNdpiflugenes)
DimPlot(object = oldNdpiflugenes, reduction = "tsne")
saveRDS(object = oldNdpiflugenes, file = "oldNdpiflugenes.rds")
oldNdpiflugenes.sub@meta.data
oldNdpiflugenes.sub$tech <- "batch 7"
oldNdpiflugenes.sub$timepoint <- "naive"

old3dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-3dpi-flugenes/filtered_gene_bc_matrices/mm10")
old3dpiflugenes <- CreateSeuratObject(counts = old3dpiflugenes)
old3dpiflugenes <- NormalizeData(object = old3dpiflugenes)
old3dpiflugenes <- FindVariableFeatures(object = old3dpiflugenes)
old3dpiflugenes <- ScaleData(object = old3dpiflugenes)
old3dpiflugenes <- RunPCA(object = old3dpiflugenes)
old3dpiflugenes <- FindNeighbors(object = old3dpiflugenes)
old3dpiflugenes <- FindClusters(object = old3dpiflugenes)
old3dpiflugenes <- RunTSNE(object = old3dpiflugenes)
DimPlot(object = old3dpiflugenes, reduction = "tsne")
saveRDS(object = old3dpiflugenes, file = "old3dpiflugenes.rds")
old3dpiflugenes.sub@meta.data
old3dpiflugenes.sub$tech <- "batch 7"
old3dpiflugenes.sub$timepoint <- "3 dpi"

old5dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-5dpi-flugenes/filtered_gene_bc_matrices/mm10")
old5dpiflugenes <- CreateSeuratObject(counts = old5dpiflugenes)
old5dpiflugenes <- NormalizeData(object = old5dpiflugenes)
old5dpiflugenes <- FindVariableFeatures(object = old5dpiflugenes)
old5dpiflugenes <- ScaleData(object = old5dpiflugenes)
old5dpiflugenes <- RunPCA(object = old5dpiflugenes)
old5dpiflugenes <- FindNeighbors(object = old5dpiflugenes)
old5dpiflugenes <- FindClusters(object = old5dpiflugenes)
old5dpiflugenes <- RunTSNE(object = old5dpiflugenes)
DimPlot(object = old5dpiflugenes, reduction = "tsne")
saveRDS(object = old5dpiflugenes, file = "old5dpiflugenes.rds")
old5dpiflugenes.sub@meta.data
old5dpiflugenes.sub$tech <- "batch 8"
old5dpiflugenes.sub$timepoint <- "5 dpi"

old7dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-7dpi-flugenes/filtered_gene_bc_matrices/mm10")
old7dpiflugenes <- CreateSeuratObject(counts = old7dpiflugenes)
old7dpiflugenes <- NormalizeData(object = old7dpiflugenes)
old7dpiflugenes <- FindVariableFeatures(object = old7dpiflugenes)
old7dpiflugenes <- ScaleData(object = old7dpiflugenes)
old7dpiflugenes <- RunPCA(object = old7dpiflugenes)
old7dpiflugenes <- FindNeighbors(object = old7dpiflugenes)
old7dpiflugenes <- FindClusters(object = old7dpiflugenes)
old7dpiflugenes <- RunTSNE(object = old7dpiflugenes)
DimPlot(object = old7dpiflugenes, reduction = "tsne")
saveRDS(object = old7dpiflugenes, file = "old7dpiflugenes.rds")
FeaturePlot(old7dpiflugenes, features = "Krt5")
old7dpiflugenes.sub@meta.data
old7dpiflugenes.sub$tech <- "batch 8"
old7dpiflugenes.sub$timepoint <- "7 dpi"

old9dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-9dpi-flugenes/filtered_gene_bc_matrices/mm10")
old9dpiflugenes <- CreateSeuratObject(counts = old9dpiflugenes)
old9dpiflugenes <- NormalizeData(object = old9dpiflugenes)
old9dpiflugenes <- FindVariableFeatures(object = old9dpiflugenes)
old9dpiflugenes <- ScaleData(object = old9dpiflugenes)
old9dpiflugenes <- RunPCA(object = old9dpiflugenes)
old9dpiflugenes <- FindNeighbors(object = old9dpiflugenes)
old9dpiflugenes <- FindClusters(object = old9dpiflugenes)
old9dpiflugenes <- RunTSNE(object = old9dpiflugenes)
DimPlot(object = old9dpiflugenes, reduction = "tsne")
saveRDS(object = old9dpiflugenes, file = "old9dpiflugenes.rds")
old9dpiflugenes.sub@meta.data
old9dpiflugenes.sub$tech <- "batch 9"
old9dpiflugenes.sub$timepoint <- "9 dpi"

old11dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-11dpi-flugenes/filtered_gene_bc_matrices/mm10")
old11dpiflugenes <- CreateSeuratObject(counts = old11dpiflugenes)
old11dpiflugenes <- NormalizeData(object = old11dpiflugenes)
old11dpiflugenes <- FindVariableFeatures(object = old11dpiflugenes)
old11dpiflugenes <- ScaleData(object = old11dpiflugenes)
old11dpiflugenes <- RunPCA(object = old11dpiflugenes)
old11dpiflugenes <- FindNeighbors(object = old11dpiflugenes)
old11dpiflugenes <- FindClusters(object = old11dpiflugenes)
old11dpiflugenes <- RunTSNE(object = old11dpiflugenes)
DimPlot(object = old11dpiflugenes, reduction = "tsne")
saveRDS(object = old11dpiflugenes, file = "old11dpiflugenes.rds")
old11dpiflugenes.sub@meta.data
old11dpiflugenes.sub$tech <- "batch 9"
old11dpiflugenes.sub$timepoint <- "11 dpi"

old14dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-14dpi-flugenes/filtered_gene_bc_matrices/mm10")
old14dpiflugenes <- CreateSeuratObject(counts = old14dpiflugenes)
old14dpiflugenes <- NormalizeData(object = old14dpiflugenes)
old14dpiflugenes <- FindVariableFeatures(object = old14dpiflugenes)
old14dpiflugenes <- ScaleData(object = old14dpiflugenes)
old14dpiflugenes <- RunPCA(object = old14dpiflugenes)
old14dpiflugenes <- FindNeighbors(object = old14dpiflugenes)
old14dpiflugenes <- FindClusters(object = old14dpiflugenes)
old14dpiflugenes <- RunTSNE(object = old14dpiflugenes)
DimPlot(object = old14dpiflugenes, reduction = "tsne")
saveRDS(object = old14dpiflugenes, file = "old14dpiflugenes.rds")
old14dpiflugenes.sub@meta.data
old14dpiflugenes.sub$tech <- "batch 10"
old14dpiflugenes.sub$timepoint <- "14 dpi"

old17dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-17dpi-flugenes /filtered_gene_bc_matrices/mm10")
old17dpiflugenes <- CreateSeuratObject(counts = old17dpiflugenes)
old17dpiflugenes <- NormalizeData(object = old17dpiflugenes)
old17dpiflugenes <- FindVariableFeatures(object = old17dpiflugenes)
old17dpiflugenes <- ScaleData(object = old17dpiflugenes)
old17dpiflugenes <- RunPCA(object = old17dpiflugenes)
old17dpiflugenes <- FindNeighbors(object = old17dpiflugenes)
old17dpiflugenes <- FindClusters(object = old17dpiflugenes)
old17dpiflugenes <- RunTSNE(object = old17dpiflugenes)
DimPlot(object = old17dpiflugenes, reduction = "tsne")
saveRDS(object = old17dpiflugenes, file = "old17dpiflugenes.rds")
old17dpiflugenes.sub@meta.data
old17dpiflugenes.sub$tech <- "batch 10"
old17dpiflugenes.sub$timepoint <- "17 dpi"

old21dpiflugenes <- Read10X("/Users/beppuan/Desktop/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/old-21dpi-flugenes/filtered_gene_bc_matrices/mm10")
old21dpiflugenes <- CreateSeuratObject(counts = old21dpiflugenes)
old21dpiflugenes <- NormalizeData(object = old21dpiflugenes)
old21dpiflugenes <- FindVariableFeatures(object = old21dpiflugenes)
old21dpiflugenes <- ScaleData(object = old21dpiflugenes)
old21dpiflugenes <- RunPCA(object = old21dpiflugenes)
old21dpiflugenes <- FindNeighbors(object = old21dpiflugenes)
old21dpiflugenes <- FindClusters(object = old21dpiflugenes)
old21dpiflugenes <- RunTSNE(object = old21dpiflugenes)
DimPlot(object = old21dpiflugenes, reduction = "tsne")
saveRDS(object = old21dpiflugenes, file = "old21dpiflugenes.rds")
old21dpiflugenes.sub@meta.data
old21dpiflugenes.sub$tech <- "batch 10"
old21dpiflugenes.sub$timepoint <- "21 dpi"

#double check meta data of all datasets as a sanity check
old21dpiflugenes.sub@meta.data
old17dpiflugenes.sub@meta.data
old14dpiflugenes.sub@meta.data
old11dpiflugenes.sub@meta.data
old9dpiflugenes.sub@meta.data
old7dpiflugenes.sub@meta.data
old5dpiflugenes.sub@meta.data
old3dpiflugenes.sub@meta.data
oldNdpiflugenes.sub@meta.data
redotracheaflugenes@meta.data
redo8mpiflugenes@meta.data
redo21dpiflugenes@meta.data
redo17dpiflugenes@meta.data
fourmpi.trachea.singlet@meta.data
flugenes.91114.singlet@meta.data
flugenes.n357.singlet@meta.data

#next step would be to merge the datasets
epitimecourse.preintegration <- merge(x = flugenes.n357.singlet, y = c(flugenes.91114.singlet, fourmpi.trachea.singlet, redo17dpiflugenes, redo21dpiflugenes, redo8mpiflugenes, redotracheaflugenes, oldNdpiflugenes.sub, old3dpiflugenes.sub, old5dpiflugenes.sub, old7dpiflugenes.sub, old9dpiflugenes.sub, old11dpiflugenes.sub, old14dpiflugenes.sub, old17dpiflugenes.sub, old21dpiflugenes.sub))
epitimecourse.preintegration@meta.data

#remove mitochondrial + non UMI
epitimecourse.preintegration[["percent.mt"]] <- PercentageFeatureSet(epitimecourse.preintegration, pattern = "^mt-")
VlnPlot(epitimecourse.preintegration, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
epitimecourse.preintegration <- subset(epitimecourse.preintegration, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

#save before integration
saveRDS(object = epitimecourse.preintegration, file = "flugenes.timecourse.preintegration.rds")

#now work through the data integration pipeline
epitimecourse.preintegration.list <- SplitObject(epitimecourse.preintegration, split.by = "tech")

for (i in 1:length(epitimecourse.preintegration.list)) {
  epitimecourse.preintegration.list[[i]] <- NormalizeData(epitimecourse.preintegration.list[[i]], verbose = FALSE)
  epitimecourse.preintegration.list[[i]] <- FindVariableFeatures(epitimecourse.preintegration.list[[i]], selection.method = "vst", 
                                                                 nfeatures = 6000, verbose = FALSE)
}

reference.list <- epitimecourse.preintegration.list[c("batch 1", "batch 2", "batch 3", "batch 4", "batch 5", "batch 6", "batch 7", "batch 8", "batch 9", "batch 10")]
epitimecourse.preintegration.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, k.filter = 150)

epitimecourse.integrated <- IntegrateData(anchorset = epitimecourse.preintegration.anchors, dims = 1:30)

# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(epitimecourse.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
epitimecourse.integrated <- ScaleData(epitimecourse.integrated, verbose = FALSE)
epitimecourse.integrated <- RunPCA(epitimecourse.integrated, npcs = 30, verbose = FALSE)
epitimecourse.integrated <- RunTSNE(epitimecourse.integrated, reduction = "pca", dims = 1:30)
epitimecourse.integrated <- RunUMAP(epitimecourse.integrated, reduction = "pca", dims = 1:30)
DimPlot(epitimecourse.integrated, group.by = "seurat_clusters")

#save
saveRDS(object = epitimecourse.integrated, file = "epitimecourse.integrated.rds")

