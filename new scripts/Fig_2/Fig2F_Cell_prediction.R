#script to produce celltype prediction scores for the spatial RNAseq data
#load in relavant dataset and libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
C1 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/C1.rds")
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#use pre-existing singlecell RNAsseq data to generate cell type prediction score for the spatial dataset
x <- SCTransform(x, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
C1.V1 <- SCTransform(C1.V1, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

# the annotation is stored in the 'anno' column of object metadata
DimPlot(x, group.by = "anno", label = TRUE)

anchors <- FindTransferAnchors(reference = x, query = C1.V1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = x$anno, prediction.assay = T,
                                  weight.reduction = C1.V1[["pca"]], dims = 1:10, k.weight = 25)

C1.V1[["predictions"]] <- predictions.assay
DefaultAssay(C1.V1) <- "predictions"
