#script used to generate Fig. 5D - Cell prediction of immune cells onto visium spot clusters
#load in relavant libraries and datasets
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
C1.V2 <- readRDS("C:/Users/Stripp Lab 2/Downloads/20220608.Cell_prediction_troubleshoot/20220608.Cell_prediction_troubleshoot/C1.V2.rds")
immune_combined_20220607 <- readRDS("C:/Users/Stripp Lab 2/Downloads/20220608.Cell_prediction_troubleshoot/20220608.Cell_prediction_troubleshoot/immune_combined_20220607.rds")
immune_combined_V34 <- immune_combined_20220607

#generate reference based on pure immune populations and timepoints n-21 days post infection
immune_combined_V35 <- subset(immune_combined_V34, idents = c("A_Mac", "Activated_CD4","Activated_CD8","B_Lymph.1","B_Lymph.2","DC.1","DC.2","DC.3","gd_T","Int_Mac","Mono.1","Mono.2","Myeloid.1","Myeloid.2","Naive_T_Lymph","Neut.","NK","Prolif_T_Lymph","Tregs"))
immune_combined_V35@meta.data$anno <- Idents(immune_combined_V35)
Idents(immune_combined_V35) <- "condition"
table(Idents(immune_combined_V35))
immune_combined_V37 <- subset(immune_combined_V35, idents = c("n","3","5","7","9","11","14","17","21"))
Idents(immune_combined_V35) <- 'anno'

#stash Ident into metadata
DimPlot(immune_combined_V36) 
immune_combined_V36@meta.data$anno <- Idents(immune_combined_V36)
table(immune_combined_V36@meta.data$anno)

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# After subsetting, we renormalize cortex
immune_combined_V37 <- SCTransform(immune_combined_V37, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

C1.V2 <- SCTransform(C1.V2, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

#transfer anchors
anchors <- FindTransferAnchors(reference = immune_combined_V37, query = C1.V2, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = immune_combined_V37$anno, prediction.assay = TRUE,
                                  weight.reduction = C1.V2[["pca"]], dims = 1:30, k.weight = 50)

C1.V2[["predictions"]] <- predictions.assay
table(C1.V2[["predictions"]])
DefaultAssay(C1.V2) <- "predictions"

#visualize
SpatialFeaturePlot(C1.V2, features = "gd-T", pt.size.factor = 2)

#save
saveRDS(C1.V2, file = "C1.V5.rds")



