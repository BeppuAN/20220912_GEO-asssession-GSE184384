#code used to annotate epithelial dataset
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated <- readRDS("/Volumes/20190816.singlecellRNAseq/20190710.secondary anaysis.epi/20190826.datasetintegration.withdemultiplexing/epitimecourse.integrated.rds")
x <- epitimecourse.integrated

#dataset contains clusters populated by Stromal/Immune/Endothelial cells which aren't of interest to this study
#These cells will be removed from the dataset
FeaturePlot(x, features = c("rna_Vim", "rna_Ptprc", "rna_Acta2")) 
Idents(x) <- "orig.ident"
plot <- DimPlot(x, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(x, cells = select.cells) <- "Stromal_Immune_Endo"
x <- subset(x, idents = "SeuratProject")

#recluster
x <- FindNeighbors(x, dims = 1:30)
x <- FindClusters(x, resolution = 0.5)
x <- RunUMAP(x, dims = 1:30)

#Perform Differential expression analysis
markers <- FindAllMarkers(x)

#group clusters with similar gene expression profiles
x <- RenameIdents(x, "1" = "TII")
x <- RenameIdents(x, "0" = "TII")
x <- RenameIdents(x, "2" = "TII")
x <- RenameIdents(x, "3" = "TII")
x <- RenameIdents(x, "9" = "TII")
x <- RenameIdents(x, "8" = "Ciliated.1")
x <- RenameIdents(x, "4" = "Ciliated.1")
x <- RenameIdents(x, "12" = "Ciliated.2")
x <- RenameIdents(x, "11" = "TI")
x <- RenameIdents(x, "10" = "Trans.")
x <- RenameIdents(x, "7" = "BC")
x <- RenameIdents(x, "6" = "Serous")
x <- RenameIdents(x, "5" = "Club")

#stash Idents
x@meta.data$anno <- Idents(x)

saveRDS(x, file = "epitimecourse.integrated.V2.rds")

