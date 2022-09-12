#script for Supplemental Fig. 2B
#Spatial feature plots of PR8 infected mouse lung (14days following PR8)
#load in relevant datasets and libraries
library(Seurat)
C1_14 <-  readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220822_spatial_31dpi/influenza_spatial.rds")

DefaultAssay(C1_14) <- 'Spatial'
SpatialDimPlot(C1_14)
SpatialFeaturePlot(C1_14, features = "Vim") #fibroblast
SpatialFeaturePlot(C1_14, features = "Ptprc") #immune
SpatialFeaturePlot(C1_14, features = "Pecam1") #smooth muscle
SpatialFeaturePlot(C1_14, features = "Hopx") #type I
SpatialFeaturePlot(C1_14, features = "Sftpc") #type II
SpatialFeaturePlot(C1_14, features = "Krt5") #BC
SpatialFeaturePlot(C1_14, features = "Scgb1a1") #Club
SpatialFeaturePlot(C1_14, features = "Scgb3a2") #Pan-secretory
SpatialFeaturePlot(C1_14, features = "Bpifa1") #Serous
SpatialFeaturePlot(C1_14, features = "Muc5ac") #goblet
SpatialFeaturePlot(C1_14, features = "Dclk1") #tuft

