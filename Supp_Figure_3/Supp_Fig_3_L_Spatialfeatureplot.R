#script to generate Supplemental Figure 3L. 
#Spatial feature plot assessing gene expression signature 31 days following PR8
#load in relavant libraries and data

C1_31 <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220822_spatial_31dpi/C1_31.rds")

#make sure correct data was loaded in 
SpatialDimPlot(C1_31)

#generate feature plot for selected genes
SpatialFeaturePlot(C1_31, features = c("Krt5", "Sftpc", "Scgb3a2", "Scgb1a1"))
