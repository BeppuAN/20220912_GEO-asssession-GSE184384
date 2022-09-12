#script for Supplemental Fig. 6B
#Gene expression of IL-22ra1 
#load in relevant datasets and libraries
influenza_timecourse_epithelial <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220822_spatial_31dpi/influenza_timecourse_epithelial.rds")

DotPlot(influenza_timecourse_epithelial, features = "rna_Il22ra1", cols = c("light blue", "red"))

