#script for Supplemental Fig. 3F
#Expression of H2-K1 within the epithelial scRNAseq dataset
#load in relavent datasets and libraries

influenza_timecourse_epithelial <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220822_spatial_31dpi/influenza_timecourse_epithelial.rds")

FeaturePlot(influenza_timecourse_epithelial, features = "rna_H2-K1")
