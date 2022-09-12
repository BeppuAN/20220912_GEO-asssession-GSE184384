#script for Supplemental Fig. 3G
#Expression of H2-K1 within epithelial cells at indicated timepoint following PR8 exposure.
#load in relevent datasets and libraries

influenza_timecourse_epithelial <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220822_spatial_31dpi/influenza_timecourse_epithelial.rds")

table(Idents(influenza_timecourse_epithelial))
IS <- subset(influenza_timecourse_epithelial, idents = "Serous")
VlnPlot(IS, features = "rna_H2-K1", group.by = "timepoint")
