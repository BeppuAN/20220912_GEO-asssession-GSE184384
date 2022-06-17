#script to generate Fig. 5E - Dot plot of type 17 activation in gamma/delta T cells
#load in relavant libraries and datasets
library(Seurat)
library(ggplot2)

immune_combined_V33 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220519.Flu_manuscript_figure_5_redo/immune_combined_V33.rds")

DotPlot(immune_combined_V33, features = c("rna_Tcrg-C1", "rna_Il17a", "rna_Il23r", "rna_Il22"), cols = c("light blue", "red"))

#end session

