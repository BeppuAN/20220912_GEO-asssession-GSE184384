#code used to visualize transcriptional landscape of lung serous cells
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

DotPlot(x, features = c("rna_Msln", "rna_Bpifa1", "rna_Ltf", "rna_Scgb3a2", "rna_Scgb1a1", "rna_Cldn10"), group.by = "anno")
