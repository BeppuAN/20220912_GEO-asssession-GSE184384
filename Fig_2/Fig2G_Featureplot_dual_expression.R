#Script to generate Feature co-expression
#load in relavant dataset and libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

C1 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220616.scRNAseq_Scripts/C1.rds")

FeaturePlot(C1, feature = c("Krt5","Scgb3a2"), blend = T, cols = c("black", "green","red"))



