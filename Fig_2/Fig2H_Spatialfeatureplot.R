#script to produce spatial featureplots
#load in relavant dataset and libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
C1 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/C1.rds")

SpatialFeaturePlot(C1, features = "Krt5")
SpatialFeaturePlot(C1, features = "Cldn10")
SpatialFeaturePlot(C1, features = "Scgb3a2")
SpatialFeaturePlot(C1, features = "Scgb1a1")
