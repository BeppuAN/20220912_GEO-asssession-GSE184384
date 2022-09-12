#code used to assess gene epression of BC genes in candidate non-BC progenitors in response to flu infection
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#Subset cells of interest
Idents(x) <- "anno"
table(Idents(x))
Serous <- subset(x, idents = c("Serous"))
Club <- subset(x, idents = c("Club"))
TII <- subset(x, idents = c("TII"))

#Generate Dotplots to assess gene expression as a function of "days post infection"
DotPlot(Serous, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Bpifa1", "rna_Ltf", "rna_Scgb3a2"), group.by = "timepoint", scale.max = 10)
DotPlot(Club, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Hp", "rna_Scgb1a1", "rna_Cldn10"), group.by = "timepoint", scale.max = 10)
DotPlot(TII, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Cd74", "rna_Chil1", "rna_Sftpc"), group.by = "timepoint", scale.max = 10)
