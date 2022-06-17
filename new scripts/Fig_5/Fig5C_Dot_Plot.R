#Script for generation of DotPlot in Fig 5C
#load in relavant libraries and datasets
library(Seurat)
library(ggplot2)
immune.combined.V32 <- readRDS("C:/Users/Stripp Lab/Downloads/20220512.immunecombined/immune.combined.V32")

#Generate DotPlot
DotPlot(immune.combined.V32, features = c("rna_Cd3g","rna_Cd8a", "rna_Cd4","rna_Top2a","rna_Tcrg-C1","rna_Foxp3","rna_Gzma","rna_Cd79a","rna_Lyz2","rna_Cd14","rna_C1qc","rna_Mrc1","rna_Flt3","rna_S100a9", "rna_Col1a1"), group.by = "anno", cols = c("light blue", "red")) + theme(axis.text.x = element_text(angle = 90))
?DotPlot

#end of session