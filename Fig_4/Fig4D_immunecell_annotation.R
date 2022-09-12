#Annotation of the immune dataset
#load in data
library(Seurat)
library(ggplot2)
immune.combined.V30 <- readRDS("C:/Users/Stripp Lab/Downloads/20220512.immunecombined/immune.combined.V30.rds")

#make sure correct data set was loaded in
DimPlot(immune.combined.V30, label = T)
DimPlot(immune.combined.V32, label = T)

#Feature plot cell type specific immune markers
FeaturePlot(immune.combined.V30, features = c("rna_Cd3g",
                                              "rna_Cd8a",
                                              "rna_Cd4",
                                              "rna_Top2a",
                                              "rna_Il2ra",
                                              "rna_Il17a",
                                              "rna_Il13", 
                                              "rna_Lyz2",
                                              "rna_Mrc1", 
                                              "rna_C1qc",
                                              "rna_Flt3", 
                                              "rna_S100a9",
                                              "rna_Cd79a", 
                                              "rna_Gzma",
                                              "rna_Foxp3", 
                                              "rna_Ly6c2", 
                                              "rna_Ccr2", 
                                              "rna_Jchain"), order = T)

#apply simple annotation 
immune.combined.V31 <-RenameIdents(immune.combined.V30, "0" = "T_Lymph.1")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "1" = "B_Lymph.1")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "2" = "Myeloid.1")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "3" = "B_Lymph.2")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "4" = "T_Lymph.2")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "5" = "T_Lymph.3")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "6" = "Myeloid.2")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "7" = "Myeloid.3")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "8" = "T_Lymph.4")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "9" = "T_Lymph.5")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "10" = "Myeloid.4")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "11" = "Myeloid.5")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "12" = "Myeloid.6")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "13" = "T_Lymph.6")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "14" = "T_Lymph.7")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "15" = "Myeloid.8")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "16" = "Unknown.1")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "17" = "Myeloid.9")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "18" = "Unknown.2")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "19" = "Unknown.3")
immune.combined.V31 <-RenameIdents(immune.combined.V31, "20" = "Myeloid.10")

#apply more complex annotation based off feature gene expression
immune.combined.V32 <-RenameIdents(immune.combined.V31, "T_Lymph.1" = "Activated_CD8")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.2" = "Activated_CD4")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.3" = "Naive_T_Lymph")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.4" = "NK")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.5" = "Prolif_T_Lymph")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.6" = "gd_T")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "T_Lymph.7" = "Tregs")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.1" = "Mono.1")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.2" = "Myeloid.1")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.3" = "Neut.")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.4" = "Mono.2")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.5" = "Int_Mac")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.6" = "DC.1")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.8" = "DC.2")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Unknown.1" = "DC.3")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.9" = "A_Mac")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Myeloid.10" = "Myeloid.2")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Unknown.3" = "Fibroblast")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Unknown.2" = "Activated_CD8")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Activated_CD8" = "Activated_CD4.")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Activated_CD4" = "Activated_CD8.")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Activated_CD8." = "Activated_CD8")
immune.combined.V32 <-RenameIdents(immune.combined.V32, "Activated_CD4." = "Activated_CD4")


DimPlot(immune.combined.V32, label = T)
#stash annotations into metadata
immune.combined.V32@meta.data$'anno' <- Idents(immune.combined.V32)
table(immune.combined.V32@meta.data$'anno')

# Define an order of cluster identities
my_levels <- c("Naive_T_Lymph","Activated_CD8","Activated_CD4", "Prolif_T_Lymph","gd_T","Tregs","NK","B_Lymph.1","B_Lymph.2","Myeloid.1","Myeloid.2","Mono.1","Mono.2","Int_Mac","A_Mac","DC.1","DC.2","DC.3","Neut.","Unknown.2","Fibroblast")

# Relevel object@ident
immune.combined.V32@meta.data$'anno' <- factor(immune.combined.V32@meta.data$'anno', levels = my_levels)

DimPlot(immune.combined.V32, group.by = 'anno', label = T)

#save
saveRDS(immune.combined.V32, file = "immune_combined_V33.rds")

#end of session