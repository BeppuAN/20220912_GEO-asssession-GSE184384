#Script to generate Fig. 5F - Violin plot of type 17 module score
#Load in relavant libraries and datasets
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
immune_combined_20220607 <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220601_cellprediction/immune_combined_20220607.rds")
influenza_timecourse_epithelial <- readRDS("C:/Users/stripplab/Desktop/20220601_analysis_quants_figures/20220609.Cell_Prediction/influenza_timecourse_epithelial.rds")

# Define an order of cluster identities
my_levels <- c("n", "3","5","7","9","11","14","17","21","120")

# Relevel object@ident
immune_combined_20220607@meta.data$condition <- factor(immune_combined_20220607@meta.data$condition, levels = my_levels)

#subset cells of interest
Basal <- subset(influenza_timecourse_epithelial, idents = "Basal")
gd_t <- subset(immune_combined_20220607, idents = "gd_T")
Neut <- subset(immune_combined_20220607, idents = "Neut.")
Activated_CD4 <- subset(immune_combined_20220607, idents = "Activated_CD4")

#check cytokine expression in a cell type specific manner
#add module score
#generate type17 module score
DefaultAssay(gd_t) <- 'RNA'

type17 <- list(c(
  'Il17a',
  'Il22',
  'Il23r'
))

gd_t <- AddModuleScore(
  object = gd_t,
  features = type17,
  ctrl = 10,
  name = 'type17'
)

#visiualize type_17 response
library(ggplot2)
install.packages('ggsignif')
library(ggsignif)

#with wilcoxin T test
VlnPlot(gd_t, features = "type171", group.by = "condition", cols = c("red","white","white","white","white","white","white","white","white","white")) +
  geom_signif(
    comparisons = list(c("n","3"), c("n","5"), c("n","7"), c("n","9"), c("n","11"), c("n","14"), c("n","17"), c("n","21"), c("n","120")),
    map_signif_level = TRUE, textsize = 3
    , y_position = c(7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3)) +
  ylim(NA, 8)

#without T test
VlnPlot(gd_t, features = "type171", group.by = "condition", cols = c("red","white","white","white","white","white","white","white","white","white"))  +
  ylim(NA, 2)

#end session