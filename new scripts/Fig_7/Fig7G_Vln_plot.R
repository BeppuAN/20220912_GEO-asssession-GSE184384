#RScrtipt for figure 7 Vlnplot
#We will perform T test using the package ggsignif
#link:https://github.com/const-ae/ggsignif
#download relavant libraries and datasets
IL22ra1KO_21dpi <- readRDS("/Volumes/20190816.singlecellRNAseq/20210909_Geo_upload_V2/rds/IL22ra1KO_21dpi.rds")

#install ggsignif
install.packages("ggsignif")
library(ggplot2)
library(ggsignif)

#make sure correct dataset was loaded in
DimPlot(IL22ra1KO_21dpi)

#subset serous and basal clusters for downstream analysis
serous_basal <- subset(IL22ra1KO_21dpi, idents = c("Serous", "Basal"))

# reorder conditions to have WT condition displayed first
# Define an order of cluster identities
my_levels <- c("WT","IL22ra1KO")
serous_basal@meta.data$condition <- factor(serous_basal@meta.data$condition, levels = my_levels)

#generate Vlnplot with "ggsignif" ammendments
#t-test used is the wilcox t test
VlnPlot(serous_basal, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Itga6", "rna_Ngfr"), group.by = "condition", cols = c("red", "black"))
VlnPlot(serous_basal, features = c("rna_Scgb3a2", "rna_Bpifa1", "rna_Ltf", "rna_Msln"), group.by = "condition", cols = c("red", "black"))

?geom_signif

P1 <- VlnPlot(serous_basal, features = c("rna_Krt5"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 5)

P2 <- VlnPlot(serous_basal, features = c("rna_Itga6"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 3)

P3 <- VlnPlot(serous_basal, features = c("rna_Trp63"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 3.5)

P4 <- VlnPlot(serous_basal, features = c("rna_Krt14"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 6)

P5 <- VlnPlot(serous_basal, features = c("rna_Ngfr"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 1.5)

P6 <- VlnPlot(serous_basal, features = c("rna_Scgb3a2"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 8)

P7 <- VlnPlot(serous_basal, features = c("rna_Bpifa1"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 10)

P8 <- VlnPlot(serous_basal, features = c("rna_Ltf"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 8)

P9 <- VlnPlot(serous_basal, features = c("rna_Msln"), group.by = "condition", cols = c("red", "white")) +
  geom_signif(
    comparisons = list(c("WT","IL22ra1KO")),
    map_signif_level = TRUE, textsize = 3
  ) +
  ylim(NA, 7)

(P1 + P3 + P4 + P7 + P8 + P9)

#end session
