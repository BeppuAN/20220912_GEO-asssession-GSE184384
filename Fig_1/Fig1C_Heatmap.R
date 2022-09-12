#code used to generate Heatmap
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#generate differential expression list
#use roc test to find clean markers
amarkers <- FindAllMarkers(x, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

top10 <- amarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_diff)
x <- ScaleData(x, verbose = F, assay = "RNA")

DoHeatmap(subset(x, downsample=100), features = top10$gene) + scale_fill_gradientn(colors = c("light blue", "white", "red"))


