#code used to assess changes in cell representation in serous vs BC
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#Subset BC and serous clusters
Idents(x) <- "anno"
table(Idents(x))
Serous <- subset(x, idents = c("Serous"))
BC <- subset(x, idents = c("BC"))

#change identity class to timepoint
#use values from table to generate graph in Prism
Idents(Serous) <- "timepoint"
Idents(BC) <- "timepoint"
table(Idents(Serous))
table(Idents(BC))

#only including cells from TII depleted datasets as per reviewer request
#ie batches 1-6 in epitimecourse.integrated.V2@meta.data$'tech'
Idents(x) <- "tech"
table(Idents(x))
x <- subset(x, idents = c("batch 1","batch 2","batch 3","batch 4","batch 5","batch 6"))
