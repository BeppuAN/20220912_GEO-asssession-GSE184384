#Script used to generate GRN. 
#load in relevant datasets and libraries.
flugenes.timecourse.postintegration.V21 <- readRDS("~/Desktop/20201030.big_scale_troubleshoot/20201111.Peter_RO1/flugenes.timecourse.postintegration.V21.rds")

library(data.table)
library(Matrix)
library("bigSCale")
library("FNN")
library("Seurat")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")
library("scater")

install.packages("umap")
library("umap")

#first subset by celltype.
DimPlot(flugenes.timecourse.postintegration.V21)
Basal <- subset(flugenes.timecourse.postintegration.V21, ident = "Basal (Krt5)")
Serous <- subset(flugenes.timecourse.postintegration.V21, ident = "Serous-like")

#subset by timepoint: We need naive and 14 dpi.
#subsetting by celltype first has the benefit of reducing the number of input cells. therefore, we wont be downsampling for this run.
Idents(Basal) <- "timepoint"
Basalnaive <- subset(Basal, ident = c("trachea"))
Basalnaive <- Basalnaive@assays$RNA@data
Basal14 <- subset(Basal, ident = c("11 dpi", "14 dpi", "17 dpi"))
Basal14 <- Basal14@assays$RNA@data

Idents(Serous) <- "timepoint"
Serousnaive <- subset(Serous, ident = c("trachea"))
Serousnaive <- Serousnaive@assays$RNA@data
Serous14 <- subset(Serous, ident = c("11 dpi", "14 dpi", "17 dpi"))
Serous14 <- Serous14@assays$RNA@data

#compute models
model_1 <- compute.network.model(expr.data = cbind(Basalnaive, Basal14))
saveRDS(model_1, file = "model_1.rds")
model_2 <- compute.network.model(expr.data = cbind(Serousnaive, Serous14))
saveRDS(model_2, file = "model_2.rds")

#import gene names. Order of the gene name should match how its list in the singlecell data matrix.
genes <- read.delim("/home/strippuser/Desktop/20201030.big_scale_troubleshoot/epi_naive/filtered_gene_bc_matrices/mm10/genes.tsv", header=FALSE)
View(genes)
genes$V1 <- NULL
View(genes)
genes <- as.matrix(genes) #the "Compute.network" function only accepts gene lists in matrix format. This line is essential to run the code

#Compute Node centralities
#Basal cells: naive and 14 dpi
results.Basalnaive <- compute.network(expr.data = Basalnaive, gene.names = genes, model = model_1)
saveRDS(results.Basalnaive, file = "results.Basalnaive.rds")
results.Basal14 <- compute.network(expr.data = Basal14, gene.names = genes, model = model_1)
saveRDS(results.Basal14, file = "results.Basal14.rds")

#Serous cells: naive and 14 dpi
results.Serousnaive <- compute.network(expr.data = Serousnaive, gene.names = genes, model = model_2)
saveRDS(results.Serousnaive, file = "results.Serousnaive.rds")
results.Serous14 <- compute.network(expr.data = Serous14, gene.names = genes, model = model_2)
saveRDS(results.Serous14,  file = "results.Serous14.rds")

#homogenize networks
output <- homogenize.networks(list(results.Basalnaive,results.Basal14))
results.Basalnaive.homo <- output[[1]]
results.Basal14.homo <- output[[2]]
output <- homogenize.networks(list(results.Serousnaive,results.Serous14))
results.Serousnaive.homo <- output[[1]]
results.Serous14.homo <- output[[2]]

#save homogenized networks
saveRDS(results.Basalnaive.homo, file = "results.Basalnaive.homo.rds")
saveRDS(results.Basal14.homo, file = "results.Basal14.homo.rds")
saveRDS(results.Serousnaive.homo, file = "results.Serousnaive.homo.rds")
saveRDS(results.Basal14.homo, file = "rresults.Basal14.homo.rds")

#compute Delta Degree centrality
comparison <- compare.centrality(list(results.Basalnaive.homo$centrality,results.Basal14.homo$centrality),c('results.Basalnaive.homo','results.Basal14.homo'))
DT::datatable(comparison$Degree)
comparison <- compare.centrality(list(results.ctl$centrality,results.t2d$centrality),c('results.Basalnaive.homo','results.Basal14.homo'))
DT::datatable(comparison$Degree)
