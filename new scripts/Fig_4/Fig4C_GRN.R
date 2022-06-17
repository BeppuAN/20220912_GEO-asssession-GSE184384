#Script used to generate GRN 
#load in relevant datasets and libraries
a <- Read10X(data.dir = "/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_naive/filtered_gene_bc_matrices/mm10")
head(a)
b <- Read10X(data.dir = "/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_11dpi/mm10")
head(b)
c <- Read10X(data.dir = "/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_14dpi/mm10")
head(c)
d <- Read10X(data.dir = "/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_17dpi/mm10")
head(d)
e <- Read10X(data.dir = "/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_21dpi/mm10")
head(e)

genesN <- read.delim("/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_naive/filtered_gene_bc_matrices/mm10/genes.tsv", header=FALSE)
View(genesN)
genes11 <- read.delim("/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_11dpi/mm10/genes.tsv", header=FALSE)
View(genes11)
genes14 <- read.delim("/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_14dpi/mm10/genes.tsv", header=FALSE)
View(genes14)
genes17 <- read.delim("/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_17dpi/mm10/genes.tsv", header=FALSE)
View(genes17)
genes21 <- read.delim("/Users/beppuan/Desktop/20201030.big_scale_troubleshoot/epi_21dpi/mm10/genes.tsv", header=FALSE)
View(genes21)

#this time around, make sure that the first column in the "gene.names" input is removed so that only the Gene symbols are left
#last time this code was run the first column was the ensembl ID, which made the final Gene regulatory network difficult to interpret
genesN$V1 <- NULL
View(genesN)
genesN <- as.matrix(genesN) #essential that the these inputs are converted into matricies

genes11$V1 <- NULL #column 1 contains ensemblIDs whereas column 2 contains gene symbols. column 1 is deleted in this context to force column2 as input for the Bigscale function
View(genes11)
genes11 <- as.matrix(genes11)

genes14$V1 <- NULL
View(genes14)
genes14 <- as.matrix(genes14)

genes17$V1 <- NULL
View(genes17)
genes17 <- as.matrix(genes17)

genes21$V1 <- NULL
View(genes21)
genes21 <- as.matrix(genes21)

#network model is generated from the dataframes (which are concatenated together).
model_1 <- compute.network.model(expr.data = cbind(a,b))
model_2 <- compute.network.model(expr.data = cbind(a,c))
model_3 <- compute.network.model(expr.data = cbind(a,d))
model_4 <- compute.network.model(expr.data = cbind(a,e))

#compute network then save
results.11 <- compute.network(expr.data = b, gene.names = genes11, model = model_1)
saveRDS(results.11, file = "results.11.V2.rds")
toCytoscape(G = results.11$graph, file.name = 'results.11.V2.json')

#compute network then save
results.14 <- compute.network(expr.data = c, gene.names = genes14, model = model_2)
saveRDS(results.14, file = "results.14.V2.rds")
toCytoscape(G = results.14$graph, file.name = 'results.14.V2.json')

#compute network then save
results.17 <- compute.network(expr.data = d, gene.names = genes17, model = model_3)
saveRDS(results.17, file = "results.17.V2.rds")
toCytoscape(G = results.17$graph, file.name = 'results.17.V2.json')

#compute network then save
results.21 <- compute.network(expr.data = e, gene.names = genes21, model = model_4)
saveRDS(results.21, file = "results.21.V2.rds")
toCytoscape(G = results.21$graph, file.name = 'results.21.V2.json')