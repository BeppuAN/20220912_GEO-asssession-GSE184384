#Script to generate Supplementary Figure 1 A,B and c - Assessment of G1_S_G2M homogeniety in serous and club cells
#Load in relavant libraries and datasets
#A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
#segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
flugenes.timecourse.postintegration.V21 <- readRDS("/Volumes/20190816.singlecellRNAseq/Flupaper.Rscripts/flugenes.timecourse.postintegration.V21.rds")
a <- flugenes.timecourse.postintegration.V21

table(Idents(a))
a <- subset(a, idents = c("Club", "Serous-like"))

a <- RunPCA(a, features = VariableFeatures(a), ndims.print = 6:10, nfeatures.print = 10)
DimHeatmap(a, dims = c(8, 9))

#assign cell cycle scores
a <- CellCycleScoring(a, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(a[[]])
a@meta.data
# Visualize the distribution of cell cycle markers across
RidgePlot(a, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
DimPlot(a, reduction = "pca")
DimPlot(a, reduction = "pca", group.by = "anno6", cols = c("black", "red"))

#Regress out cell cycle scores during data scaling
a <- ScaleData(a, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(a))
a <- RunPCA(a, features = VariableFeatures(a), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
a <- RunPCA(a, features = c(s.genes, g2m.genes))

DimPlot(a, reduction = "pca", group.by = "anno6", cols = c("black", "red"))

#PCA should look homogenized after regression 
#try alternative workflow
a$CC.Difference <- a$S.Score - a$G2M.Score
a <- ScaleData(a, vars.to.regress = "CC.Difference", features = rownames(a))

# cell cycle effects strongly mitigated in PCA
a <- RunPCA(a, features = VariableFeatures(a), nfeatures.print = 10)
DimPlot(a, reduction = "pca")

#either way, regression of cell-cycle genes doesn;t seem to affect clustering of serous-like and club cells at all
#these are likely different cells types

#next thing we want to do is to assess if the porportion of G1/G2M/S genes changes between conditions
#subset cells based on subset
Idents(a) <- "anno6"
table(Idents(a))
DimPlot(a, reduction = "pca")

b <- subset(a, idents = "Serous-like")
c <- subset(a, idents = "Club")

#make sure the cells have been properly subsetted
DimPlot(b, reduction = "pca")
DimPlot(c, reduction = "pca")

#Set identity back to cell-cycle score
Idents(b) <- "Phase"
Idents(c) <- "Phase"

#use prism to generate stacked barplot
table(Idents(b))
table(Idents(c))

#end session

