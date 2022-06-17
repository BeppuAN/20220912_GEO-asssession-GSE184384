#script to generate Supp figure 2A. - Feature gene expression between serous and club cells
#load in relavant libraries and datasets
library("Seurat")
influenza_timecourse_epithelial <- readRDS("/Volumes/20190816.singlecellRNAseq/20210909_Geo_upload_V2/rds/influenza_timecourse_epithelial.rds")

#subset cells of interest
table(Idents(influenza_timecourse_epithelial))
serous_club <- subset(influenza_timecourse_epithelial, idents = c("Club", "Serous"))

#recluster
serous_club <- RunUMAP(serous_club, dims = 1:30)
DimPlot(serous_club)

#Generate featureplot
FeaturePlot(serous_club, features = "rna_Scgb1a1", min.cutoff = 5)
FeaturePlot(serous_club, features = "rna_Scgb3a2", min.cutoff = 2)
FeaturePlot(serous_club, features = "rna_Bpifa1")
FeaturePlot(serous_club, features = "rna_Msln")
FeaturePlot(serous_club, features = "rna_Ltf")
FeaturePlot(serous_club, features = "rna_Krt5")

#end session
