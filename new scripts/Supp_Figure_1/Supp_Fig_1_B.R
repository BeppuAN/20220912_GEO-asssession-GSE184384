#script to generate Supp figure 1B. - Temporal gene expression of antimicrobial factors
#load in relavant libraries and datasets

influenza_timecourse_epithelial <- readRDS("/Volumes/20190816.singlecellRNAseq/20210909_Geo_upload_V2/rds/influenza_timecourse_epithelial.rds")

#subset serous
influenza_timecourse_epithelial@meta.data
serous <- subset(influenza_timecourse_epithelial, idents = "Serous")

#Vlnplot order 
# Define an order of cluster identities
table(serous@meta.data$timepoint)
my_levels <- c("trachea", "naive", "5 dpi", "7 dpi", "9 dpi", "11 dpi", "14 dpi", "17 dpi", "21 dpi", "4 mpi", "8 mpi")
my_levels <- c("8 mpi", "4 mpi", "21 dpi", "17 dpi", "14 dpi", "11 dpi", "9 dpi", "7 dpi", "5 dpi", "naive", "trachea")

# Relevel object@ident
serous@meta.data$timepoint <- factor(serous@meta.data$timepoint, levels = my_levels)

VlnPlot(serous, features = c("rna_Ifitm1"), group.by = "timepoint")
VlnPlot(serous, features = c("rna_Ifitm3"), group.by = "timepoint")
VlnPlot(serous, features = c("rna_Ltf"), group.by = "timepoint")
VlnPlot(serous, features = c("rna_Bpifa1"), group.by = "timepoint")

#end session