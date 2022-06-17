#code used to assess gene set enrichment for cell cycle genes as a function of time
#load in relevant datasets and libraries
library(Seurat)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#subset serous cells
table(Idents(x))
Serous <- subset(x, idents = c("Serous"))

#set ident to timepoint
Idents(Serous) <- "timepoint"

#generate deseq list (one of two inputs for the geneset enrichment analysis)
c <- FindMarkers(Serous, ident.1 = "5 dpi", ident.2 = "trachea")

#pipeline DE list into FGSE package
#gene symbols need to be copied into there own columns - use "setDT" function in "data.table" 
library(data.table)
c <- setDT(c, keep.rownames = T)
head(c)

#One of the first things you'll have to do is create a mapping table - a list of ensemble Id's and their associated symbols
#we will generate this object with the "biomart" package
library(biomaRt)
library(dplyr)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()
head(bm)
bm

#we have "mt-" in front of every gene 
#remove character strings by position using the gsub function
bm$external_gene_name <- gsub('mt-', '', bm$external_gene_name)
bm$hsapiens_homolog_associated_gene_name <- gsub('MT-', '', bm$hsapiens_homolog_associated_gene_name)
bm

#combine your DE lists with the mapping table made in the previous step
head(c)
c <- inner_join(c, bm, by=c("rn"="external_gene_name"))
head(c)

#generate ranks (1 of 2 inputs for FGSEA)
c <- c %>% 
  dplyr::select(hsapiens_homolog_associated_gene_name, avg_logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  summarize(avg_logFC=mean(avg_logFC))
head(c)

#this step is essential - genes need to be ranked from highest to lowest expression  
#load package 'tibble' to gain access to deframe command line
library(tibble)
c <- deframe(c)
c <- sort(c, decreasing=T)

#the second input for GSEA are the gene sets that are going to be tested. 
#for this analysis, we will use MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/
#download the gene set and place the .txt file into your workspace
#you won't be able to load in the pathway with out command line "gmtpathways" which is part of the "FGSEA" package
library(fgsea)
pathways.hallmark <- gmtPathways("h.all.v7.1.symbols.gmt")

#let's take a look at what genes correspond to which pathway
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

#Now, run the fgsea algorithm with 1000 permutations:
fgsearanks <- fgsea(pathways=pathways.hallmark, stats=c, nperm=1000)

#Tidy the results:
fgsearanksTidy <- fgsearanks %>%
  as_tibble() %>%
  arrange(desc(NES))

#generate table
fgsearanksTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
