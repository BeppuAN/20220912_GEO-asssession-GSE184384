# 20220912_GEO-asssession-GSE184384 
scRNA script used to generate data for GEO asssession: GSE184384

# System Requirements

## Hardware Requirements
Hardware Requirements

Processor: Intel® Core™ i5-5575R Processor
RAM: 16+ GB

Data were generated using the aforementioned computer specs  

No non-standard hardware was used
  
## Software Requirements

### OS Requirements
Mac OSX: 10.11.6
Windows: 10

All packages used were compatiable with Mac and Windows Operating systems.

# Installation Guide

### Package Installation

From an `R` session, type and run:

```
install.packages('Seurat')
library(Seurat)
```

Installation of the package should take approcimately 20 seconds on a computer with recomennded specifications

# Demo

Scripts for generation of scRNAseq data are organized in this repository based upon their corresponding figure in the manuscript. Processed data is found on GSE184384. For this demo, we will be running the script labeled 'Fig1F_Dotplot.R'.

```
#code used to assess gene epression of BC genes in candidate non-BC progenitors in response to flu infection
#load in relevant datasets and libraries
library(Seurat)
library(ggplot2)
library(cowplot)
epitimecourse.integrated.V2 <- readRDS("/Volumes/20190816.singlecellRNAseq/20220112.scRNAseq_Scripts/epitimecourse.integrated.V2.rds")
x <- epitimecourse.integrated.V2

#Subset cells of interest
Idents(x) <- "anno"
table(Idents(x))
Serous <- subset(x, idents = c("Serous"))
Club <- subset(x, idents = c("Club"))
TII <- subset(x, idents = c("TII"))

#Generate Dotplots to assess gene expression as a function of "days post infection"
DotPlot(Serous, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Bpifa1", "rna_Ltf", "rna_Scgb3a2"), group.by = "timepoint", scale.max = 10)
DotPlot(Club, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Hp", "rna_Scgb1a1", "rna_Cldn10"), group.by = "timepoint", scale.max = 10)
DotPlot(TII, features = c("rna_Krt5", "rna_Krt14", "rna_Trp63", "rna_Cd74", "rna_Chil1", "rna_Sftpc"), group.by = "timepoint", scale.max = 10)
```

Data should show cell-type specific changes in gene experssion of common epitheliall cell marker as a function of time following PR8 influenza exposure. 

Generation of Dotplots should take approcimately 20 seconds on a computer with recomennded specifications.

# Instructions for use

No special considerations on how to run the software on the data
