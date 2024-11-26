#Load libraries
library(dplyr)
if (!require('Seurat')) install.packages('Seurat'); library(Seurat)
library(ggplot2)
library(patchwork)
if (!require('sctransform')) install.packages('sctransform'); library(sctransform)
if (!require('BiocManager')) install.packages('BiocManager'); library('BiocManager')
if (!require('multtest')) BiocManager::install("multtest"); library('multtest')
if (!require('SeuratData')) {install.packages("remotes"); remotes::install_github("satijalab/seurat-data")}
library(SeuratData)
options(bitmapType='cairo')

output <- readRDS("~/FinalClustering.rds")
levels(output)
DefaultAssay(output) <- "RNA"

#Differential Gene Expression analysis
#Positive values are higher in VaD
clusters <- levels(output)
print(clusters)
for (x in clusters) {
  DGE <- FindMarkers(output, ident.1 = "VaD", ident.2 = "Control", verbose = TRUE, group.by="group", subset.ident = x, logfc.threshold = 0, min.pct = 0)
  file = paste("~/DEGs/", x, ".csv", sep = "")
  write.csv(DGE, file = file)
}