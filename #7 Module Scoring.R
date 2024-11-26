library(dplyr)
if (!require('Seurat')) install.packages('Seurat'); library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(RColorBrewer)
if (!require('sctransform')) install.packages('sctransform'); library(sctransform)
if (!require('BiocManager')) install.packages('BiocManager'); library('BiocManager')
if (!require('multtest')) BiocManager::install("multtest"); library('multtest')
if (!require('SeuratData')) {install.packages("remotes"); remotes::install_github("satijalab/seurat-data")}
library(SeuratData)
if (!require('HGNChelper')) install.packages('HGNChelper'); library('HGNChelper')
require(gridExtra)
options(bitmapType='cairo')

sansing <- readRDS(file = "~/FinalClustering.rds")
new.cluster.ids <- c('Oligodendrocytes', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPCs', 'Neurons', 'Oligodendrocytes', 'Endothelial Cells')
names(new.cluster.ids) <- levels(sansing)
sansing <- RenameIdents(sansing, new.cluster.ids)
cardapi <- readRDS("~/MitroiDAPIFinalClustering.rds")
cardapi <- subset(x = cardapi, idents = c("Oligodendrocytes", "Astrocytes", "OPCs", "Neurons"))
carerg <- readRDS("~/MitroiERGFinalClustering.rds")
new.cluster.ids <- c('Microglia', 'Oligodendrocytes', 'Microglia', 'Microglia', 'Microglia', 'Astrocytes', 'Microglia', 'T-Cells', 'Endothelial Cells', 'OPCs', 'Macrophages')
names(new.cluster.ids) <- levels(carerg)
carerg <- RenameIdents(carerg, new.cluster.ids)

#Load the datasets from Reactome
if (!require('msigdbr')) install.packages("msigdbr"); library("msigdbr")
all_gene_sets = msigdbr(species = "Homo sapiens")

###############################################################
names <- c("REACTOME_CELLULAR_SENESCENCE", 
           "REACTOME_CELL_DEATH_SIGNALLING_VIA_NRAGE_NRIF_AND_NADE", 
           "REACTOME_HSF1_ACTIVATION",
           "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR")

#Implement Module Scoring
source <- c("Sansing", "CARDAPI", "CARERG")
seurat <- c(sansing, cardapi, carerg)

for (x in names){
  list1 <- subset(all_gene_sets, gs_name == x)
  list1 <- as.vector(list1$gene_symbol)
  a <- AddModuleScore(seurat[[1]], features = list(list1), name=x)
  b <- AddModuleScore(seurat[[2]], features = list(list1), name=x)
  c <- AddModuleScore(seurat[[3]], features = list(list1), name=x)
  feats <- paste0(x, "1")
  minimal <- c(min(a[[feats]]), min(b[[feats]]), min(c[[feats]]))
  minimal <- min(minimal)
  maximal <- c(max(a[[feats]]), max(b[[feats]]), max(c[[feats]]))
  maximal <- max(maximal)
  obj.list1 <- SplitObject(a, split.by = "group")
  obj.list2 <- SplitObject(b, split.by = "group")
  obj.list3 <- SplitObject(c, split.by = "group")
  gg1 <- FeaturePlot(obj.list1$Control, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "Control")
  gg2 <- FeaturePlot(obj.list1$VaD, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "VaD")
  gg3 <- FeaturePlot(obj.list2$Control, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "Control")
  gg4 <- FeaturePlot(obj.list2$VaD, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "VaD")
  gg5 <- FeaturePlot(obj.list3$Control, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "Control")
  gg6 <- FeaturePlot(obj.list3$VaD, features = feats, label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PuOr")), limits = c(minimal, maximal)) + labs(title = "VaD")
  wp <- wrap_plots(gg1, gg2, gg3, gg4, gg5, gg6, guides = 'collect', ncol = 2) + plot_annotation(x, theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")))
  pdf(paste("~/ModulePlots/", x, ".pdf", sep = ""), width = 8, height = 10)
  print(wp)
  dev.off()
}