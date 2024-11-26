if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if(!require('reticulate')) install.packages('reticulate'); library('reticulate')
library(dplyr)
if (!require('Seurat')) remotes::install_version('Seurat', version='5.0.0'); library(Seurat)
library(ggplot2)
library(patchwork)
if (!require('sctransform')) install.packages('sctransform'); library(sctransform)
if (!require('BiocManager')) install.packages('BiocManager'); library('BiocManager')
if (!require('multtest')) BiocManager::install("multtest"); library('multtest')
if (!require('SeuratData')) {install.packages("remotes"); remotes::install_github("satijalab/seurat-data")}
library(SeuratData)
options(bitmapType='cairo')
if (!require('ComplexHeatmap')) devtools::install_github("jokergoo/ComplexHeatmap"); library('ComplexHeatmap')
if (!require('NMF')) install.packages('NMF'); library('NMF')
if (!require('circlize')) devtools::install_github("jokergoo/circlize"); library('circlize')
if (!require('BiocNeighbors')) devtools::install_github("LTLA/BiocNeighbors"); library('BiocNeighbors')
if (!require('BiocGenerics')) BiocManager::install("BiocGenerics"); library('BiocGenerics')
if (!require('Biobase')) BiocManager::install("Biobase"); library('Biobase')
if (!require('presto')) devtools::install_github("immunogenomics/presto"); library('presto')
if (!require('CellChat')) devtools::install_github("jinworks/CellChat"); library('CellChat')
if (!require ('nloptr')) install.packages('nloptr'); library('nloptr')
if (!require('CellChat')) devtools::install_github("jinworks/CellChat"); library('CellChat')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library(RColorBrewer)
if (!require('ggtext')) install.packages('ggtext'); library(ggtext)
source("~/Modified Bubble Function RdBu.R")
options(stringsAsFactors = FALSE)

#Create CellChat Objects for Control and VaD
object <- readRDS("~/FinalClustering.rds")
object <- subset(x = object, idents = c('Oligo_1', 'Oligo_2', 'Astrocytes', 'Microglia', 'OPCs', 'Oligo_3'))
new.cluster.ids <- c('Oligodendrocytes', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPCs', 'Oligodendrocytes')
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)
DimPlot(object, reduction = "umap", label = TRUE)

obj.list <- SplitObject(object, split.by = "group")
CTRL <- obj.list[["Control"]]
VAD <- obj.list[["VaD"]]
rm(obj.list)

control.data.input <- GetAssayData(CTRL, assay = "RNA", layer = "data")
CTRL@meta.data$clusters <- CTRL@active.ident
labels <- Idents(CTRL)
meta <- data.frame(group = labels, row.names = names(labels))
cellchatCT <- createCellChat(object = control.data.input, meta = meta, group.by = "group")

vad.data.input <- GetAssayData(VAD, assay = "RNA", layer = "data")
VAD@meta.data$clusters <- VAD@active.ident
labels <- Idents(VAD)
meta <- data.frame(group = labels, row.names = names(labels))
cellchatVD <- createCellChat(object = vad.data.input, meta = meta, group.by = "group")

#####CellChat R-L Analysis
cellchat <- cellchatCT; s <- "Control"
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use #Default, uses entire database
cellchat <- CellChat::subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, population.size = TRUE, raw.use = FALSE) #Receptor Ligand interactions
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchatCT <- cellchat

cellchat <- cellchatVD; s <- "VaD"
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use #Default, uses entire database
cellchat <- CellChat::subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, population.size = TRUE, raw.use = FALSE) #Receptor Ligand interactions
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchatVD <- cellchat

#Merge CellChat Objects
cellchatCT <- netAnalysis_computeCentrality(cellchatCT, slot.name = "netP")
cellchatVD <- netAnalysis_computeCentrality(cellchatVD, slot.name = "netP")
cclist <- list(CTRL = cellchatCT, VAD = cellchatVD)
cellchatM <- mergeCellChat(cclist, add.names = names(cclist))
rm(cellchatCT, cellchatVD)

#######Visualize Results
#Bar Chart: Total Counts CT v VD (Figure 1E)
gg1 <- compareInteractions(cellchatM, show.legend = F, group = c(1,2), color.use = c('deepskyblue3', 'red1'))
gg1 #Figure 1E

#Heatmap: Total Counts CT v VD by Cell Type (Figure 1F)
gg1 <- netVisual_heatmap(cellchatM)
gg1 #Figure 1F


#Signaling Pathways (manually selected):
x <- 'Astrocytes' #Figure 2C
s <- c('ADGRB', 'ADGRL', 'ANGPT', 'ANGPTL', 'ApoE', 'APP', 'CALCR', 'CNTN', 'FGF', 'FLRT', 'GALECTIN', 'LAMININ', 'MIF', 'NCAM', 'NECTIN', 'NGL', 'NRXN', 'PECAM2', 'PROS', 'PTN', 'PTPR', 'SEMA3', 'SEMA4', 'SEMA6', 'SEMA7', 'SLITRK', 'TENASCIN', 'UNC5')
x <- 'Microglia' #Figure 3C
s <- c('ADGRB', 'ADGRL', 'ANGPT', 'ANGPTL', 'ApoE', 'APP', 'CD45', 'CDH', 'CNTN', 'EGF', 'EPHA', 'GALECTIN', 'GAP', 'LAMININ', 'MIF', 'NRG', 'PECAM1', 'PECAM2', 'PROS', 'PTPR', 'SEMA3', 'SEMA4', 'SEMA7', 'SLIT', 'SLITRK', 'TENASCIN', 'THBS', 'UNC5')
x <- 'Oligodendrocytes' #Figure 4C
s <- c('ADGRL', 'ADGRB', 'ANGPT', 'ANGPTL', 'ApoE', 'APP', 'CD45', 'CDH', 'CLDN', 'CNTN', 'CSF', 'EGF', 'EPHA', 'FGF', 'GALECTIN', 'GAP', 'MAG', 'MIF', 'NECTIN', 'Netrin', 'NGL', 'NRG', 'NRXN', 'PROS', 'PTN', 'PTPR', 'PTPRM', 'SEMA3', 'SEMA4', 'SEMA6', 'SEMA7', 'SLIT', 'SLITRK', 'SPP1', 'TENASCIN', 'THBS', 'UNC5')
x <- 'OPCs' #Figure 5C
s <- c('ADGRL', 'ANGPT', 'ANGPTL', 'ApoE', 'APP', 'CADM', 'CALCR', 'CD45', 'CNTN', 'COMPLEMENT', 'EPHA', 'FGF', 'FLRT', 'GALECTIN', 'LAMININ', 'MIF', 'NCAM', 'NECTIN', 'NEGR', 'Netrin', 'NGL', 'NRG', 'NRXN', 'PROS', 'PTPR', 'SEMA3', 'SEMA4', 'SEMA6', 'SEMA7', 'SLIT', 'SLITRK', 'SPP1', 'TENASCIN', 'THBS', 'UNC5')

t1 <- paste("Signals From ", x, sep ="")
t2 <- paste("Signals Received By ", x, sep = "")
gg1 <- rankNet(cellchatM, mode = "comparison", stacked = T, do.stat = FALSE, color.use = c('#7281FC', '#FC553E'), sources.use = x, title = t1, signaling = s) + theme(axis.text.y = element_text(color = "black"))
gg2 <- rankNet(cellchatM, mode = "comparison", stacked = T, do.stat = FALSE, color.use = c('#7281FC', '#FC553E'), targets.use = x, title = t2, signaling = s) + theme(axis.text.y = element_text(color = "black"))
gg1 + gg2


#R-L Interactions
#Figure 2D
a <- "Astrocytes"
signaling <- c('ANGPT', 'ApoE', 'MIF', 'PECAM2', 'PROS')
b <- c('Oligodendrocytes', 'Microglia', 'Astrocytes', 'OPCs')
CellChatDB <- CellChatDB.human
x <- signaling
df <- subset(AllPaths, subset = pathway_name %in% x)
col <- brewer.pal(n = 10, name = "RdBu")
col <- rev(col)
#Cluster signaling targeting all other cells
gg1 <- modBubble(cellchatM, sources.use = a, targets.use = b,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg1 <- gg1 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90))
gg1
#All signaling received by cluster
gg2 <- modBubble(cellchatM, sources.use = b, targets.use = a,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg2 <- gg2 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90)) 
gg2


#Figure 3D
a <- "Microglia"
#signaling <- c('ApoE', 'GALECTIN', 'MIF', 'PROS')
LRList <- c('APOE_TREM2_TYROBP', 'LGALS9_CD45', 'LGALS9_HAVCR2', 'LGALS9_CD44', 'MIF_CD74_CXCR4', 'MIF_CD74_CD44', 'MIF_ACKR3', 'PROS1_AXL', 'PROS1_TYRO3', 'PROS1_MERTK')
b <- c('Oligodendrocytes', 'Microglia', 'Astrocytes', 'OPCs')
CellChatDB <- CellChatDB.human
AllPaths <- as.data.frame(CellChatDB$interaction)
x <- LRList
df <- subset(AllPaths, subset = interaction_name %in% x)
col <- brewer.pal(n = 10, name = "RdBu")
col <- rev(col)
#Cluster signaling targeting all other cells
gg1 <- modBubble(cellchatM, sources.use = a, targets.use = b,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg1 <- gg1 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90))
gg1
#All signaling received by cluster
gg2 <- modBubble(cellchatM, sources.use = b, targets.use = a,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg2 <- gg2 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90)) 
gg2


#Figure 4D
a <- "Oligodendrocytes"
LRList <- c('CLDN11_CLDN11', 'MAG_MAG', 'NRG2_ERBB3', 'NRG2_ERBB2_ERBB3', 'NRG2_ERBB4', 'NRG2_ERBB2_ERBB4', 'NRG3_ERBB4', 'NRG3_ERBB2_ERBB4')
b <- c('Oligodendrocytes', 'Microglia', 'Astrocytes', 'OPCs')
CellChatDB <- CellChatDB.human
AllPaths <- as.data.frame(CellChatDB$interaction)
x <- LRList
df <- subset(AllPaths, subset = interaction_name %in% x)
col <- brewer.pal(n = 10, name = "RdBu")
col <- rev(col)
#Cluster signaling targeting all other cells
gg1 <- modBubble(cellchatM, sources.use = a, targets.use = b,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg1 <- gg1 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90))
gg1
#All signaling received by cluster
gg2 <- modBubble(cellchatM, sources.use = b, targets.use = a,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg2 <- gg2 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90)) 
gg2


#Figure 5D
a <- "OPCs"
LRList <- c('FGF1_FGFR1', 'FGF1_FGFR2', 'FGF1_FGFR3', 'FGF2_FGFR1', 'FGF2_FGFR2', 'FGF2_FGFR3')
b <- c('Oligodendrocytes', 'Microglia', 'Astrocytes', 'OPCs')
CellChatDB <- CellChatDB.human
AllPaths <- as.data.frame(CellChatDB$interaction)
x <- LRList
df <- subset(AllPaths, subset = interaction_name %in% x)
col <- brewer.pal(n = 10, name = "RdBu")
col <- rev(col)
#Cluster signaling targeting all other cells
gg1 <- modBubble(cellchatM, sources.use = a, targets.use = b,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg1 <- gg1 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90))
gg1
#All signaling received by cluster
gg2 <- modBubble(cellchatM, sources.use = b, targets.use = a,  comparison = c(1, 2), color.heatmap = "Spectral", pairLR.use = df['interaction_name'])
gg2 <- gg2 + theme(axis.text.x = element_text(color = c("blue3","red1"), size = 7, angle = 90)) 
gg2