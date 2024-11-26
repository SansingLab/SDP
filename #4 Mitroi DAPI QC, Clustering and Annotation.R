#Load libraries
library(dplyr)
if (!require('Seurat')) remotes::install_version('Seurat', version='5.0.0'); library(Seurat)
library(ggplot2)
library(patchwork)
if (!require('sctransform')) install.packages('sctransform'); library(sctransform)
if (!require('BiocManager')) install.packages('BiocManager'); library('BiocManager')
if (!require('multtest')) BiocManager::install("multtest"); library('multtest')
if (!require('SeuratData')) {install.packages("remotes"); remotes::install_github("satijalab/seurat-data")}
library(SeuratData)
if (!require('HGNChelper')) install.packages('HGNChelper'); library('HGNChelper')
options(bitmapType='cairo')
getwd()

#Load raw data counts (TSV files):
dapi.ctrl.data <- Read10X(data.dir = "~/DAPI/CTRL/2526/")
dapi.ctrl.data2 <- Read10X(data.dir = "~/DAPI/CTRL/4320/")
dapi.ctrl.data3 <- Read10X(data.dir = "~/DAPI/CTRL/4582/")
dapi.ctrl.data4 <- Read10X(data.dir = "~/DAPI/CTRL/4621/")
dapi.ctrl.data5 <- Read10X(data.dir = "~/DAPI/CTRL/4858/")
dapi.vad.data <- Read10X(data.dir = "~/DAPI/VAD/3505F/")
dapi.vad.data2 <- Read10X(data.dir = "~/DAPI/VAD/3761F/")
dapi.vad.data3 <- Read10X(data.dir = "~/DAPI/VAD/4143F/")
dapi.vad.data4 <- Read10X(data.dir = "~/DAPI/VAD/HBEMA001F/")
dapi.vad.data5 <- Read10X(data.dir = "~/DAPI/VAD/UCLAF/")

# Initialize the Seurat object with the raw (non-normalized) data.
DCTRL <- CreateSeuratObject(counts = dapi.ctrl.data, project = "DAPIControl1", min.cells = 3, min.features = 200)
DCTRL$group <- "Control"
DCTRL2 <- CreateSeuratObject(counts = dapi.ctrl.data2, project = "DAPIControl2", min.cells = 3, min.features = 200)
DCTRL2$group <- "Control"
DCTRL3 <- CreateSeuratObject(counts = dapi.ctrl.data3, project = "DAPIControl3", min.cells = 3, min.features = 200)
DCTRL3$group <- "Control"
DCTRL4 <- CreateSeuratObject(counts = dapi.ctrl.data4, project = "DAPIControl4", min.cells = 3, min.features = 200)
DCTRL4$group <- "Control"
DCTRL5 <- CreateSeuratObject(counts = dapi.ctrl.data5, project = "DAPIControl5", min.cells = 3, min.features = 200)
DCTRL5$group <- "Control"

DVAD <- CreateSeuratObject(counts = dapi.vad.data, project = "DAPIVaD1", min.cells = 3, min.features = 200)
DVAD$group <- "VaD"
DVAD2 <- CreateSeuratObject(counts = dapi.vad.data2, project = "DAPIVaD2", min.cells = 3, min.features = 200)
DVAD2$group <- "VaD"
DVAD3 <- CreateSeuratObject(counts = dapi.vad.data3, project = "DAPIVaD3", min.cells = 3, min.features = 200)
DVAD3$group <- "VaD"
DVAD4 <- CreateSeuratObject(counts = dapi.vad.data4, project = "DAPIVaD4", min.cells = 3, min.features = 200)
DVAD4$group <- "VaD"
DVAD5 <- CreateSeuratObject(counts = dapi.vad.data5, project = "DAPIVaD5", min.cells = 3, min.features = 200)
DVAD5$group <- "VaD"

#Filtering based on mitochondrial percent and number of counts
DCTRL[["percent.mt"]] <- PercentageFeatureSet(DCTRL, pattern = "^MT-")
DCTRL2[["percent.mt"]] <- PercentageFeatureSet(DCTRL2, pattern = "^MT-")
DCTRL3[["percent.mt"]] <- PercentageFeatureSet(DCTRL3, pattern = "^MT-")
DCTRL4[["percent.mt"]] <- PercentageFeatureSet(DCTRL4, pattern = "^MT-")
DCTRL5[["percent.mt"]] <- PercentageFeatureSet(DCTRL5, pattern = "^MT-")

DVAD[["percent.mt"]] <- PercentageFeatureSet(DVAD, pattern = "^MT-")
DVAD2[["percent.mt"]] <- PercentageFeatureSet(DVAD2, pattern = "^MT-")
DVAD3[["percent.mt"]] <- PercentageFeatureSet(DVAD3, pattern = "^MT-")
DVAD4[["percent.mt"]] <- PercentageFeatureSet(DVAD4, pattern = "^MT-")
DVAD5[["percent.mt"]] <- PercentageFeatureSet(DVAD5, pattern = "^MT-")

# Filter by QC Metrics
DCTRL <- subset(DCTRL, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DCTRL2 <- subset(DCTRL2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DCTRL3 <- subset(DCTRL3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DCTRL4 <- subset(DCTRL4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DCTRL5 <- subset(DCTRL5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

DVAD <- subset(DVAD, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DVAD2 <- subset(DVAD2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DVAD3 <- subset(DVAD3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DVAD4 <- subset(DVAD4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
DVAD5 <- subset(DVAD5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

obj.list <- c(DCTRL2, DCTRL3, DCTRL4, DCTRL5, DVAD, DVAD2, DVAD3, DVAD4, DVAD5)
obj <- merge(DCTRL, y = obj.list)

#Normalization (Using scTransform) and Integration
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- SCTransform(obj, vars.to.regress = "percent.mt")
obj <- RunPCA(obj)
obj <- IntegrateLayers(object = obj, method = CCAIntegration, normalization.method = "SCT")
ktest <- FindNeighbors(obj, reduction = "integrated.dr", dims = 1:40)
ktest <- FindClusters(ktest, resolution = 0.1)
ktest <- RunUMAP(ktest, dims = 1:40, reduction = "integrated.dr")
DimPlot(ktest, reduction = "umap")
DefaultAssay(ktest) <- "RNA"
ktest <- JoinLayers(ktest)

#Cluster Annotation
#Manual: Visualize characterized pre-chosen cell-specific markers
DefaultAssay(ktest) <- "RNA"

#Microglia
VlnPlot(ktest, features = c("CX3CR1","SALL1", "P2RY12", "IL18", "CD74", "FYB1"))
VlnPlot(ktest, features = c("C3", "SIGLEC1", "CD44", "MYB"))
VlnPlot(ktest, features = c("CD163", "CD14", "PTPRC", "CD86"))
#T-Cells
VlnPlot(ktest, features = c("CD3E", "CD3D", "CD3G", "LCK", "IL7R"))
VlnPlot(ktest, features = c('CD8A', 'TRAC', 'UBASH3A', 'GPR171', 'CD69'))
#Astrocytes
VlnPlot(ktest, features = c("ALDH1L1", "GFAP", "AQP4", "FGFR3"))
VlnPlot(ktest, features = c("SOX2", "FAM107A", "EGFR"))
#Oligodendrocytes
VlnPlot(ktest, features = c("PLP1", "MBP", "CNP", "MOG"))
#Oligo-Precursor Cells
VlnPlot(ktest, features = c("PDGFRA", "TNR", 'GPR17', 'PCDH15', 'COL9A1'))
#Endothelial Cells
VlnPlot(ktest, features = c("FLT1", "VWF", "PECAM1"))
#Neurons 
VlnPlot(ktest, features = c("TMEM130", "CNTN4", "RBFOX3", "SYP", "SNAP25"))
#Inh Neurons
VlnPlot(ktest, features = c('MEPE', 'OTOF', 'VAX1', 'AKAIN1', 'GPR149'))
#ExNeurons
VlnPlot(ktest, features = c('MCHR2', 'CPNE4', 'ADAMTSL1', 'CBLN2'))
#More Neuronal Markers
VlnPlot(ktest, features = c("SLC17A7", "SATB2", "GAD2"))
#Fibroblasts
VlnPlot(ktest, features = c("GFAP", "AQP4", "COL1A1", "CCL11", "KERA", "MYOC", "VEGFD"))

#Cluster Annotation
#Rename each cluster manually based on marker expression
new.cluster.ids <- c('Oligodendrocytes', 'Oligodendrocytes', 'Oligodendrocytes', 'Oligodendrocytes', 'Astrocytes', 'Oligodendrocytes', 'OPCs', 'Neurons', 'Oligodendrocytes', 'Oligodendrocytes', 'Unidentified')
names(new.cluster.ids) <- levels(ktest)
ktest <- RenameIdents(ktest, new.cluster.ids)
DimPlot(ktest, reduction = "umap") #Figure S1

#Save the clustered Seurat object
saveRDS(ktest, file = "~/MitroiDAPIFinalClustering.rds")