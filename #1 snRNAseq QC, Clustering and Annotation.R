#Pre-Processing based on Seurat tutorials
library(dplyr)
if (!require('Seurat')) remotes::install_version('Seurat', version='5.0.0'); library(Seurat)
library(ggplot2)
library(patchwork)
if (!require('sctransform')) install.packages('sctransform'); library(sctransform)
if (!require('BiocManager')) install.packages('BiocManager'); library('BiocManager')
if (!require('multtest')) BiocManager::install("multtest"); library('multtest')
if (!require('SeuratData')) {install.packages("remotes"); remotes::install_github("satijalab/seurat-data")}
if (!require('HGNChelper')) install.packages('HGNChelper'); library('HGNChelper')
library(SeuratData)
options(bitmapType='cairo')
getwd()

#Load processed data counts (TSV files):
ctrl.data <- Read10X(data.dir = "/results_control/outs/filtered_feature_bc_matrix/")
vad.data <- Read10X(data.dir = "/results_VAD/outs/filtered_feature_bc_matrix/")
ctrl.data2 <- Read10X(data.dir = "/results_control2/outs/filtered_feature_bc_matrix/")
vad.data2 <- Read10X(data.dir = "/results_VAD2/outs/filtered_feature_bc_matrix/")
vad.data3 <- Read10X(data.dir = "/results_VAD3Fem/outs/filtered_feature_bc_matrix/")
ctrl.data4 <- Read10X(data.dir = "/results_CTRL4Fem/outs/filtered_feature_bc_matrix/")
ctrl.data5 <- Read10X(data.dir = "/results_CTRL5/outs/filtered_feature_bc_matrix/")
vad.data5 <- Read10X(data.dir = "/results_VAD5/outs/filtered_feature_bc_matrix/")

#Create Seurat Objects
CTRL <- CreateSeuratObject(counts = ctrl.data, project = "Control1", min.cells = 3, min.features = 200)
CTRL$group <- "Control"
VAD <- CreateSeuratObject(counts = vad.data, project = "VaD1", min.cells = 3, min.features = 200)
VAD$group <- "VaD"
CTRL2 <- CreateSeuratObject(counts = ctrl.data2, project = "Control2", min.cells = 3, min.features = 200)
CTRL2$group <- "Control"
VAD2 <- CreateSeuratObject(counts = vad.data2, project = "VaD2", min.cells = 3, min.features = 200)
VAD2$group <- "VaD"
VAD3 <- CreateSeuratObject(counts = vad.data3, project = "VaD3", min.cells = 3, min.features = 200)
VAD3$group <- "VaD"
CTRL4 <- CreateSeuratObject(counts = ctrl.data4, project = "Control4", min.cells = 3, min.features = 200)
CTRL4$group <- "Control"
CTRL5 <- CreateSeuratObject(counts = ctrl.data5, project = "Control5", min.cells = 3, min.features = 200)
CTRL5$group <- "Control"
VAD5 <- CreateSeuratObject(counts = vad.data5, project = "VaD5", min.cells = 3, min.features = 200)
VAD5$group <- "VaD"

#Quantify Mitochondrial count percentage
CTRL[["percent.mt"]] <- PercentageFeatureSet(CTRL, pattern = "^MT-")
VAD[["percent.mt"]] <- PercentageFeatureSet(VAD, pattern = "^MT-")
CTRL2[["percent.mt"]] <- PercentageFeatureSet(CTRL2, pattern = "^MT-")
VAD2[["percent.mt"]] <- PercentageFeatureSet(VAD2, pattern = "^MT-")
VAD3[["percent.mt"]] <- PercentageFeatureSet(VAD3, pattern = "^MT-")
CTRL4[["percent.mt"]] <- PercentageFeatureSet(CTRL4, pattern = "^MT-")
CTRL5[["percent.mt"]] <- PercentageFeatureSet(CTRL5, pattern = "^MT-")
VAD5[["percent.mt"]] <- PercentageFeatureSet(VAD5, pattern = "^MT-")

# Remove cells with total genes < 500, total genes > 6000 OR mito percentage > 5%
CTRL <- subset(CTRL, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VAD <- subset(VAD, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
CTRL2 <- subset(CTRL2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VAD2 <- subset(VAD2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VAD3 <- subset(VAD3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
CTRL4 <- subset(CTRL4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
CTRL5 <- subset(CTRL5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VAD5 <- subset(VAD5, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)

obj.list <- c(VAD, CTRL2, VAD2, VAD3, CTRL4, CTRL5, VAD5)
obj <- merge(CTRL, y = obj.list)

#Normalization (scTransform) and Integration
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- SCTransform(obj, vars.to.regress = "percent.mt", seed.use = 1448145)
obj <- RunPCA(obj, npcs = 50, seed.use = 42)
obj <- IntegrateLayers(object = obj, orig.reduction = "pca", method = CCAIntegration, normalization.method = "SCT")
ktest <- FindNeighbors(obj, reduction = "integrated.dr", dims = 1:40)
ktest <- FindClusters(ktest, resolution = 0.1, random.seed = 0)
ktest <- RunUMAP(ktest, dims = 1:40, reduction = "integrated.dr", seed.use = 42)
DimPlot(ktest, reduction = "umap")
DefaultAssay(ktest) <- "RNA"
ktest <- JoinLayers(ktest)
rm('obj')

#Cluster Annotation
#Visualize the expression of cell type-specific markers
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
#Renamed each cluster manually based on marker expression
new.cluster.ids <- c('Oligo_1', 'Oligo_2', 'Astrocytes', 'Microglia', 'OPCs', 'Neurons', 'Oligo_3', 'Endothelial Cells')
names(new.cluster.ids) <- levels(ktest)
ktest <- RenameIdents(ktest, new.cluster.ids)

#Save the clustered Seurat object
saveRDS(ktest, file = "~/FinalClustering.rds")

ktest <- readRDS("~/FinalClustering.rds")

#Figure 1B: UMAP
DimPlot(ktest, reduction = "umap", label = TRUE, cols = c('red', 'orange', 'yellow', 'green3', 'deepskyblue3', 'violet', 'maroon3', 'palevioletred1'))

#Figure 1C: DotPlot of Markers
dots <- ktest
new.cluster.ids <- c('Oligodendrocytes', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPCs', 'Neurons', 'Oligodendrocytes', 'Endothelial Cells')
names(new.cluster.ids) <- levels(dots)
dots <- RenameIdents(dots, new.cluster.ids)
features <- c('FLT1', 'ANGPT2', 'EGFL7', 'VCAM1', 'RBFOX3', 'SYP', 'SATB2', 'SLC17A7', 'TNR', 'GPR17', 'PCDH15', 'COL9A1', 'FYB1', 'P2RY12', 'CX3CR1', 'CD86', 'GFAP', 'AQP4', 'ALDH1L1', 'FAM107A', 'MOG', 'PLP1', 'MBP', 'CNP')
plot <- DotPlot(dots, features = features) + theme(axis.text.x = element_text(angle = 45, hjust=1))
plot

#Figure 1D:
#Cell count by cluster and sample
tbl <- table(ktest@active.ident, ktest@meta.data$group)
tbl <- as.data.frame(tbl)
ggplot(tbl, aes(fill = Var1, y = Freq, x = Var2)) + geom_bar(position="fill", stat="identity") + labs(x = "Sample", y = "Population Fraction", fill = "Cluster") + scale_fill_manual(values = c('red', 'orange', 'yellow', 'green3', 'deepskyblue3', 'violet', 'maroon3', 'palevioletred1'))