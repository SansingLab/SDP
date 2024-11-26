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

#Load raw data counts (TSV files):
ctrl.data <- Read10X(data.dir = "/results_control/outs/filtered_feature_bc_matrix/")
hypo.data <- Read10X(data.dir = "/results_hypoxia/outs/filtered_feature_bc_matrix/")
ogd.data <- Read10X(data.dir = "/results_OGD/outs/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized) data.
CTRL <- CreateSeuratObject(counts = ctrl.data, project = "Control", min.cells = 3, min.features = 200)
CTRL$group <- "Control"
HYPO <- CreateSeuratObject(counts = hypo.data, project = "Hypoxia", min.cells = 3, min.features = 200)
HYPO$group <- "Hypoxia"
OGD <- CreateSeuratObject(counts = ogd.data, project = "OGD", min.cells = 3, min.features = 200)
OGD$group <- "OGD"

#Filtering based on mitochondrial percent and number of counts
CTRL[["percent.mt"]] <- PercentageFeatureSet(CTRL, pattern = "^MT-")
HYPO[["percent.mt"]] <- PercentageFeatureSet(HYPO, pattern = "^MT-")
OGD[["percent.mt"]] <- PercentageFeatureSet(OGD, pattern = "^MT-")

#Filter by QC Metrics
CTRL <- subset(CTRL, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & percent.mt < 15)
HYPO <- subset(HYPO, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & percent.mt < 15)
OGD <- subset(OGD, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & percent.mt < 15)

obj.list <- list(HYPO, OGD)
obj <- merge(CTRL, y = obj.list)

#Normalization (Using scTransform) and Integration
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- SCTransform(obj, vars.to.regress = "percent.mt")
obj <- RunPCA(obj)
obj <- IntegrateLayers(object = obj, method = CCAIntegration, normalization.method = "SCT")
ktest <- FindNeighbors(obj, dims = 1:40, k.param = 20, graph.name = c('graph1', 'graph2'))
ktest <- FindClusters(ktest, graph.name = 'graph2', resolution = 0.1, algorithm = 1)
ktest <- RunUMAP(ktest, reduction = "pca", dims = 1:40); DimPlot(ktest, reduction = "umap")
DimPlot(ktest, reduction = "umap", label = TRUE)
DefaultAssay(ktest) <- "RNA"
ktest <- JoinLayers(ktest)
rm('obj')

#Cluster Annotation
#Visualize the expression of cell type-specific markers
DefaultAssay(ktest) <- "RNA"

#Microglia/Macrophage Markers
VlnPlot(ktest, features = c('CX3CR1', 'HEXB', 'AIF1', 'PTPRC', 'NFKB1', 'CD86'))
#Stem Cell Markers
VlnPlot(ktest, features = c('FUT4', 'CD34', 'POU5F1', 'NANOG'))
#T-Cells
VlnPlot(ktest, features = c('CD3E', 'CD3D', 'CD3G', 'TRAC'))
#Others
VlnPlot(ktest, features = c('GFAP', 'MOG', 'FLT1', 'SNAP25'))

#Cluster Annotation
#Rename each cluster manually based on marker expression
new.cluster.ids <- c('Micro_1', 'Micro_2', 'Micro_3', 'Micro_4', 'Micro_5', 'Micro_6', 'Micro_7')
names(new.cluster.ids) <- levels(ktest)
ktest <- RenameIdents(ktest, new.cluster.ids)

#Save the clustered Seurat object
saveRDS(ktest, file = "~/iMGLFinalClustering.rds")
ktest <- readRDS("~/iMGLFinalClustering.rds")

#Figure S3A
DimPlot(ktest, reduction = "umap", label = TRUE)

#Figure S3B
manual.markers = c("CX3CR1", "P2RY12", "HEXB", "AIF1", "TGFB1", "TGFBR1", 'FUT4', 'CD34', 'POU5F1', 'NANOG', 'CD3E', 'GFAP', 'MOG', 'PCDH15', 'VWF', 'SNAP25')
VlnPlot(ktest, features = manual.markers, stack = TRUE, flip = TRUE) + NoLegend()

#Figure S3C
tbl <- table(ktest@active.ident, ktest@meta.data$group)
tbl <- as.data.frame(tbl)
ggplot(tbl, aes(fill = Var1, y = Freq, x = Var2)) + geom_bar(position="fill", stat="identity") + labs(x = "Sample", y = "Population Fraction", fill = "Cluster")

#Figure S3D
DimPlot(ktest, reduction = "umap", split.by = 'group')

#Figure S3E
if (!require('biomaRt')) BiocManager::install("biomaRt"); library('biomaRt')
ccmarkers <- read.csv("~/CellCycleGeneset.csv")
sgenes <- ccmarkers %>% dplyr::filter(phase == "S") %>% pull("geneID")
g2mgenes <- ccmarkers %>% dplyr::filter(phase == "G2/M") %>% pull("geneID")
DefaultLayer(ktest[["RNA"]]) <- 'data'
ktest = CellCycleScoring(ktest, g2m.features = g2mgenes, s.features = sgenes)
VlnPlot(ktest, features = "S.Score")

#Figure S3F
ktest <- SetIdent(ktest, value = "seurat_clusters")
new.cluster.ids <- c('Micro_1', 'Micro_2', 'Micro_3', 'Micro_4', 'Micro_5', 'Micro_6', 'Micro_7')
names(new.cluster.ids) <- levels(ktest)
ktest <- RenameIdents(ktest, new.cluster.ids)
if (!require('limma')) BiocManager::install('limma'); library('limma')
markers <- FindAllMarkers(ktest, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
tbl <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(as_tibble(tbl), n = 100)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
ktest.small <- subset(ktest, downsample = 1000)
DoHeatmap(ktest.small, features = top10$gene) + NoLegend()