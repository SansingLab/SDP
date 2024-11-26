if (!require("BiocManager", quietly = TRUE)) BiocManager::install(version = "3.18")
library(scater)
if (!require('Seurat')) remotes::install_version('Seurat', version='5.0.0'); library(Seurat)
library(tidyverse)
library(cowplot)
if (!require(Matrix.utils)) remotes::install_github("cvarrichio/Matrix.utils"); library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(purrr)

#Load Clustered Seurat Object
seurat <- readRDS("~/FinalClustering.rds")
seurat$active.ident <- seurat@active.ident

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat, layer = "counts", assay = "RNA")
metadata <- seurat@meta.data
Idents(object = seurat) <- "SCT_snn_res.0.1"
metadata$cluster_id <- metadata$active.ident
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

## Remove lowly expressed genes which have less than 10 cells with any counts
dim(sce)
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
nk <- length(kids)
sids <- purrr::set_names(levels(as.factor(sce$orig.ident)))
ns <- length(sids)

# Generate sample level metadata
table(sce$orig.ident)
n_cells <- table(sce$orig.ident) %>%  as.vector()
names(n_cells) <- names(table(sce$orig.ident))
m <- match(names(n_cells), sce$orig.ident)
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("orig.ident", "group", "n_cells")
kable(ei)

# Aggregate the counts per sample_id and cluster_id
groups <- colData(sce)[, c("cluster_id", "orig.ident")]
groups$orig.ident <- factor(groups$orig.ident)
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
pb[1:8, 1:8]
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",n = 2), `[`, 1)
pb <- split.data.frame(pb,factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_", "", rownames(u))))

# Print out the table of cells in each cluster-sample group
options(width = 100)
kable(table(sce$cluster_id, sce$orig.ident))

#Filter out any Clusters where cells per sample < 30 across multiple samples
keepClusters <- as.character(unique(sce$cluster_id))
(interestingClusters <- SingleCellExperiment(assays = pb[keepClusters]))

# compute MDS
mds <-  lapply(as.list(assays(interestingClusters)), function(a){
  DGEList(a, remove.zeros = TRUE) %>% 
    calcNormFactors %>% 
    plotMDS.DGEList(plot = FALSE)
})
cnames <- paste("Cluster", keepClusters) 
for (m in 1:length(mds)){
  mds[[m]]$cluster <- cnames[m]
}
plots <- lapply(mds, function(m){
  gg_df <- data.frame(m[c("x", "y")],
                      sample_id = sids,
                      group_id = ei$group,
                      cluster_id = rep(m$cluster, length(m$x)))})
plotFunc <- function(x) {
  ggplot(x, aes(x, y, col = group_id)) + 
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
    ggtitle(unique(x$cluster_id)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed() 
}
do.call(grid.arrange,c(lapply(plots, plotFunc)))

#Voom LM Fit Model
(design <- model.matrix(~ 0 + ei$group) %>% 
    set_rownames(ei$orig.ident) %>% 
    set_colnames(levels(factor(ei$group))))

(contrast <- makeContrasts("VaD-Control", levels = design))

min.count = 0.1
total.count = 15

dflist <- list()
fit <- list()
efit <- list()
lcpm <- list()
for (i in 1:length(keepClusters)) {
  dflist[[i]] = DGEList(assays(interestingClusters)[[i]])
  dflist[[i]] = calcNormFactors(dflist[[i]], method = "RLE")
  genes.use1 = filterByExpr(dflist[[i]], design = design, min.count = min.count, min.total.count = total.count)
  dflist[[i]] = dflist[[i]][ genes.use1, keep.lib.sizes=FALSE]
  fit[[i]] <- voomLmFit(dflist[[i]],  design, plot=FALSE, sample.weights=TRUE)
  fit[[i]] <- contrasts.fit(fit[[i]], contrast)
  efit[[i]] <- eBayes(fit[[i]])
}

require(limma)
ranks = list()
for (x in 1:length(efit)) {
    test = as.data.frame(topTable(efit[[x]], coef = NULL, number = Inf))
    ranks[[x]] =
      test %>% rownames_to_column("gene") %>%
      arrange(desc(t)) %>%
      select(gene, t) %>%
      column_to_rownames("gene") %>%
      t() %>%
      unlist(use.names = T)
    ranks[[x]] = ranks[[x]][1, ]
}

#DGE Analysis
res <- lapply(keepClusters, function(k) {
  y <- assays(interestingClusters)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  keep <- filterByExpr(y, design = design, group = NULL, lib.size = NULL, min.count = min.count, min.total.count = total.count)
  y <- y[keep,]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, contrast = contrast)
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})

#Write results as CSVs:
for(cluster in 1:length(keepClusters)){
  filePath <- paste0("~/TStatRanks", keepClusters[cluster])
  out <- as.table(ranks[[cluster]])
  write.csv(out, file = paste0(filePath, "VaDvsCTRL.csv"), quote=F, row.names = F)
}