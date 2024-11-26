if (!require("BiocManager")) install.packages("BiocManager")
library(scater)
if (!require('Seurat')) install.packages('Seurat'); library(Seurat)
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
if (!require("RRHO")) BiocManager::install("RRHO"); library('RRHO')

#For iPSC comparison we use FC instead of TStats
#Load Pre-Ranked Lists
SansingMicroglia <- read.csv("~/DEGs/WilcoxonDiazPerezMicroglia.csv")
SansingMicroglia <- SansingMicroglia[SansingMicroglia$p_val_adj <= 0.05,]
SansingMicroglia <- SansingMicroglia[,c(1, 3)]
names(SansingMicroglia)[names(SansingMicroglia) == 'avg_log2FC'] <- 'Sansing'
ERGMicroglia <- read.csv("~/DEGs/WilcoxonMitroiMicroglia.csv")
ERGMicroglia <- ERGMicroglia[ERGMicroglia$p_val_adj <= 0.05,]
ERGMicroglia <- ERGMicroglia[,c(1, 3)]
names(ERGMicroglia)[names(ERGMicroglia) == 'avg_log2FC'] <- 'CAR'
iMGL <- read.csv("~/DEGs/WilcoxoniMGL.csv")
iMGL <- iMGL[iMGL$p_val_adj <= 0.05,]
iMGL <- iMGL[,c(1, 3)]
names(iMGL)[names(iMGL) == 'avg_log2FC'] <- 'iMGL'

#Combine lists:
SvM <- merge(SansingMicroglia, iMGL, by = "X")
CvM <- merge(ERGMicroglia, iMGL, by = "X")
SvC <- merge(SansingMicroglia, ERGMicroglia, by = "X")
rm(iMGL, ERGMicroglia, SansingMicroglia)

SvM <- SvM[order(SvM$Sansing, decreasing = TRUE),]
SvM$rankSansing <- 1:nrow(SvM)
SvM <- SvM[order(SvM$iMGL, decreasing = TRUE),]
SvM$rankiMGL <- 1:nrow(SvM)
SvM <- SvM[,c(1, 4, 5, 2, 3)]

CvM <- CvM[order(CvM$CAR, decreasing = TRUE),]
CvM$rankCAR <- 1:nrow(CvM)
CvM <- CvM[order(CvM$iMGL, decreasing = TRUE),]
CvM$rankiMGL <- 1:nrow(CvM)
CvM <- CvM[,c(1, 4, 5, 2, 3)]

SvC <- SvC[order(SvC$Sansing, decreasing = TRUE),]
SvC$rankSansing <- 1:nrow(SvC)
SvC <- SvC[order(SvC$CAR, decreasing = TRUE),]
SvC$rankCAR <- 1:nrow(SvC)
SvC <- SvC[,c(1, 4, 5, 2, 3)]


###################Implement RRHO
k <- list(CvM, SvC, SvM)
names <- c('CARvsiMGL', 'SansingvsCAR', 'SansingvsiMGL')

for (x in 1:length(k)){
  df1 <- k[[x]][,c(1,4)]
  dfname1 <- paste0(names[x], "1")
  assign(dfname1, df1)
  df2 <- k[[x]][,c(1,5)]
  dfname2 <- paste0(names[x], "2")
  assign(dfname2, df2)
  dir = paste0("~/RRHO/", names[x])
  labels = c(colnames(df1)[2], colnames(df2)[2])
  RRHO(df1, df2, alternative = "two.sided", plots = TRUE, outputdir = dir, labels = labels, BY = TRUE, log10.ind = TRUE)
}