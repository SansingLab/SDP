if (!require("BiocManager")) install.packages("BiocManager")
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
if (!require("RRHO")) BiocManager::install("RRHO"); library('RRHO')

#Load Pre-Ranked Lists from EdgeR
list <- c("Oligodendrocytes", "Microglia", "Astrocytes", "OPCs")
for (x in list){
  filePath <- paste0("~/DEGs/DiazPerez/TStatRanks", x, "VaDvsCTRL.csv")
  a <- read.csv(filePath)
  new_name <- paste0("Sansing", x)
  assign(new_name, a)
}

#DAPI Pre-Ranked Lists
list <- c("Oligodendrocytes", "Astrocytes", "OPCs")
for (x in list){
  filePath <- paste0("~/DEGs/Mitroi/TStatRanks", x, "VaDvsCTRL.csv")
  a <- read.csv(filePath)
  new_name <- paste0("DAPI", x)
  assign(new_name, a)
}

#ERG Micro
ERGMicroglia <- read.csv("~/DEGs/Mitroi/TStatRanksMicrogliaVaDvsCTRL.csv")

#Combine respective lists between ERG and DAPI by: 
#1 Rename Freq column

complex <- list(SansingAstrocytes, SansingMicroglia, SansingOligodendrocytes, SansingOPCs, DAPIAstrocytes, ERGMicroglia, DAPIOligodendrocytes, DAPIOPCs)
rm(SansingAstrocytes, SansingMicroglia, SansingOligodendrocytes, SansingOPCs, DAPIAstrocytes, ERGMicroglia, DAPIOligodendrocytes, DAPIOPCs)
names <- c('SansingAstrocytes', 'SansingMicroglia', 'SansingOligodendrocytes', 'SansingOPCs', 'DAPIAstrocytes', 'ERGMicroglia', 'DAPIOligodendrocytes', 'DAPIOPCs')
names2 <- c('Sansing', 'Sansing', 'Sansing', 'Sansing', 'CAR', 'CAR', 'CAR', 'CAR')

for (x in 1:length(complex)){
  var_name <- names2[x]
  file_name <- names[x]
  names(complex[[x]])[names(complex[[x]]) == 'Freq'] <- var_name
  assign(file_name, as.data.frame(complex[[x]]))
}

rm(complex)

#2 Removing rows where genes are not shared

Astros <- merge(SansingAstrocytes, DAPIAstrocytes, by = "Var1")
Micro <- merge(SansingMicroglia, ERGMicroglia, by = "Var1")
Oligos <- merge(SansingOligodendrocytes, DAPIOligodendrocytes, by = "Var1")
OPCs <- merge(SansingOPCs, DAPIOPCs, by = "Var1")
rm(SansingAstrocytes, SansingMicroglia, SansingOligodendrocytes, SansingOPCs, DAPIAstrocytes, ERGMicroglia, DAPIOligodendrocytes, DAPIOPCs)

#3 Add rows ranking the genes based on T-Stat values of ERG and DAPI
#4 Save as csv and download for use with web-based app: https://systems.crump.ucla.edu/rankrank/

complex <- list(Astros, Micro, Oligos, OPCs)
names <- c('Astros', 'Micro', 'Oligos', 'OPCs')
for (x in 1:length(complex)){
  newdf <- complex[[x]]
  newdf <- newdf[order(newdf$Sansing, decreasing = TRUE),]
  newdf$rankSansing <- NA
  newdf$rankSansing <- 1:nrow(newdf)
  newdf <- newdf[order(newdf$CAR, decreasing = TRUE),]
  newdf$rankCAR <- NA
  newdf$rankCAR <- 1:nrow(newdf)
  newdf <- newdf[,c(1, 4, 5, 2, 3)]
  file_name <- names[x]
  assign(file_name, newdf)
}

###################Implement RRHO
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("RRHO")) BiocManager::install("RRHO")
library(RRHO)

k <- list(Astros, Micro, Oligos, OPCs)
names <- c('Astros', 'Micro', 'Oligos', 'OPCs')

for (x in 1:length(k)){
  df1 <- k[[x]][,c(1,4)]
  dfname1 <- paste0(names[x], "1")
  assign(dfname1, df1)
  df2 <- k[[x]][,c(1,5)]
  dfname2 <- paste0(names[x], "2")
  assign(dfname2, df2)
  dir = paste0("~/RRHO/", names[x])
  labels = c(colnames(df1)[2], colnames(df2)[2])
  RRHO(df1, df2, alternative = "two.sided", labels = labels, plots = TRUE, outputdir = dir, BY = TRUE, log10.ind = TRUE)
}