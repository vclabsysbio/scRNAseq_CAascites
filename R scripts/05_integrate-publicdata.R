# scRNAseq_CAascites Project
# Part 5 : integrate public data
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This part covers data analysis that generate figure 5B

# ------------------------------------------
# Load required libraries and set up basic parameters
# ------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

wk_dir = "../sc_analysis/github/"
setwd(wk_dir)

save_dir = paste0(wk_dir,"rds/")

# Downloaded GSE125970 processed data in this dir
data_dir = "../sc_analysis/GSE125970/"

# ------------------------------------------
# Read GSE data
# ------------------------------------------

GSE125970 <- read.table(paste0(data_dir,"GSE125970_raw_UMIcounts.txt.gz"), sep = "\t", header = T)
rownames(GSE125970) <- GSE125970$GENE
GSE125970$GENE <- NULL

GSE125970 <- CreateSeuratObject(GSE125970, project = "GSE125970")

GSE125970_cellinfo <- read.table(paste0(data_dir,"GSE125970_cell_info.txt.gz"), sep = "\t", header = T)
GSE125970$cell_anno <- GSE125970_cellinfo$CellType
GSE125970$sample <- GSE125970_cellinfo$Sample_ID
GSE125970$data <- "GSE125970"

GSE125970[["percent.mt"]] <- PercentageFeatureSet(GSE125970, pattern = "^MT-")

# ------------------------------------------
# Read scCC_Epi data
# ------------------------------------------

SCCC <- readRDS(file = paste0(save_dir,"ADC-integrated_Epi.rds"))

# ------------------------------------------
# Data integration
# ------------------------------------------

data_list <- c(SplitObject(GSE125970, split.by = "sample"),
               SplitObject(SCCC, split.by = "orig.ident"))

for (i in 1:length(data_list)) {
  DefaultAssay(data_list[[i]]) <- "RNA"
  suppressWarnings({
    data_list[[i]] <- SCTransform(data_list[[i]], vars.to.regress = "percent.mt", verbose = F)
  })
}

data_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = data_features, verbose = FALSE)
data_anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = data_features, verbose = FALSE)

scdata <- IntegrateData(anchorset = data_anchors, normalization.method = "SCT", verbose = FALSE)

DefaultAssay(scdata) <- "integrated"
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:30, n.neighbors = 30, min.dist = 0.05)

# ------------------------------------------
# Fig 5B
# ------------------------------------------

DimPlot(scdata,
        label = FALSE,
        pt.size = 0.1,
        group.by = "cell_anno",
        na.value = "lightgrey",
        split.by = "data") + labs(title = "Integrated data")

scdata$new.celltype.cluster <- factor(scdata$new.celltype.cluster, levels = str_sort(unique(scdata$new.celltype.cluster), numeric = TRUE))
DimPlot(scdata,
        label = FALSE,
        pt.size = 0.1,
        group.by = "new.celltype.cluster",
        na.value = "lightgrey",
        split.by = "data") + labs(title = "Integrated data")
