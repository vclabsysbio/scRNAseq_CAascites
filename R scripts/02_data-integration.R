# scRNAseq_CAascites Project
# Part 2 : data integration
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This shows how we integrate data before downstream analysis in our publication
# Please noted that processed file deposited on GEO are filtered cells and thus designed for part3 script

# ------------------------------------------
# Load required libraries and set up basic parameters
# ------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# set up basic parameter
wk_dir = "../sc_analysis/github/"
setwd(wk_dir)

# I will save processed/intermediate rds files in this dir
# This dir also contain processed files from part 1
save_dir = paste0(wk_dir,"rds/")

# ------------------------------------------
# Read processed data from part1
# ------------------------------------------

sample_list <- list("ADC-Mech-Pre","ADC-Enz-Pre","ADC-Mech-Post","ADC-Mech-Post") # list of names that we used for processing each lib

data_list <- list()

for (i in seq_along(sample_list)) {
  data_list[[i]] <- readRDS(paste0(save_dir,sample_list[[i]],"_DFDR.rds"))
  DefaultAssay(data_list[[i]]) <- "SCT"
}

data_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = data_features, verbose = FALSE)
data_anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = data_features, verbose = FALSE)

scdata <- IntegrateData(anchorset = data_anchors, normalization.method = "SCT", verbose = FALSE)

# factorise some metadata columns
scdata$timepoint <- factor(scdata$timepoint, levels = c("Pre","Post"))
scdata$orig.ident <- factor(scdata$orig.ident, levels = c("ADC-Mech-Pre","ADC-Mech-Post","ADC-Enz-Pre","ADC-Enz-Post"))
scdata$protocol <- factor(scdata$protocol, levels = c("Mech","Enz"))

# clean up some unused metadata
unused_meta <- c(grep("DF", colnames(scdata@meta.data), value = TRUE),
                 grep("pANN", colnames(scdata@meta.data), value = TRUE),
                 grep("SCT", colnames(scdata@meta.data), value = TRUE),
                 "seurat_clusters")
for (i in 1:length(unused_meta)) { scdata[[unused_meta[i]]] <- NULL}

# Normalise data on RNA assay
DefaultAssay(scdata) <- "RNA"
scdata <- NormalizeData(scdata, verbose = FALSE) %>% ScaleData(vars.to.regress = c("nCount_RNA"))

# Process data on integrated assay for visualisation
DefaultAssay(scdata) <- "integrated"
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:30, umap.method = "umap-learn")
scdata <- FindNeighbors(scdata, dims = 1:30)
scdata <- FindClusters(scdata, resolution = 0.6, algorithm = "leiden")

# Data visualusation
DimPlot(scdata, label = TRUE, pt.size = 0.1) + labs(title = "ADC-integrated") # Original leiden clustering
DimPlot(scdata, label = FALSE, pt.size = 0.1, group.by = "orig.ident") + labs(title = "ADC-integrated")
DimPlot(scdata, label = FALSE, pt.size = 0.1, split.by = "orig.ident", ncol = 2) + labs(title = "ADC-integrated")
DimPlot(scdata, label = FALSE, pt.size = 0.1, split.by = "protocol") + labs(title = "ADC-integrated")
DimPlot(scdata, label = FALSE, pt.size = 0.1, split.by = "timepoint") + labs(title = "ADC-integrated")

# Canonical gene markers for cluster annotation
DefaultAssay(scdata) <- "RNA"
FeaturePlot(scdata,
            features = c("EPCAM","KRT8","PTPRC","CD3E","SPARC","COL3A1","S100A8","CD14","NKG7"),
            cols = c("grey","red"),
            ncol = 3)

# ------------------------------------------
# Save
# ------------------------------------------

# save file for next step
saveRDS(scdata, file = paste0(save_dir,"ADC-integrated","_DFDR.rds"))
