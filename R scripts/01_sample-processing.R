# scRNAseq_CAascites Project
# Part 1 : sample preprocessing
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This script cover sample QC and preprocessing (by each library)
# Step 1 is running SoupX in order to correct ambient RNA
# Step 2 uses result from SoupX to create Seurat object and basic data QC and normalisation
# Step 3 uses DoubletFinder to flag potential doublets
# Then processed data is saved before proceed to next step (data integration)
# Please not that we didn't deposit cellranger output folders nor any individually processed files to GEO
# Examples of one output files can be found at XXXXXXXXXXXXXXXXX

# ------------------------------------------
# Load required libraries and set up basic parameters
# ------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(SoupX)

# set up basic parameter
wk_dir = "../sc_analysis/github/"
setwd(wk_dir)

sample_name = "ADC-Mech-Pre" # This is just a name, can be modified. I named this to match my cellranger output dir

save_dir = paste0(wk_dir,"rds/") # I will save processed/intermediate rds files in this dir
data_dir = paste0(wk_dir,"data/",sample_name) # Directory of 10x output files those caontain filtered_ and raw_ subdirectory


# ------------------------------------------
# Read data for SoupX
# ------------------------------------------

# read 10x files using soupx
scl <- load10X(data_dir, channelNames = sample_name)
colnames(scl[["toc"]]) <- paste0(sample_name,"___",colnames(scl[["toc"]]))

# I used fixed contamination fraction here because I experienced overcorrection if I follow reccommended steps.
scl = setContaminationFraction(scl, 0.05)
out = adjustCounts(scl)

# ------------------------------------------
# Seurat basic processing
# ------------------------------------------

# create Seurat object from corrected counts
scdata <- CreateSeuratObject(out, min.cells = 5, min.features = 200, project = sample_name)

scdata[["protocol"]] <- sample_name %>% str_split(pattern = "-", simplify = T) %>% .[,2]
scdata[["timepoint"]] <- sample_name %>% str_split(pattern = "-", simplify = T) %>% .[,3]
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
scdata[["percent.rp"]] <- PercentageFeatureSet(scdata, pattern = "^RP(S|L)")

cut.off.mito = 20 # set cut-off for %MT genes. I use 20% in this case
scdata <- subset(scdata, subset = nFeature_RNA > 200 & percent.mt < cut.off.mito)

scdata <- SCTransform(scdata, vars.to.regress = "percent.mt", assay = "RNA", verbose = FALSE)
scdata <- RunPCA(scdata) %>%
  RunUMAP(dims = 1:25) %>%
  FindNeighbors(dims = 1:25) %>%
  FindClusters(resolution = 0.6, algorithm = "leiden")

# Dimplot show cell clustered using leiden
# From this point we roughly get sense of clusters' identitdy by visualising some caninical genes e.g, EPCAM, SPARC, PTPRC using FeaturePlot or VlnPlot
DimPlot(scdata, label = TRUE, pt.size = 0.7) + labs(title = sample_name)

# ------------------------------------------
# DoubletFinder
# ------------------------------------------

library(DoubletFinder)

sweep.list <- paramSweep_v3(scdata, PCs = 1:25, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]
pK <- as.numeric(as.character(pK$pK))
homotypic.prop <- modelHomotypic(scdata@active.ident)
nExp_poi <- round(0.04*length(colnames(scdata)))  ## Assuming 4.0% doublet formation rate according to 10x protocol

scdata <- doubletFinder_v3(scdata, PCs = 1:25, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

DF_metadata <- paste("DF.classifications_0.25",pK,nExp_poi,sep="_")
table(scdata[[DF_metadata]])

# Dimplot shows doublets/singlets
DimPlot(scdata, group.by=DF_metadata, cols = c("black","red"), order = c("Doublet","Singlet"))

# Remove doublets from further downstream analysis
scdata <- SubsetData(scdata, ident.use = "Singlet")

# Repeat data processing steps after removing doublets
scdata <- SCTransform(scdata, vars.to.regress = "percent.mt", assay = "RNA", verbose = FALSE)
scdata <- RunPCA(scdata) %>%
  RunUMAP(dims = 1:25) %>%
  FindNeighbors(dims = 1:25) %>%
  FindClusters(resolution = 0.6, algorithm = "leiden")

DimPlot(scdata, label = TRUE, pt.size = 0.7) + labs(title = sample_name)

# ------------------------------------------
# Save
# ------------------------------------------

# save file for next step
saveRDS(scdata, file = paste0(save_dir,sample_name,"_DFDR.rds"))
