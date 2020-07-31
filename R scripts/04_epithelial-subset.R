# scRNAseq_CAascites Project
# Part 4 : Epithelial cells subset analysis
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This part covers
# 1. how to start from deposited processed files on GEO
# 2. Data analysis that generate figure 4 and 5A

# ------------------------------------------
# Load required libraries and set up basic parameters
# ------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

wk_dir = "../sc_analysis/github/"
setwd(wk_dir)

save_dir = paste0(wk_dir,"rds/")
data_dir = paste0(wk_dir,"data/ADC-integrated_Epi.txt")
metadata_dir = paste0(wk_dir,"data/ADC-integrated_Epi_cellinfo.txt")

# ------------------------------------------
# Read processed data
# ------------------------------------------

data <- read.delim(data_dir)
colnames(data) <- colnames(data) %>% str_replace_all(pattern = "\\.", replacement = "-")
metadata <- read.delim(metadata_dir)
metadata$celltype <- metadata$new.celltype.cluster

scdata <- CreateSeuratObject(data, meta.data = metadata)
head(scdata@meta.data)

# split object by original identity
data_list <- SplitObject(scdata, split.by = "orig.ident")

# data normalisation
for (i in seq_along(data.list)) {
  DefaultAssay(data_list[[i]]) <- "RNA"
  data_list[[i]][["percent.mt"]] <- PercentageFeatureSet(data_list[[i]], pattern = "^MT-")
  data_list[[i]] <- SCTransform(data_list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}

# data integration
data_features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = data_features, verbose = FALSE)
data_anchors <- FindIntegrationAnchors(object.list = data_list, normalization.method = "SCT", anchor.features = data_features, verbose = FALSE)

scdata <- IntegrateData(anchorset = data_anchors, normalization.method = "SCT", verbose = FALSE)

# factorise some metadata columns
scdata$timepoint <- factor(scdata$timepoint, levels = c("Pre","Post"))
scdata$orig.ident <- factor(scdata$orig.ident, levels = c("ADC-Mech-Pre","ADC-Mech-Post","ADC-Enz-Pre","ADC-Enz-Post"))
scdata$protocol <- factor(scdata$protocol, levels = c("Mech","Enz"))

# Normalise data on RNA assay
DefaultAssay(scdata) <- "RNA"
scdata <- NormalizeData(scdata, verbose = FALSE) %>% ScaleData(vars.to.regress = c("nCount_RNA"))

# Process data on integrated assay for visualisation
DefaultAssay(scdata) <- "integrated"
scdata <- RunPCA(scdata)
scdata <- RunUMAP(scdata, dims = 1:30, umap.method = "umap-learn")
scdata <- FindNeighbors(scdata, dims = 1:30)
scdata <- FindClusters(scdata, resolution = 0.4, algorithm = "leiden")

DimPlot(scdata, label = TRUE, pt.size = 0.1)

# ------------------------------------------
# Save
# ------------------------------------------

# save file
saveRDS(scdata, file = paste0(save_dir,"ADC-integrated_Epi.rds"))

# ------------------------------------------
# Fig 4A
# ------------------------------------------

DimPlot(scdata, label = TRUE, pt.size = 0.1)

# ------------------------------------------
# Fig 4B
# ------------------------------------------

library(dendextend)

Idents(scdata) <- "celltype"

pop_tb <- table(scdata$new.celltype.cluster,scdata$orig.ident)
pop_tb <- as.data.frame.matrix(pop_tb)
pop_tb$new.celltype.cluster <- rownames(pop_tb)

pop_plot <- tidyr::gather(pop_tb, key = sample, value = counts, -new.celltype.cluster)

dend <- as.dendrogram(hclust(dist(with(pop_plot, tapply(counts, new.celltype.cluster, mean)))))
dend %>% plot()
dend %>% plot(leaflab = "none", edgePar = list(lwd = 2))
plot_order <- dend %>% labels

pop_plot$new.celltype.cluster <- factor(pop_plot$new.celltype.cluster, levels = plot_order)
pop_plot$sample <- factor(pop_plot$sample, levels = c("ADC-Mech-Pre","ADC-Enz-Pre","ADC-Mech-Post","ADC-Enz-Post"))

pop_plot <- pop_plot %>%
  group_by(new.celltype.cluster,sample) %>%
  summarise(counts = sum(counts)) %>%
  mutate(percent = counts*100/sum(counts)) %>%
  mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
  mutate(position = 100-(cumsum(percent) - (0.5*percent)))

ggplot(pop_plot, aes(x= new.celltype.cluster, y = percent)) +
  geom_bar(aes(fill = sample),stat = "identity") +
  geom_text(aes(label = percent.text, y = position)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(label = "Percentage") + xlab(label = "") +
  scale_fill_manual(values = c("lightsalmon","lightsalmon3","lightskyblue2","skyblue3"))

# ------------------------------------------
# Fig 4C
# ------------------------------------------

scdata$new.celltype.cluster <- factor(scdata$new.celltype.cluster, levels = str_sort(unique(scdata$new.celltype.cluster), numeric = TRUE))
Idents(scdata) <- "new.celltype.cluster"

scdata_all_markers <- FindAllMarkers(scdata, assay = "RNA", verbose = FALSE)

scdata_all_markers <- scdata_all_markers %>%
  mutate(pct.diff = pct.1-pct.2) %>%
  mutate(p_val_adj_spf = sprintf("%2.2e",p_val_adj)) %>%
  mutate(log10p_val = -log10(p_val)) %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(log2_avgFC = log2(expm1(avg_logFC)+1))

min_pct <- 0.25
scdata_top5_markers <- scdata_all_markers %>% filter(pct.1 >= min_pct) %>% group_by(cluster) %>% top_n(5, log2_avgFC)


Idents(scdata) <- "new.celltype.cluster"
avg.data <- AverageExpression(scdata, verbose = FALSE, assays = "RNA")

emat <- avg.data$RNA %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(gene %in% scdata_top5_markers$gene) %>%
  arrange(factor(gene, levels = scdata_top5_markers$gene, labels = scdata_top5_markers$gene))
rownames(emat) <- emat$gene
emat$gene <- NULL
emat <- log1p(emat)
custom_row_order <- intersect(scdata_top5_markers$gene, rownames(emat))

emat <- t(scale(t(emat))) # scale and center rows

Heatmap(emat,
        name = "row z score",
        col = circlize::colorRamp2(c(-2, 0, 2), c("magenta", "black", "yellow")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 0.2),
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_side = "top",
        row_order = custom_row_order
)

# ------------------------------------------
# Fig 5A
# ------------------------------------------

library(fgsea)
set.seed(42)

# Load MSigDB geneset. Can be found at https://www.gsea-msigdb.org/
pathway_gmt <- gmtPathways("../msigdb/h.all.v6.2.symbols.gmt")

Idents(scdata) <- "new.celltype.cluster"

# This can take long
scdata_all_markers <- FindAllMarkers(scdata, assay = "RNA", min.pct = 0, logfc.threshold = 0)
# saveRDS(scdata_all_markers, file = paste0(save_dir,"ADC-integrated_Epi_allmarkers.rds"))

scdata_all_markers <- readRDS(file = paste0(save_dir,"ADC-integrated_Epi_allmarkers.rds"))

cluster_list <- levels(scdata@active.ident)

hallmark_nes_df <- data.frame("pathway" = names(pathway_gmt))
hallmark_padj_df <- data.frame("pathway" = names(pathway_gmt))

for (n in seq_along(levels(scdata@active.ident))) {
  cluster.name <- cluster_list[n]
  ranks <- scdata_all_markers %>%
    filter(cluster==cluster.name) %>%
    dplyr::select(gene, avg_logFC) %>%
    na.omit() %>%
    distinct() %>%
    group_by(gene) %>%
    tibble::deframe() %>% sort(decreasing = TRUE)
  fgseaRes <- fgsea(pathways=pathway_gmt, stats=ranks, maxSize = 500, nperm=1000)

  tmp.padj <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme, -pval, -size) %>%
    arrange(padj) %>%
    dplyr::select(pathway, padj)
  colnames(tmp.padj) <- c("pathway",cluster.name)
  hallmark_padj_df <- hallmark_padj_df %>% full_join(tmp.padj)

  tmp.nes <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme, -pval, -size) %>%
    arrange(padj) %>%
    dplyr::select(pathway, NES)
  colnames(tmp.nes) <- c("pathway",cluster.name)
  hallmark_nes_df <- hallmark_nes_df %>% full_join(tmp.nes)

}

hallmark_padj_df <- hallmark_padj_df %>% filter_all(any_vars(. < 0.05))
hallmark_nes_df <- hallmark_nes_df %>% filter(pathway %in% hallmark_padj_df$pathway)

rownames(hallmark_nes_df) <- hallmark_nes_df$pathway
hallmark_nes_df$pathway <- NULL

rownames(hallmark_padj_df) <- hallmark_padj_df$pathway
hallmark_padj_df$pathway <- NULL

Heatmap(hallmark_nes_df,
        name = "NES",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        rect_gp = gpar(col = "white", lwd = 1),
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(hallmark_padj_df[i,j]) & hallmark_padj_df[i,j] < 0.05 & hallmark_nes_df[i,j] > 0){
            grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 16, col = "white", fontface = "bold"))
          }
        },
        row_names_max_width = max_text_width(rownames(hallmark_nes_df)),
        row_names_gp = gpar(fontsize = 10)
)
