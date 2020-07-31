# scRNAseq_CAascites Project
# Part 3 : DE analysis
# version 1.0 (July 2020)
# Tiraput Poonpanichakul, Faculty of Science, Mahidol University

# ------------------------------------------
# README
# ------------------------------------------
# This part covers
# 1. how to start from deposited processed files on GEO
# 2. Data analysis that generate figure 2 and 3

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
data_dir = paste0(wk_dir,"data/ADC-integrated.txt")
metadata_dir = paste0(wk_dir,"data/ADC-integrated_cellinfo.txt")

# ------------------------------------------
# read files and create seurat object
# ------------------------------------------

data <- read.delim(data_dir)
colnames(data) <- colnames(data) %>% str_replace_all(pattern = "\\.", replacement = "-")
metadata <- read.delim(metadata_dir)
metadata$celltype_cluster <- paste0(metadata$celltype,"_",metadata$seurat_clusters)

scdata <- CreateSeuratObject(data, meta.data = metadata)
head(scdata@meta.data)

# split into objects by original sample
data_list <- SplitObject(scdata, split.by = "orig.ident")

# data normalisation
for (i in seq_along(data_list)) {
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

# ------------------------------------------
# Save
# ------------------------------------------

# save file
saveRDS(scdata, file = paste0(save_dir,"ADC-integrated.rds"))
# scdata <- readRDS(file = paste0(save_dir,"ADC-integrated.rds"))

# ------------------------------------------
# Fig 2A
# ------------------------------------------

plot_title = "ADC-integrated"

table(scdata$celltype, scdata$orig.ident) #freq by celltype x samples
DimPlot(scdata, label = TRUE, pt.size = 0.7, group.by = "celltype") + labs(title = plot_title)
DimPlot(scdata, label = TRUE, pt.size = 0.7, group.by = "celltype_cluster") + labs(title = plot_title)

# ------------------------------------------
# Fig 2B
# ------------------------------------------

DefaultAssay(scdata) <- "RNA"
FeaturePlot(scdata,
            features = c("EPCAM","KRT8","PTPRC","CD3E","SPARC","COL3A1","S100A8","CD14","NKG7"),
            cols = c("grey","red"),
            pt.size = 0.1,
            ncol = 3
            )

FeaturePlot(scdata,
            features = c("PTPRC","CD3E","CD79A","NKG7"),
            cols = c("grey","red"),
            pt.size = 0.1,
            ncol = 2
)

# ------------------------------------------
# Fig 2C
# ------------------------------------------

library(ComplexHeatmap)

Idents(scdata) <- "celltype"
scdata_all_markers <- FindAllMarkers(scdata, assay = "RNA", verbose = FALSE)

scdata_all_markers <- scdata_all_markers %>%
  mutate(pct.diff = pct.1-pct.2) %>%
  mutate(p_val_adj_spf = sprintf("%2.2e",p_val_adj)) %>%
  mutate(log10p_val = -log10(p_val)) %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(log2_avgFC = log2(expm1(avg_logFC)+1))

min_pct <- 0.25
scdata_top20_markers <- scdata_all_markers %>% filter(pct.1 >= min_pct) %>% group_by(cluster) %>% top_n(20, log2_avgFC)


Idents(scdata) <- "celltype"
avg.data <- AverageExpression(scdata, verbose = FALSE, assays = "RNA")

emat <- avg.data$RNA %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(gene %in% scdata_top20_markers$gene) %>%
  arrange(factor(gene, levels = scdata_top20_markers$gene, labels = scdata_top20_markers$gene))
rownames(emat) <- emat$gene
emat$gene <- NULL
emat <- log1p(emat)
custom_row_order <- intersect(scdata_top20_markers$gene,rownames(emat))

emat <- t(scale(t(emat))) # scale and center rows

Heatmap(emat,
        name = "row z score",
        col = circlize::colorRamp2(c(-2, 0, 2), c("magenta", "black", "yellow")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 0.2),
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_side = "top",
        row_order = custom_row_order,
)

# ------------------------------------------
# Fig 3A
# ------------------------------------------

pop_tb <- table(scdata$orig.ident,scdata$celltype)
pop_tb <- as.data.frame.matrix(pop_tb)
pop_tb$sample <- rownames(pop_tb)

pop_plot <- tidyr::gather(pop_tb, key = celltype, value = counts, -sample)

dend <- as.dendrogram(hclust(dist(with(pop_plot, tapply(counts, sample, mean)))))
dend %>% plot()
dend %>% plot(leaflab = "none", edgePar = list(lwd = 2))
plot_order <- dend %>% labels

pop_plot$sample <- factor(pop_plot$sample, levels = plot_order)

pop_plot <- pop_plot %>%
  group_by(sample,celltype) %>%
  summarise(counts = sum(counts)) %>%
  mutate(percent= counts*100/sum(counts)) %>%
  mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
  mutate(position = 100-(cumsum(percent) - (0.5*percent)))

ggplot(pop_plot, aes(x= sample, y = percent)) +
  geom_bar(aes(fill = celltype),stat = "identity", width = 0.5, color = "darkgrey") +
  ggrepel::geom_label_repel(aes(label = percent.text, y = position, colour = celltype),
                            size = 6,
                            nudge_x = 0.3,
                            direction = "y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab(label = "Percentage") + xlab(label = "")

# ------------------------------------------
# Fig 3C
# ------------------------------------------

Idents(scdata) <- "celltype"
ident_subset = "Epi"
scdata <- subset(scdata, idents = ident_subset)

Idents(scdata) <- "protocol"

min_pct <- 0.25
logfc_threshold <- 0

MvE_markers <- FindMarkers(scdata, ident.1 = "Mech", ident.2 = "Enz",
                           min.pct = min_pct, logfc.threshold = logfc_threshold, assay = "RNA"
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(pct.diff = pct.1-pct.2) %>%
  arrange(-avg_logFC) %>%
  mutate(cluster = "MvE") %>%
  mutate(p_val_adj_spf = sprintf("%2.2e",p_val_adj)) %>%
  mutate(log10p_val = -log10(p_val)) %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(log2_avgFC = log2(expm1(avg_logFC)+1))

#sort table based on fold change and adj p val
MvE_markers <- MvE_markers %>%
  arrange(-log2_avgFC*log10p_val_adj)

#factor threshold criteria
avglfc_threshold <- 0.5
padj_threshold <- 0.001
MvE_markers$threshold = as.factor(abs(MvE_markers$log2_avgFC) > avglfc_threshold & MvE_markers$p_val_adj < padj_threshold)

#label gene name for top abs(lof change) that were sorted priori
ngenes_label <- 15
last_row <- nrow(MvE_markers)
MvE_markers$genelabel = ""
MvE_markers[1:ngenes_label,]$genelabel <- MvE_markers[1:ngenes_label,]$gene
MvE_markers[(last_row-ngenes_label):last_row,]$genelabel <- MvE_markers[(last_row-ngenes_label):last_row,]$gene

#find max x value in order to create balanced x axis plot
max_axis_x <- MvE_markers$log2_avgFC %>% abs() %>% max()

#volcano plot
ggplot(MvE_markers, aes(x=log2_avgFC, y=log10p_val_adj, label = genelabel)) +
  geom_vline(xintercept = c(-avglfc_threshold,avglfc_threshold), colour="lightgrey", linetype = "longdash") +
  geom_hline(yintercept = -log10(padj_threshold), colour="lightgrey", linetype = "longdash") +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(color = threshold)) +
  coord_cartesian(xlim = c(-max_axis_x,max_axis_x)) +
  scale_colour_manual(values=c("TRUE"="red","FALSE"="grey40")) +
  xlab("log2_avgFC") + ylab("-log10 p val adj") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle(label = "DEGs: Mech vs Enz")


# ------------------------------------------
# Fig 3D
# ------------------------------------------

Idents(scdata) <- "protocol"
DefaultAssay(scdata) <- "RNA"

enz_bc <- WhichCells(scdata, idents = "Enz", slot = "count")
mech_bc <- WhichCells(scdata, idents = "Mech", slot = "count")

df <- data.frame(gene = character(),
                 protocol = character(),
                 expr = character(),
                 pct = numeric(),
                 mean.exp = numeric(),
                 mean.all = numeric()
)

for (gene in c("HES1","IER3","ID1","CD81","PPP1CB")) {

  mean.all <- FetchData(object = scdata, vars = gene, slot = "data") %>% filter(.data[[gene]] > 0) %>% .[[gene]] %>% expm1() %>% quantile(.,0.25) %>% log1p()

  mean.enz <- FetchData(object = scdata, cells = enz_bc, vars = gene, slot = "data") %>% filter(.data[[gene]] > 0) %>% .[[gene]] %>% ExpMean()
  mean.mech <- FetchData(object = scdata, cells = mech_bc, vars = gene, slot = "data") %>% filter(.data[[gene]] > 0) %>% .[[gene]] %>% ExpMean()

  enz <- FetchData(object = scdata, cells = enz_bc, vars = gene, slot = "data") %>% filter(.data[[gene]] > mean.all) %>% nrow()
  mech <- FetchData(object = scdata, cells = mech_bc, vars = gene, slot = "data") %>% filter(.data[[gene]] > mean.all) %>% nrow()

  pct.enz <- enz/length(enz_bc)*100
  pct.mech <- mech/length(mech_bc)*100

  df <- rbind(df, data.frame(gene = gene,
                             protocol = "Enz",
                             expr = "exp",
                             pct = pct.enz,
                             mean.exp = mean.enz,
                             mean.all = mean.all)
  )

  df <- rbind(df, data.frame(gene = gene,
                             protocol = "Enz",
                             expr = "nexp",
                             pct = 100-pct.enz,
                             mean.exp = mean.enz,
                             mean.all = mean.all)
  )

  df <- rbind(df, data.frame(gene = gene,
                             protocol = "Mech",
                             expr = "exp",
                             pct = pct.mech,
                             mean.exp = mean.mech,
                             mean.all = mean.all)
  )

  df <- rbind(df, data.frame(gene = gene,
                             protocol = "Mech",
                             expr = "nexp",
                             pct = 100-pct.mech,
                             mean.exp = mean.mech,
                             mean.all = mean.all)
  )

}


ggplot(df, aes(x = protocol, y= pct, fill = expr)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(gene)) +
  scale_fill_manual(values=c("red3", "darkgrey")) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Preparation protocol") + ylab("Percent") +
  ggtitle("% cells expression > Q25")

# ------------------------------------------
# Fig 3E
# ------------------------------------------

Idents(scdata) <- "timepoint"

min_pct <- 0.25
logfc_threshold <- 0
avglfc_threshold <- 0.5
padj_threshold <- 0.001

pre_post_markers <- FindMarkers(scdata, ident.1 = "Pre", ident.2 = "Post",
                             min.pct = min_pct, logfc.threshold = logfc_threshold, assay = "RNA"
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(pct.diff = pct.1-pct.2) %>%
  arrange(-avg_logFC) %>%
  mutate(cluster = "Pre_Post") %>%
  mutate(p_val_adj_spf = sprintf("%2.2e",p_val_adj)) %>%
  mutate(log10p_val = -log10(p_val)) %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(log2_avgFC = log2(expm1(avg_logFC)+1)) %>%
  mutate(threshold = as.factor(abs(log2_avgFC) > avglfc_threshold & p_val_adj < padj_threshold))

max_row <- nrow(pre_post_markers)
DE_table <- rbind(pre_post_markers[1:20,],
                  pre_post_markers[(max_row-20):max_row,])

Idents(scdata) <- "orig.ident"
avg_scdata <- AverageExpression(scdata, assays = "RNA")

emat <- avg_scdata$RNA %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(gene %in% DE_table$gene) %>%
  arrange(factor(gene, levels = DE_table$gene))
rownames(emat) <- emat$gene
emat$gene <- NULL
emat <- log1p(emat)

emat <- t(scale(t(emat))) # scale and center rows
Heatmap(emat,
        name = "log1p expr\n(row z-score)",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lwd = 0.5),
        column_title_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_rot  = 45,
        column_names_side = "top",
        column_order = c("ADC-Mech-Pre","ADC-Enz-Pre","ADC-Mech-Post","ADC-Enz-Post")
)
