# scRNAseq_CAascites (UNDER CONSTRUCTION)

## Description
This is a repository of scripts used for single-cell RNA-sequencing analysis in the publication entitled *"Capturing tumor heterogeneity in pre- and post-chemotherapy colorectal cancer ascites-derived cells using single-cell RNA-sequencing"*

## Data
Data used in this publication was deposited under ancession number GSEXXXXXXX

## Raw data processing
(Cell Ranger)[https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger] version 3.0.1 was used for reads mappings and UMI quantification using the following command

```
cellranger count --id=<sample-id> --transcriptome=<refdata-cellranger-GRCh38-1.2.0> --fastqs=<fastq_path> --expect-cells=5000 --localcores=36 --localmem=256
```

## R Analysis overview
The bioinformatic analyses of this work is divided into xxx steps:
  **1. Data preprocessing:** This part described how individual sample was QC and processed before data integration.
  **2. Data integration:** This part described how samples were integrated
  **3. Downstream analysis:** This part included downstream analysis in Fig. 2 & 3
  **4. Epithelial cells subset analysis:** This part shows epithelial cells subset analysis in Fig. 4 and 5A
  **5. Integrated with public dataset:** This part shows data integration with normal intestine single-cell dataset as shown in Fig.5B


## Packages used in this analysis

R 3.6.0
* Seurat (3.1.0)[https://github.com/satijalab/seurat]
* sctransform (0.2.0)[https://github.com/ChristophH/sctransform]
* SoupX (xxx)
* DoubletFinder (2.0.1)[https://github.com/chris-mcginnis-ucsf/DoubletFinder]
* ComplexHeatmap (2.0.0)[https://github.com/jokergoo/ComplexHeatmap]
* dplyr (0.8.3)[https://dplyr.tidyverse.org/]
* stringr (1.4.0)[https://stringr.tidyverse.org/]
* ggplot2 (3.2.0)[https://ggplot2.tidyverse.org/]
* ggrepel (0.8.1)[https://github.com/slowkow/ggrepel]
* dendextend (1.12.0)[https://github.com/talgalili/dendextend]
* fgsea (1.10.0)[https://github.com/ctlab/fgsea]
* reticulate (1.12)[https://rstudio.github.io/reticulate/]

Python 3.6.9
* numpy (1.17.4)[https://numpy.org/]
* umap-learn (0.3.10)[https://umap-learn.readthedocs.io/]
* leidenalg (0.7.0)[https://github.com/vtraag/leidenalg]
* python_igraph (0.7.1)[https://igraph.org/python/]


Data analyses were performed on Ubuntu 16.04.6 LTS
