library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(data.table)
library(glue)
library(DT)
library(stringr)
library(kableExtra)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(Seurat)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(plot1cell)
library(tidyverse)
library(purrr)
library(patchwork)
library(reshape2)
library(pheatmap)
library(Matrix)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
setwd("/Volumes/seqData/seq\ data/BR\ seq/code/20250617")

load("seurat.Rda")

DimPlot(seurat,group.by = "layer3")
marker_list <- list(
  "B Cells" = c("Cd19", "Ms4a1", "Cd79a"),
  "CD4 Naive" = c("Cd4", "Sell", "Tcf7"),
  "CD4 Tcm" = c("Il7r", "Ccr7", "Tcf7"),
  "CD4 Tem" = c("Il2ra", "Gzmk", "Cxcr3"),
  "CD4 Trm" = c("Cxcr6", "Itgae", "Zfp683"),
  "CD8 Effector" = c("Gzma", "Gzmb", "Prf1"),
  "CD8 Exhausted" = c("Tox", "Pdcd1", "Lag3"),
  "CD8 Naive" = c("Cd8a", "Sell", "Tcf7"),
  "CD8 Tem" = c("Il2rb", "Gzmk", "Cxcr3"),
  "DC" = c("Itgax", "H2-Aa", "Cd74"),
  "Macrophage" = c("Adgre1", "Cd68", "Mertk"),
  "Monocyte" = c("Ly6c2", "Ccr2", "Cd14"),
  "NK Cells" = c("Nkg7", "Klrd1", "Gzmb"),
  "Plasma" = c("Mzb1", "Sdc1", "Prdm1"),
  "Tfh" = c("Cxcr5", "Bcl6", "Pdcd1"),
  "Th1" = c("Tbx21", "Ifng", "Cxcr3"),
  "Treg" = c("Foxp3", "Il2ra", "Ctla4")
)

seen <- character()
unique_genes_ordered <- c()
for (cell in names(marker_list)) {
  for (gene in marker_list[[cell]]) {
    if (!(gene %in% seen)) {
      seen <- c(seen, gene)
      unique_genes_ordered <- c(unique_genes_ordered, gene)
    }
  }
}

cell_order <- names(marker_list)
Idents(seurat) <- "layer3"
seurat$layer3 <- factor(seurat$layer3, levels = cell_order)

dotdata_raw <- DotPlot(seurat, features = unique_genes_ordered, group.by = "layer3")$data

dotdata_raw$id <- factor(dotdata_raw$id, levels = cell_order)
dotdata_raw$features.plot <- factor(dotdata_raw$features.plot, levels = rev(unique_genes_ordered))

ggplot(dotdata_raw, aes(x = id, y = features.plot)) +
  geom_point(aes(size = avg.exp.scaled, color = pct.exp)) +
  scale_color_gradient(low = "grey90", high = "red") +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
