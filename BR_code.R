# libraries from CRAN-------------------------------------------
library(here)
library(SingleR)
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
library(reticulate)
library(SeuratDisk)
library(sceasy)
mouseRNA <- MouseRNAseqData()
mouseImmu <- ImmGenData()
scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")
ad <- import("anndata")
loompy <- reticulate::import('loompy')
setwd("~/seq data/BR seq/")

#read 10x-------------------------------------------
seurat_data <- Read10X("/LBYY-20241223-ScRNA-BR-SH/Gq-BIBN4096.matrix/")
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                 project = "BR",
                                 min.features = 200,
                                 min.cells = 3)
sce <- NormalizeData(seurat_obj, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
rm(seurat_data,seurat_obj)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce),verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:17)
sce <- FindClusters(sce, resolution = 0.7)  
sce <- RunTSNE(sce, dims = 1:17)
DimPlot(sce, reduction = "umap",label=T)

# include the short version of the most recent git repository SHA
save(seurat, file = file.path("/LBYY-20241223-ScRNA-BR-SH/Gq-BIBN4096.matrix/", "seurat.Rda"))

#inte---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load libraries
library(here)
library(readr)
library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(ggplot2)
library(Seurat)
library(harmony)

set.seed(100)

here::here()

dir.create("figures")
dir.create("output")

message("@ loading data...")

setwd("~/BR seq/code")

# Prepare data -----------------------------------------------------------------
# 1. load the individual objects
# 2. apply the sample-group association to the metadata
# 3. apply specified filters as they're loaded
message("@ loading Seurat objects...")
# desired behaviour: if the processed joint object already exists, do not
# regenerate it, to allow for cases where the job might fail in the second half
# of the script due to mem/time limits, and allow for restarting part way
if (file.exists("output/seurat_joint.Rda")) {
  
  message("@ found joint object at output/seurat_joint.Rda\n",
          "@ loading instead of recalculating...\n",
          "@ to recalculate, clear old results by deleting output/figures folders")
  
  load("output/seurat_joint.Rda")
  
} else {
  
  # initialize empty list
  seurat_indiv <- list()
  
  for (row in 1:nrow(info_samples)) {
    
    # since these are .Rda objects (as opposed to .Rds), we need to use this
    # paradigm of get(load(.)) in order to save the contents of the .Rda into the
    # variable called `seurat`. (Used defensively in case the object in the .Rda
    # is not named seurat)
    seurat <- get(load(info_samples[row, ]$Path))
    
    seurat@project.name <- info_samples[row, ]$Sample
    seurat$Sample <- seurat@project.name
    print(seurat@project.name)
    
    # remove certain columns from the individual seurat objects,
    # for cleanliness, when columns from the individual space would not apply
    # in the joint space
    seurat@meta.data <- seurat@meta.data[, ! colnames(seurat@meta.data) %in% info_experiment$vars_drop]
    
    # put the covariate/group info in the metadata; it will be the same for all cells
    # since this is sample-level group info
    for (covariate in names(info_groups)) {
      metadata= matrix(rep(unlist(info_samples[row, covariate]),each = length(seurat@meta.data$orig.ident)))
      rownames(metadata)=rownames(seurat@meta.data)
      seurat <- AddMetaData(seurat,
                            metadata =  metadata,
                            col.name = covariate)
    }
    
    # subset to malignant cells only
    #if (!is.null(info_experiment$malignant_only) & info_experiment$malignant_only) {
    
    #    message("@ subsetting to malignant cells only...")
    
    #    seurat <- subset(seurat,
    #                     subset = Malignant_normal_consensus_Jessa2022 %in% c("Malignant", "Likely malignant"))
    
    #}
    
    # populate the list
    seurat_indiv[[row]] <- seurat
    
    # clean up
    rm(seurat)
    
  }
  
  # Preprocess data ------------------------------------------------------------
  
  message("@ preprocessing data...")
  
  # merge into a single Seurat object, normalize, scale, and run PCA
  seurat_joint <- merge(x = seurat_indiv[[1]],
                        y = seurat_indiv[2:length(seurat_indiv)],
                        merge.data = FALSE)
  
  # clean up
  rm(seurat_indiv)
  
  seurat_joint <- seurat_joint %>% 
    # all samples should have been normalized the same way, but just in case, 
    # re-run it here
    Seurat::NormalizeData(verbose = info_experiment$verbose) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = .@var.genes, npcs = info_experiment$n_pcs, verbose = info_experiment$verbose) %>% 
    RunTSNE(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, check_duplicates = FALSE) %>%
    RunUMAP(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100)
  
  # perform clustering
  seurat_joint <- seurat_joint %>% 
    FindNeighbors(seurat,
                  reduction = "pca",
                  dims      = 1:info_experiment$n_pcs,
                  verbose   = TRUE,
                  nn.eps    = 0.5) %>% 
    FindClusters(verbose = info_experiment$verbose,
                 n.start = 10,
                 random.seed = 100, 
                 resolution = 0.5)
  
  # Save -----------------------------------------------------------------------
  
  message("@ saving joined data...")
  
  seurat_joint_dr <- list(
    "pca"  = seurat_joint@reductions$pca@cell.embeddings[, c(1, 2)],
    "tsne" = seurat_joint@reductions$tsne@cell.embeddings[, c(1, 2)],
    "umap" = seurat_joint@reductions$umap@cell.embeddings[, c(1, 2)])
  saveRDS(seurat_joint_dr, file = "output/dimred.Rds")
  
  seurat_joint_meta <- seurat_joint@meta.data
  saveRDS(seurat_joint_meta, file = "output/metadata.Rds")
  
  save(seurat_joint, file = "output/seurat_joint.Rda")
  
  # Covariates ----
  
  # colour the low-dimensional embeddings by the covariates
  iwalk(info_groups, function(palette, covariate) {
    
    plot_fun <- purrr::partial(DimPlot,
                               object = seurat_joint,
                               group.by = covariate,
                               cols = palette)
    
    ggsave(plot = plot_fun(reduction = "pca"),  filename = glue("figures/PCA_{covariate}.png"),  width = 10, height = 8)
    ggsave(plot = plot_fun(reduction = "tsne"), filename = glue("figures/tSNE_{covariate}.png"), width = 10, height = 8)
    ggsave(plot = plot_fun(reduction = "umap"), filename = glue("figures/UMAP_{covariate}.png"), width = 10, height = 8)
    
  })
  
  # colour by clusters
  ggsave(plot = DimPlot(seurat_joint, reduction = "pca"),  filename = "figures/PCA_clusters.png",  width = 10, height = 8)
  ggsave(plot = DimPlot(seurat_joint, reduction = "tsne"), filename = "figures/tSNE_clusters.png", width = 10, height = 8)
  ggsave(plot = DimPlot(seurat_joint, reduction = "umap"), filename = "figures/UMAP_clusters.png", width = 10, height = 8)
  
}

# Integrate with Harmony -------------------------------------------------------

if (info_experiment$integrate) {
  
  message("@ integrating data with harmony...")
  
  # run the integration across the selected variables using Harmony; print the 
  # convergence plot for inspection
  png(filename = glue("figures/convergence.png"), width = 500, height = 400)
  seurat_joint_harmony <- seurat_joint %>% 
    RunHarmony(group.by.vars    = "orig.ident",
               
               reduction.use        = "pca",
               dims.use         = 1:info_experiment$n_pcs,
               plot_convergence = TRUE,
               verbose          = info_experiment$verbose)
  dev.off()
  
  # clean up
  rm(seurat_joint)
  
  # Downstream analysis ----------------------------------------------------------
  
  message("@ performing downstream analysis...")
  
  # run tSNE, UMAP, clustering
  seurat_joint_harmony <- seurat_joint_harmony %>% 
    RunTSNE(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony") %>%
    RunUMAP(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony")
  
  seurat_joint_harmony <- seurat_joint_harmony %>% 
    FindNeighbors(seurat,
                  reduction = "harmony",
                  dims      = 1:info_experiment$n_pcs,
                  verbose   = TRUE,
                  nn.eps    = 0.5) %>% 
    FindClusters(verbose = info_experiment$verbose,
                 n.start = 10,
                 random.seed = 100, 
                 resolution = 0.5)
  
  # Save ----
  
  message("@ saving integrated data...")
  
  seurat_joint_harmony_dr <- list(
    "pca"  = seurat_joint_harmony@reductions$pca@cell.embeddings[, c(1, 2)],
    "tsne" = seurat_joint_harmony@reductions$tsne@cell.embeddings[, c(1, 2)],
    "umap" = seurat_joint_harmony@reductions$umap@cell.embeddings[, c(1, 2)])
  saveRDS(seurat_joint_harmony_dr, file = "output/dimred.harmony.Rds")
  
  save(seurat_joint_harmony, file = "output/seurat_joint.harmony.Rda")
  
  # Covariates -------------------------------------------------------------------
  
  iwalk(info_groups, function(palette, covariate) {
    
    plot_fun <- purrr::partial(DimPlot,
                               object = seurat_joint_harmony,
                               group.by = covariate,
                               cols = palette)
    
    ggsave(plot = plot_fun(reduction = "pca"),  filename = glue("figures/PCA_{covariate}.harmony.png"),  width = 10, height = 8)
    ggsave(plot = plot_fun(reduction = "tsne"), filename = glue("figures/tSNE_{covariate}.harmony.png"), width = 10, height = 8)
    ggsave(plot = plot_fun(reduction = "umap"), filename = glue("figures/UMAP_{covariate}.harmony.png"), width = 10, height = 8)
    
  })
  
  # colour by clusters
  ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "pca"),  filename = "figures/PCA_clusters.harmony.png",  width = 10, height = 8)
  ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "tsne"), filename = "figures/tSNE_clusters.harmony.png", width = 10, height = 8)
  ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "umap"), filename = "figures/UMAP_clusters.harmony.png", width = 10, height = 8)
  
} else {
  
  message("@ Skipping harmony integration.")
  
}



# Session info -----------------------------------------------------------------
sessionInfo()
