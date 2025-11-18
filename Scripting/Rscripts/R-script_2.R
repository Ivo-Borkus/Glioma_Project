library(Seurat)
library(dplyr)
library(ggplot2)
library(here)
library(tidyr)
library(patchwork)
library(UCell)
library(qs2)
library(miloR)
library(SingleCellExperiment)
library(scater)
source(here("Brain-Met_V2/functions.R"))
exc_path <- here("Glioma_Project/01_data/Excel")
rmd_path <- here("Glioma_Project/02_knits/figs")
obj_path <- here("Glioma_Project/01_data")
fig_path <- here("Glioma_Project/04_figs")

.run_milo <- \(Processed_name, n_dims = 50, n_ks = 100, proportion = 0.1, sample_col = "sample_name"){
    print("Processing: ")
    print(Processed_name)
    print("Paramaters: ")
    print(paste0("Number of Dimensions: ", n_dims))
    print(paste0("Number of minimum cells per neighborhood: ", n_ks))
    print(paste0("Sampling proportion: ", proportion))
    print(paste0("Sample column: ", sample_col))

    seurat_obj_path <- ifelse(length(sub_dir) > 0, here(obj_path, sub_dir, paste0(Processed_name, ".qs2")), here(obj_path, paste0(Processed_name, ".qs2")))
    print(seurat_obj_path)

    obj <- qs_read(seurat_obj_path)
    print("object is read")
    .size(obj)

    UMAP_dim <- paste0("UMAP.HARMONY.", toupper(Processed_name))
    PCA_dim <- paste0("HARMONY_PCA.", toupper(Processed_name))

    sce.obj <- as.SingleCellExperiment(obj, assay = "RNA")
    milo.obj <- Milo(sce.obj)
    reducedDim(milo.obj, "UMAP") <- reducedDim(sce.obj, UMAP_dim)
    reducedDim(milo.obj, "PCA") <- reducedDim(sce.obj, PCA_dim)
    print("Milo object generated")

    milo.obj <- buildGraph(
        milo.obj,
        k = n_ks,
        d = n_dims,
        transposed = FALSE,
        get.distance = FALSE,
        reduced.dim = "PCA"
    )
    print("Graph is built")
    milo.obj <- makeNhoods(milo.obj, prop = proportion, k = n_ks, d = n_dims, refined = TRUE)
    plotNhoodSizeHist(milo.obj)
    milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), sample = sample_col)
    head(nhoodCounts(milo.obj))
    milo.obj <- calcNhoodDistance(milo.obj, d = n_dims)
    milo.obj <- buildNhoodGraph(milo.obj)

    return(milo.obj)
}



Processed_name <- "Myeloid_Harmony_1000"
sample_col <- "sample_name"
condition <- "tumour_full"
sub_dir <- NULL

d <- 50
k <- 200
prop <- 0.1
milo_obj <- .run_milo(Processed_name, n_dims = d, n_ks = k, proportion = prop, sample_col = sample_col)
milo_obj_path <- ifelse(length(sub_dir) > 0, here(obj_path, sub_dir, paste0(Processed_name, "_Milo.qs2")), here(obj_path, paste0(Processed_name, "_Milo.qs2")))
qs_save(milo_obj, milo_obj_path)
milo_obj <- qs_read(milo_obj)

traj_design <- data.frame(colData(milo_obj))[, c(sample_col, condition)]
traj_design <- distinct(traj_design)
dim(traj_design)

da_results <- testNhoods(milo_obj, design = reformulate(condition), design.df = traj_design)

da_results %>%
    arrange(-SpatialFDR) %>%
    head()

plotUMAP(milo_obj) + plotNhoodGraphDA(milo_obj, da_results, alpha = 0.05) +
    plot_layout(guides = "collect")
ggsave(paste0(fig_path, "Milo_overview_", processed_name, ".png"))

sessionInfo()
