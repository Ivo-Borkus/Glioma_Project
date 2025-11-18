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
    # print("reduction name: ")
    # print("Inhouse_Harmony")
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

    UMAP_dim <- paste0("UMAP.HARMONY.", toupper("Inhouse_Harmony"))
    PCA_dim <- paste0("HARMONY_PCA.", toupper("Inhouse_Harmony"))

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
    milo.obj <- makeNhoods(milo.obj, prop = proportion, k = n_ks, d = n_dims, refined = TRUE, refinement_scheme = "graph") # https://github.com/MarioniLab/miloR/issues/249
    plotNhoodSizeHist(milo.obj) + geom_vline(aes(xintercept = 100))
    milo.obj <- countCells(milo.obj, meta.data = data.frame(colData(milo.obj)), sample = sample_col)
    head(nhoodCounts(milo.obj))
    milo.obj <- buildNhoodGraph(milo.obj)
    return(milo.obj)
}

seurat_obj$IDH_status

Processed_name <- "Harmony"
sample_col <- "sample_name"
condition <- "IDH_status"
sub_dir <- NULL

d <- 50
k <- 100
prop <- 0.1
for (k in range(k_values)) {
    milo.obj <- .run_milo(Processed_name, n_dims = d, n_ks = k, proportion = prop, sample_col = sample_col)
}
milo_obj_path <- ifelse(length(sub_dir) > 0, here(obj_path, sub_dir, paste0(Processed_name, "_", k, "_Milo.qs2")), here(obj_path, paste0(Processed_name, "_Milo.qs2")))
qs_save(milo.obj, milo_obj_path)

traj_design <- data.frame(colData(milo.obj))[, c(sample_col, condition)]
traj_design <- distinct(traj_design)
dim(traj_design)
rownames(traj_design) <- traj_design[, "sample_name"]
traj_design[[condition]] <- factor(traj_design[[condition]])


contrast.1 <- c("tumour_fullastrocytoma - tumour_fulloligodendroglioma") # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>
reformulate(paste0(" 0 + ", condition))
# we need to use the ~ 0 + Variable expression here so that we have all of the levels of our variable as separate columns in our model matrix
da_results <- testNhoods(milo.obj,
    design = reformulate(paste0(" 0 + ", condition)),
    design.df = traj_design, model.contrasts = contrast.1,
    fdr.weighting = "graph-overlap", norm.method = "logMS"
)


contrast.1 <- c("astrocytoma - tumour_fulloligodendroglioma") # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>

da_results <- testNhoods(milo.obj,
    design = reformulate(paste0(" 0 + ", condition)),
    design.df = traj_design, model.contrasts = contrast.1,
    fdr.weighting = "graph-overlap", norm.method = "logMS"
)



# da_results <- testNhoods(milo_obj, design = reformulate(condition), design.df = traj_design, fdr.weighting = "graph-overlap")
milo_obj <- qs_read(milo_obj)

traj_design <- data.frame(colData(milo_obj))[, c(sample_col, condition)]
traj_design <- distinct(traj_design)
dim(traj_design)

da_results <- testNhoods(milo_obj, design = reformulate(condition), design.df = traj_design)
da_results %>% colnames()
da_results %>%
    arrange(SpatialFDR) %>%
    head(n = 10)


milo.obj <- buildNhoodGraph(milo.obj)

reduction_pca <- paste0("harmony_pca.", "Inhouse_Harmony")
reduction_umap <- paste0("umap.harmony.", "Inhouse_Harmony")
dittoSeq::dittoDimPlot(seurat_obj, var = "Ann_lvl_1", reduction.use = reduction_umap, color.panel = cell_colors)
# dittoSeq::dittoDimPlot(seurat_obj_T, var = "Ann_lvl_2T", color.panel = cell_colors, reduction.use = "umap.harmony.T_cells_harmony_1000", do.label = T) & NoAxes()
ggsave(here(fig_path, paste0("Ann_lvl_2T.png")), dpi = 400)

output_figure <- DimPlot(seurat_obj, group.by = "Ann_lvl_1H", cols = cell_colors_2, reduction = reduction_umap) + NoAxes() + plotNhoodGraphDA(milo.obj, da_results, alpha = 0.10) +
    plot_layout(guides = "collect") + plot_annotation(title = "Milo DA results (spatialFDR: 10%)", subtitle = contrast.1)

ggsave(paste0(fig_path, "Milo_", processed_name, ".png"))

ggsave(paste0(fig_path, "Milo_", processed_name, ".png"))

ggplot(da_results, aes(PValue)) +
    geom_histogram(bins = 50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

da_results <- annotateNhoods(milo.obj, da_results, coldata_col = "Ann_lvl_2")
ggplot(da_results, aes(Ann_lvl_2_fraction)) +
    geom_histogram(bins = 50)

da_results$celltype <- ifelse(da_results$Ann_lvl_2_fraction < 0.7, "Mixed", da_results$Ann_lvl_2)
plotDAbeeswarm(da_results, group.by = "celltype")
ggsave(paste0(fig_path, "Milo_", processed_name, ".png"))

sessionInfo()
da_results %>% colnames()
