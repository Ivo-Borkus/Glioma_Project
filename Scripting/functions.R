library(cowplot)
theme_set(theme_bw())


.extract_sample_names <- \(dir, list_of_files) {
    unique_samples <- vector()

    for (file_name in list_of_files) {
        file_name <- as.character(file_name) # I think this might be not necessary
        cat(file_name) # Just fun
        if (grepl("barcodes", file_name)) { # grab one of the files per sample
            sample_name <- gsub("*.barcodes.*", replacement = "", x = file_name) # extract everything before barcode
            unique_samples <- append(unique_samples, sample_name) # Make this into a list
        }
    }
    return(unique_samples)
}

.extracting_data <- \(dir){
    if (!dir.exists(file.path(dir))) {
        cat("directory does not exist")
    } else {
        list_of_files <- list.files(dir)
        unique_samples <- .extract_sample_names(dir, list_of_files)
        print(unique_samples)
        seurat_obj_list <- lapply(unique_samples, function(sample) {
            cat(sample)
            sample_files <- paste0(dir, list_of_files[grepl(sample, list_of_files)])
            barcodes_dir <- sample_files[grepl("barcodes", sample_files)]
            matrix_dir <- sample_files[grepl("matrix", sample_files)]
            features_dir <- sample_files[grepl("features", sample_files)]
            if (any(!c(file.exists(barcodes_dir), file.exists(matrix_dir), file.exists(features_dir)))
            ) {
                stop("warning directories don't exist")
            } else {
                matrix_read <- ReadMtx(mtx = matrix_dir, cells = barcodes_dir, features = features_dir)
                seurat_object <- CreateSeuratObject(counts = matrix_read)
                seurat_object[["percent_mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
                seurat_object$Sample <- sample
                return(seurat_object)
            }
        })
    }
    return(seurat_obj_list)
}

.size <- \(obj) {
    print(object.size(obj), units = "Gb")
}



.barplot <- \(obj, x_var, y_var, fill_var = NULL, ...){
    plot <- ggplot(
        data = obj,
        mapping = aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var))
    ) +
        geom_col(colour = NA) +
        theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size = 15), axis.title = element_text(size = 20), title = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(fill = fill_var) +
        list(...)

    return(plot)
}

.cell_count <- \(object, filling_var = NULL, sorting = TRUE){
    # This function expects a merged seurat object
    # It should contain meta data of cells with column name: Sample
    # the filling_var should be available in the metadata and be the same for all cells in a sample
    # Sorting orders the samples by cell count
    meta_data <- object@meta.data

    cell_data <- as.data.frame(table(meta_data$Sample))
    colnames(cell_data) <- c("Sample", "Cells")
    if (sorting) {
        print("Sorting the samples based on count")
        cell_data %>% arrange(Cells) -> cell_data
        cell_data$Sample <- factor(cell_data$Sample, levels = cell_data$Sample)
    }
    if (length(filling_var) == 0) {
        filling_var <- "Sample"
        pos <- "none"
    } else {
        treatment <- sapply(unique(meta_data$Sample), function(x) {
            unique(meta_data[[filling_var]][meta_data$Sample == x])
        })
        treatment <- sapply(treatment, function(x) x[1])

        treatment <- as.data.frame(treatment)
        treatment$Sample <- rownames(treatment)
        colnames(treatment) <- c(filling_var, "Sample")
        cell_data <- merge(cell_data, treatment, by = "Sample")
        pos <- "right"
    }
    # print(head(cell_data))
    .barplot(
        obj = cell_data, x_var = "Sample", y_var = "Cells", fill_var = filling_var,
        scale_y_continuous(expand = c(0, 0)),
        theme(legend.position = pos, axis.text.x = element_text(angle = 45, hjust = 1)),
        geom_hline(yintercept = 1000),
        coord_flip()
    )
}


.scrublet <- function(seurat_obj) {
    suppressPackageStartupMessages(library(rscrublet, character.only = TRUE)) # Need to install it again
    count_matrix <- t(as(seurat_obj@assays$RNA$counts, "TsparseMatrix"))
    print(count_matrix[1:2, 1:3])
    scrr <- scrub_doublets(E_obs = count_matrix, expected_doublet_rate = 0.06, min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30)
    scrr <- call_doublets(scrr)
    # plot_doublet_histogram(scrr)
    seurat_obj$doublet.score <- scrr$doublet_scores_obs
    seurat_obj$predicted.doublets <- scrr$predicted_doublets
    # print(FeaturePlot(seurat_obj, features = "doublet.score", cols = c("gray", "red")))
    return(seurat_obj)
}

.violin_plot <- \(metadata, x_var, y_var, fill_var = NULL, sorting = T, log_transf = T, int_low = NULL, int_high = NULL, ...){
    if (sorting) {
        metadata %>%
            group_by(!!sym(x_var)) %>%
            dplyr::summarize(average_count = median(!!sym(y_var))) %>%
            arrange(average_count) -> new_levels
        metadata[[x_var]] <- factor(metadata[[x_var]], levels = new_levels[[x_var]])
        print("Sorted")
    }
    if (log_transf) {
        print("Log transforming")
        metadata[[y_var]] <- log(metadata[[y_var]])
        if (length(int_low) > 0) {
            int_low <- log(int_low)
        }
        if (length(int_high) > 0) {
            int_high <- log(int_high)
        }
    }
    if (length(fill_var) == 0) {
        fill_var <- x_var
    }
    print(y_var)
    plot <- ggplot(
        data = metadata,
        mapping = aes(x = !!sym(x_var), y = !!sym(y_var), fill = !!sym(fill_var))
    ) +
        ggrastr::geom_violin_rast(trim = F, raster.dpi = 300, colour = NA) +
        geom_boxplot(width = 0.1, alpha = 0.5, colour = "black", show.legend = F) +
        theme(
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text = element_text(size = 15), axis.title = element_text(size = 20), title = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        {
            if (length(int_low) > 0) geom_hline(aes(yintercept = int_low, colour = "red"), show.legend = F)
        } +
        {
            if (length(int_high) > 0) geom_hline(aes(yintercept = int_high, colour = "red"), show.legend = F)
        } +
        labs(fill = fill_var) +
        list(...)
    return(plot)
}

.ncount_plot <- \(obj, fill_variable = NULL, min_umi = 500, max_umi = Inf){
    metadata <- obj@meta.data
    plot <- .violin_plot(metadata,
        x_var = "Sample",
        y_var = "nCount_RNA",
        fill_var = fill_variable,
        sorting = T,
        log_transf = T, int_low = min_umi, int_high = max_umi, coord_flip()
    )
    return(plot)
}

.nmt_plot <- \(obj, low_mt = 20, fill_variable = NULL){
    metadata <- obj@meta.data
    plot <- .violin_plot(metadata,
        x_var = "Sample",
        y_var = "percent_mt",
        fill_var = fill_variable,
        sorting = T,
        log_transf = F, int_low = low_mt, int_high = NULL, coord_flip()
    )
    return(plot)
}

.nfeat_plot <- \(obj, fill_variable = NULL, min_count = 500, max_count = 10000){
    metadata <- obj@meta.data
    plot <- .violin_plot(metadata,
        x_var = "Sample",
        y_var = "nFeature_RNA",
        fill_var = fill_variable,
        sorting = T,
        log_transf = T, int_low = min_count, int_high = max_count, coord_flip()
    )
    return(plot)
}

.scatter_plot <- \(metadata) {
    plot <- ggplot(
        data = metadata,
        mapping = aes(
            x = nCount_RNA, y = nFeature_RNA,
            color = percent_mt,
            size = percent_mt
        )
    ) +
        geom_point() +
        geom_hline(yintercept = 300, color = "red") +
        geom_hline(yintercept = 10000, color = "red") +
        geom_vline(xintercept = 500, color = "red") +
        # geom_vline(xintercept = maxUMI, color = "red") +
        scale_x_log10() +
        scale_y_log10() +
        labs(
            x = "UMI count", y = "Unique genes"
        )
    return(plot)
}

.qc_plot <- \(object){
    return(plot)
}





generating_metadata <- function(sample_id) {
    csv_list <- purrr::map(sample_id, function(sample) {
        # corresponding sample Id with folder name
        index <- match(sample, sample_id) # sample_fullid[index]
        ## Extracting the CSV from cellranger outputs given the Sample and directory names.

        csv <- read.csv(here("01_data", sample_fullid[index], sample, "outs/per_sample_outs", sample, "metrics_summary.csv"))
        csv$Metric.Value <- as.character(csv$Metric.Value) # To make sure they are all the same class
        csv <- csv %>%
            group_by(Metric.Name) %>%
            summarise(across(everything(), first), .groups = "drop") ## Get rid of redundant information
        csv$sample_name <- sample
        csv$sample_fullname <- sample_fullid[index]
        csv <- csv %>%
            select(sample_name, sample_fullname, Metric.Name, Metric.Value) %>%
            pivot_wider(names_from = Metric.Name, values_from = Metric.Value)
        ## Get the data in a way that contains a sample per row and useful info on each column name, instead of a long db
    })
    meta_data <- bind_rows(csv_list) # Combine all of the datasets (which were 1 row each) together.
    ## Getting out the percentage signs in the columns
    percentage_list <- sapply(meta_data, function(column) {
        any(grepl("%$", column))
    })
    meta_data[percentage_list] <- lapply(meta_data[percentage_list], function(x) as.numeric(sub("%", "", x)))
    conv_data <- type_convert(meta_data) # Convert all of the characters back into numeric etc.
    print(knitr::kable(conv_data[, 1:4], format = "markdown"))

    return(conv_data)
}

plot_settings <- theme_bw() + theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 20), title = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)
)

# y_limits <- scale_y_continuous(expand = c(0, 0), limits = c(0, max(conv_data$Cells) + 500))
.filtering_percent <- \(object, ncount = c(0, Inf), nfeat = c(0, Inf), nmt = 100){
    object$Sample %>% table() -> before
    seurat_obj_subs <- subset(object, nCount_RNA >= ncount[[1]] &
        nCount_RNA <= ncount[[2]] &
        nFeature_RNA >= nfeat[[1]] &
        nFeature_RNA <= nfeat[[2]] &
        percent_mt <= nmt)

    seurat_obj_subs$Sample %>% table() -> after
    percent_list <- ((before - after) / before) * 100
    return(percent_list)
}




#### Seurat processing


# umap_by_sample_primary <- DimPlot(
#     object = seurat_obj,
#     reduction = paste0("umap.", reduction_name),
#     group.by = "sample",
#     pt.size = 0.5,
#     label = F, cols = primary_colors_sample,
# ) & theme(plot.title = element_text(size = 10)) &
#     NoAxes() & labs(title = paste0("Data: ", reduction_name, " plot of sample/primary"))
# ggsave(umap_by_sample_primary, file = paste0(output_dir, reduction_name, "_Umap_sample_primary.png"), height = 20, width = 20)


# umap_by_annotation <- DimPlot(
#     object = seurat_obj,
#     reduction = paste0("umap.", reduction_name),
#     group.by = "general_annotation",
#     pt.size = 0.1,
#     label = T,
# ) & theme(plot.title = element_text(size = 10)) &
#     NoAxes() & labs(title = paste0("Data: ", reduction_name, " plot of annotation"))
# ggsave(umap_by_annotation, file = paste0(output_dir, reduction_name, "_Umap_annotation.png"), height = 20, width = 20)

# if (reduction_name == "merged_subset_immune") {
#     print("This is the immune subset")
#     umap_by_annotation <- DimPlot(
#         object = seurat_obj,
#         reduction = paste0("umap.", reduction_name),
#         group.by = "specific_annotation_immune",
#         pt.size = 0.1,
#         label = T,
#     ) & theme(plot.title = element_text(size = 10)) & labs(title = paste0("Data: ", reduction_name, " plot of specific immune annotation")) & NoAxes()
#     ggsave(umap_by_annotation, file = paste0(output_dir, reduction_name, "_Umap_specific_annotation.png"), height = 20, width = 20)
# } else {
#     print("this is not an immune subset")
#     print(reduction_name)
# }

.excel_sheet <- function(markers, output_dir, name) {
    library(writexl)
    print(paste0("Output will be put in: ", output_dir, paste0(name, ".xlsx")))
    if (file.exists(output_dir)) {
        markers %>%
            arrange(cluster, desc(avg_log2FC)) %>% # Arrange within each cluster
            group_by(cluster) %>%
            select(cluster, pct.1, pct.2, p_val, p_val_adj, avg_log2FC, gene) %>%
            group_split() %>% # Split into list by 'cluster'
            setNames(unique(markers$cluster)) %>% # Name list elements
            writexl::write_xlsx(here(output_dir, paste0(name, ".xlsx")))
    } else {
        stop("Directory does not exist")
    }
}

.regress_PCA <- \(seurat_obj, reduction_name, harmony, n_PCS, features){
    if (harmony) {
        reduction <- paste0("harmony_pca.", reduction_name)
    } else {
        reduction <- paste0("pca_", reduction_name)
    }
    matrix <- as.matrix(Embeddings(seurat_obj, reduction = reduction)[, 1:n_PCS])
    results <- list()
    for (feature in features) {
        temp_res <- sapply(colnames(matrix), function(pc) {
            model <- lm(matrix[, pc] ~ seurat_obj@meta.data[[feature]])
            summary(model)$r.squared
        })
        results[[feature]] <- temp_res
    }
    results <- as.data.frame(results)
    results$PCs <- paste0("PC_", 1:n_PCS)

    results$PCs <- factor(
        results$PCs,
        levels = paste0("PC_", 1:n_PCS)
    )
    df_corr <- as.data.frame(results) %>% tidyr::pivot_longer(cols = -c("PCs"), names_to = "QC", values_to = "cor")
    return(df_corr)
}

plotting_PCA_regress <- \(df_corr, limit_n = 0.5){
    plot_list <- list()
    features <- df_corr$QC %>% unique()
    for (feature in features) {
        plot_list[[feature]] <- ggplot(df_corr %>% filter(QC == feature), aes(x = PCs, y = cor, group = 1)) +
            geom_line() + # line plot
            geom_point() + # optional: show points on the line
            theme_classic() +
            ylim(c(0, limit_n)) +
            theme(axis.text.x = element_text(
                angle = 90, vjust = 1, hjust = 1, size = 6
            )) +
            labs(x = "PCs", y = "Correlation", title = "Correlation by PCs Faceted by QC")
    }
    return(plot_list)
}


.volcano_plotting <- function(marker_list, ident.1 = "", ident.2 = "") {
    marker_list$genes <- row.names(marker_list)
    marker_list$diffexpressed <- "NO"
    marker_list$diffexpressed[marker_list$avg_log2FC > 1.5 & marker_list$p_val < 0.05] <- "UP"
    marker_list$diffexpressed[marker_list$avg_log2FC < -1.5 & marker_list$p_val < 0.05] <- "DOWN"
    marker_list$diffexpressed[marker_list$avg_log2FC > 1.5 & marker_list$p_val < 0.05] <- "UP"
    marker_list$diffexpressed[marker_list$avg_log2FC < -1.5 & marker_list$p_val < 0.05] <- "DOWN"
    # marker_list$delabel <- NA
    # marker_list$delabel[marker_list$diffexpressed != "NO"] <- row.names(marker_list)[marker_list$diffexpressed != "NO"]

    marker_list$delabel <- NA
    marker_list$delabel[marker_list$diffexpressed != "NO"] <- row.names(marker_list)[marker_list$diffexpressed != "NO"]
    marker_list %>%
        arrange(desc(avg_log2FC)) %>%
        dplyr::select(avg_log2FC) %>%
        head(, n = 15) %>%
        row.names() -> labels_1
    marker_list %>%
        arrange(avg_log2FC) %>%
        dplyr::select(avg_log2FC) %>%
        head(, n = 15) %>%
        row.names() -> labels_2
    labels <- c(labels_1, labels_2)
    marker_list$delabel <- NA
    marker_list$delabel[marker_list$genes %in% labels] <- marker_list$genes[marker_list$genes %in% labels]
    min_above <- min(marker_list$p_val[marker_list$p_val > 0])
    marker_list$p_val <- ifelse(marker_list$p_val == 0, min_above, marker_list$p_val)
    volcano_plot <- ggplot(data = marker_list, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
        geom_point() +
        ggrepel::geom_text_repel(max.overlaps = Inf) +
        scale_color_manual(values = c("blue", "black", "red")) +
        geom_vline(xintercept = c(-1.5, 1.5), col = "red") +
        geom_hline(yintercept = -log10(0.05), col = "red") & labs(title = paste0("Comparing: ", ident.1, " versus the ", ident.2))
    return(volcano_plot)
}


.standard_processing <- \(obj, reduction_name, harmony = F, subset_col = NULL, cells_of_interest = NULL){
    if (length(subset_col) > 0) {
        Idents(seurat_obj) <- subset_col
        seurat_obj <- subset(x = seurat_obj, idents = cells_of_interest, invert = F)
    }
    seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj@meta.data[["Sample"]])
    if (harmony) {
        reduction_pca <- paste0("harmony_pca.", reduction_name)
        reduction_umap <- paste0("umap.harmony.", reduction_name)
    } else {
        reduction_pca <- paste0("pca_", reduction_name)
        reduction_umap <- paste0("umap_", reduction_name)
    }


    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
    # hvf_info <- HVFInfo(seurat_obj)
    # top_variable_genes <- hvf_info[order(hvf_info$variance.standardized, decreasing = TRUE), ]
    # top_2000_genes <- rownames(top_variable_genes)[1:2000]
    # VariableFeatures(seurat_obj) <- top_2000_genes
    # all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(object = seurat_obj, vars.to.regress = c("percent_mt"), verbose = FALSE)
    seurat_obj <- RunPCA(
        object = seurat_obj, nfeatures.print = 5, ndims.print = 1:2, reduction.name = paste0("pca_", reduction_name), npcs = 200
    )
    # ElbowPlot(seurat_obj, reduction = paste0("pca_", reduction_name), ndims = 100)
    # VizDimLoadings(seurat_obj, dims = 1:5, reduction = paste0("pca_", reduction_name), nfeatures = 15)

    if (harmony) {
        seurat_obj <- IntegrateLayers(
            object = seurat_obj, method = HarmonyIntegration,
            orig.reduction = paste0("pca_", reduction_name), new.reduction = paste0("harmony_pca.", reduction_name),
            verbose = TRUE
        )
    }

    seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_pca, dims = 1:25)
    seurat_obj <- RunUMAP(seurat_obj, reduction = reduction_pca, dims = 1:25, reduction.name = reduction_umap)

    seurat_obj <- JoinLayers(seurat_obj)
    return(seurat_obj)
}


.adjust_settings <- \(seurat_obj, minUMI, maxUMI, minfeat, maxfeat, maxmt){
    # Extract actual values
    actual_min_UMI <- min(seurat_obj@meta.data$nCount_RNA)
    actual_max_UMI <- max(seurat_obj@meta.data$nCount_RNA)
    actual_min_features <- min(seurat_obj@meta.data$nFeature_RNA)
    actual_max_features <- max(seurat_obj@meta.data$nFeature_RNA)
    actual_max_mito <- max(seurat_obj@meta.data$percent_mt)

    minUMI <- ifelse(is.na(minUMI), min(seurat_obj@meta.data$nCount_RNA), minUMI)
    maxUMI <- ifelse(is.na(maxUMI), max(seurat_obj@meta.data$nCount_RNA), maxUMI)
    minfeat <- ifelse(is.na(minfeat), min(seurat_obj@meta.data$nFeature_RNA), minfeat)
    maxfeat <- ifelse(is.na(maxfeat), max(seurat_obj@meta.data$nFeature_RNA), maxfeat)
    maxmt <- ifelse(is.na(maxmt), max(seurat_obj@meta.data$percent_mt), maxmt)

    if (minUMI < actual_min_UMI) {
        warning(paste("minUMI is lower than the actual minimum value. Adjusting minUMI to", actual_min_UMI))
        minUMI <- actual_min_UMI
    }


    if (maxUMI > actual_max_UMI) {
        warning(paste("maxUMI is greater than the actual maximum value. Adjusting maxUMI to", actual_max_UMI))
        maxUMI <- actual_max_UMI
    }

    if (minfeat < actual_min_features) {
        warning(paste("minfeat is lower than the actual minimum feature count. Adjusting minFeatures to", actual_min_features))
        minfeat <- actual_min_features
    }


    if (maxfeat > actual_max_features) {
        warning(paste("maxfeat is higher than the actual maximum feature count. Adjusting maxfeat to", actual_max_features))
        maxfeat <- actual_max_features
    }

    if (maxmt > actual_max_mito) {
        warning(paste("maxmt is higher than the actual maximum mitochondrial %. Adjusting maxmt to", actual_max_mito))
        maxmt <- actual_max_mito
    }
    return(c(minUMI, maxUMI, minfeat, maxfeat, maxmt))
}



## Seurat processing

.processing <- \(obj, reduction_name, Processed_name, subsetting = c(NULL, NULL),
    norm_meth_scale = c(TRUE, "LogNormalize", 10000), vf = c("vst", 2000, FALSE, FALSE), scaling = TRUE, var.to.regress = NULL, harmony = F, n_PCS = 25,
    clustering = c(TRUE, 4)){
    cat(paste0("Reduction name is: ", reduction_name, "\n"))
    cat(paste0("File name is: ", Processed_name, "\n"))

    cat(paste0("Normalising data: ", norm_meth_scale[1], "\n"))
    if (norm_meth_scale[1]) {
        cat(paste0("Normalisation method: ", norm_meth_scale[2], "\n"))
        cat(paste0("Normalisation scale factor: ", norm_meth_scale[3], "\n"))
    } else {
        cat("Not normalising data, as it is already processed before")
    }
    cat(paste0("Variable feature method: ", vf[1], "\n"))
    cat(paste0("Variable features: ", as.numeric(vf[2]), "\n"))
    cat(paste0("Perform rigid vf calling: ", vf[3], "\n"))
    cat(paste0("Quiet VDJ genes: ", vf[4], "\n"))
    cat(paste0("Perform scaling: ", scaling, "\n"))
    cat(paste0("Perform regression on variable(s): ", var.to.regress, "\n"))
    cat(paste0("Using PCs: ", as.numeric(n_PCS), "\n"))
    cat(paste0("Perform harmony integration: ", harmony, "\n"))
    cat(paste0("Perform clustering: ", clustering[1], "\n"))
    cat(paste0("Clustering algorithm: ", clustering[2], "\n"))

    if (length(subsetting[1]) > 0) {
        Idents(obj) <- subsetting[1]
        if (!length(subsetting[-1] > 0)) {
            warning("Please provide the cell of interest to subset")
        }
        cat(paste0("Subsetting seurat for ", subsetting[-1], " on column", subsetting[1], "\n"))

        obj <- subset(x = obj, idents = subsetting[-1], invert = F)
    }
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[["Sample"]])
    if (harmony) {
        reduction_pca <- paste0("harmony_pca.", reduction_name)
        reduction_umap <- paste0("umap.harmony.", reduction_name)
    } else {
        reduction_pca <- paste0("pca_", reduction_name)
        reduction_umap <- paste0("umap_", reduction_name)
    }
    if (norm_meth_scale[1]) {
        obj <- NormalizeData(
            object = obj,
            normalization.method = norm_meth_scale[2],
            scale.factor = as.numeric(norm_meth_scale[3])
        )
        cat(paste0("###################", "\n"))
        cat(paste0("Normalisation complete", "\n"))
    }
    obj <- FindVariableFeatures(obj, selection.method = vf[1], nfeatures = as.numeric(vf[2]), verbose = F)
    if (vf[3]) {
        # Calling variable features with more batch effect presence
        ## It calls the HVF based on most variance over all samples. Which means it might be skewed by batch effects.
        ## If F: uses seurat method, which takes the variables shared between samples and highest variance.
        hvf_info <- HVFInfo(obj)
        top_variable_genes <- hvf_info[order(hvf_info$variance.standardized, decreasing = TRUE), ]
        top_2000_genes <- rownames(top_variable_genes)[1:vf[2]]
        VariableFeatures(obj) <- top_2000_genes
        all.genes <- rownames(obj)
    }
    if (vf[4]) {
        ## Exclude VDJ genes
        # require(scRepertoire)
        # obj <- scRepertoire::quietVDJgenes(obj)
        VariableFeatures(obj) <- .remove_VDJ(VariableFeatures(obj))
    }
    cat(paste0("###################", "\n"))
    cat("FindVariableFeatures complete")
    if (scaling) {
        # Such as c("percent_mt")
        obj <- ScaleData(object = obj, features = VariableFeatures(obj), vars.to.regress = var.to.regress, verbose = FALSE)
        cat(paste0("###################", "\n"))
        cat("scaling complete \n")
    } else {
        warning("Data was not scaled")
    }

    obj <- RunPCA(
        object = obj, features = VariableFeatures(obj),
        nfeatures.print = 5, ndims.print = 1:2,
        reduction.name = paste0("pca_", reduction_name),
        npcs = 200
    )
    cat(paste0("###################", "\n"))
    cat("RunPCA complete \n")

    # ElbowPlot(obj, reduction = paste0("pca_", reduction_name), ndims = 100)
    # VizDimLoadings(obj, dims = 1:5, reduction = paste0("pca_", reduction_name), nfeatures = 15)

    if (harmony) {
        obj <- IntegrateLayers(
            object = obj, method = HarmonyIntegration,
            orig.reduction = paste0("pca_", reduction_name), new.reduction = paste0("harmony_pca.", reduction_name),
            verbose = TRUE
        )
        cat(paste0("###################", "\n"))
        cat("Running Harmony complete \n")
    } else {
        cat("Harmony integration was skipped \n")
    }

    obj <- FindNeighbors(obj, reduction = reduction_pca, dims = 1:as.numeric(n_PCS))
    if (clustering[1]) {
        res_values <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.5)
        obj <- FindClusters(obj, resolution = res_values, algorithm = clustering[2])
    }
    obj <- RunUMAP(obj, reduction = reduction_pca, dims = 1:as.numeric(n_PCS), reduction.name = reduction_umap)
    cat("Reduction created for PCA \n")
    cat(reduction_pca)

    cat("Reduction created for UMAP \n")
    cat(reduction_umap)
    obj <- JoinLayers(obj)

    seurat_obj_path <- ifelse(length(sub_dir) > 0, here(obj_path, sub_dir, paste0(Processed_name, ".qs2")), here(obj_path, paste0(Processed_name, ".qs2")))

    cat("saving seurat object into: \n")
    cat(seurat_obj_path)
    qs_save(obj, seurat_obj_path)
    .seurat_history()

    return(obj)
}


.seurat_dimplot <- \(seurat_obj, reduction = reduction_UMAP, group = "Sample", pt.size = 0.5, label = F, title = "dimplot", ...){
    plot <- DimPlot(
        object = seurat_obj,
        reduction = reduction,
        group.by = group,
        pt.size = pt.size,
        label = label,
        raster = TRUE
    ) &
        theme(plot.title = element_text(size = 10)) &
        NoAxes() & list(...)
    return(plot)
}


.feature_plot <- \(seurat_obj, reduction, features = umap_feature_vec, ncol = 3){
    features <- FeaturePlot(seurat_obj,
        features = features,
        reduction = reduction,
        ncol = ncol, order = TRUE, raster = TRUE
    ) & NoAxes() # & labs(title = paste0("Data: ", reduction_name, " plot of features "))
    return(features)
}


.Running_plots <- function(
    seurat_obj, reduction_name = "No_reduction_name_Given", output_dir = fig_path,
    harmony = FALSE, primary = c(FALSE, "primary_column"), umap_feature_vec = feature_vec, batch_cols) {
    plot_list <- list()
    metadata_list <- list()

    if (harmony) {
        reduction_UMAP <- paste0("umap.harmony.", reduction_name)
    } else {
        reduction_UMAP <- paste0("umap_", reduction_name)
    }
    print("running plots")

    resolution_columns <- grep("^RNA_snn_res\\.",
        colnames(seurat_obj@meta.data),
        value = TRUE
    )
    plot_list[["umap_resolution"]] <- .seurat_dimplot(seurat_obj, reduction = reduction_UMAP, group = resolution_columns, label = TRUE) & NoLegend()
    size_vec <- length(umap_feature_vec)
    a <- 0
    n <- 6
    for (i in seq(from = 1, to = size_vec, by = n)) {
        a <- a + 1
        if ((i + (n - 1)) > size_vec) {
            plot <- .feature_plot(seurat_obj, reduction = reduction_UMAP, features = umap_feature_vec[i:size_vec]) & NoLegend()
        } else {
            plot <- .feature_plot(seurat_obj, reduction = reduction_UMAP, features = umap_feature_vec[i:(i + (n - 1))]) & NoLegend()
        }
        plot_list[[paste0("Features: ", i, "-", (i + (n - 1)))]] <- plot
    }
    # plot_list[["Features"]] <- .feature_plot(seurat_obj, reduction = reduction_UMAP, features = umap_feature_vec) & NoLegend()


    for (batch_column in batch_cols) {
        metadata_list[[batch_column]] <- .seurat_dimplot(seurat_obj, reduction = reduction_UMAP, group = batch_column)
    }

    confounders <- c("nCount_RNA", "nFeature_RNA", "percent_mt")
    metadata_list[["Confounders_plot"]] <- .feature_plot(seurat_obj, reduction = reduction_UMAP, features = confounders, ncol = 1) & NoAxes()

    # ggsave(PCA_elbow, file = here(output_dir, paste0(reduction_name, "_PCA_elbow.png")))
    # ggsave(Genes_influence_PCA, file = here(output_dir, paste0(reduction_name, "_PCA_loadings.png")))
    # ggsave(umap_resolution_combined, file = here(output_dir, paste0(reduction_name, "_Harmony_", harmony, "_Umap_res.png")))
    # ggsave(features, file = here(output_dir, paste0(reduction_name, "_Harmony_", harmony, "_Umap_features.png")), width = 18, height = 20)
    # ggsave(umap_by_sample, file = here(output_dir, paste0(reduction_name, "_Harmony_", harmony, "_Umap_sample.png")), height = 20, width = 20)
    # ggsave(umap_by_batch, file = here(output_dir, paste0(reduction_name, "_Harmony_", harmony, "_Umap_batches.png")), height = 20, width = 20)
    # ggsave(confounders, file = here(output_dir, paste0(reduction_name, "_Harmony_", harmony, "_Umap_confounders.png")), width = 18, height = 20)
    plots <- list(plot_list, metadata_list)
    return(plots)
}

.tabsetter <- \(plot_list, overall_name = "Plots", path_tmp = rmd_path) {
    cat("## ", overall_name, " {.tabset} \n")
    cat("\n")

    size <- length(plot_list)
    for (i in 1:size) {
        header_name <- names(plot_list)[i]

        cat("### ", header_name, "\n")
        cat("\n")
        cat("\n")
        file_name <- paste0(path_tmp, paste0("/tmp_file_", overall_name, "_", i, ".png"))
        ggsave(filename = file_name, plot = plot_list[[i]], dpi = 300)
        cat(sprintf("![](%s)\n\n", file_name))
    }
}

.seurat_history <- \(){
    # try( exist_df <- read.table(here(obj_path, "seurat_history.csv"), sep = "/"))
    if (length(subsetting) == 0) {
        subsetting <- NA
        subsetted_cols <- NA
    } else {
        subsetted_cols <- do.call(paste, as.list(subsetting[-1]))
    }
    if (length(var.to.regress) == 0) {
        var.to.regress <- NA
    }

    df <- data.frame(
        input_seurat = object_name,
        reduction_name = reduction_name,
        sample_col = "Sample",
        Batch_col = "Batch",
        subset_col = subsetting[1],
        cells_of_interest = subsetted_cols,
        normalising = norm_meth_scale[1],
        normalise_method = norm_meth_scale[2],
        normalise_scale = norm_meth_scale[3],
        HVF_M = vf[1],
        HVF_N = vf[2],
        HVF_old = vf[3],
        HVF_Quiet = vf[4],
        scaling = scaling,
        Regressed_var = var.to.regress,
        harmony = harmony,
        Number_PC = as.numeric(n_PCS),
        Clustering = clustering[1],
        Clustering_alg = clustering[2],
        Resolution = NA,
        Annotation_level = NA,
        row.names = Processed_name
    )
    res <- try(read.table(here(obj_path, "seurat_history.csv"), sep = ","), silent = TRUE)
    if (class(res) == "try-error") {
        write.table(df, here(obj_path, "seurat_history.csv"), sep = ",")
    } else {
        exist_df <- read.table(here(obj_path, "seurat_history.csv"), sep = ",")
        exist_df <- readr::type_convert(exist_df) # Convert all of the characters back into numeric etc.
        df <- readr::type_convert(df)
        combined_data <- bind_rows(exist_df, df)
        write.table(combined_data, here(obj_path, "seurat_history.csv"), sep = ",")
    }
    cat("Info put into: seurat_history.csv")
}



.lisi_running <- function(seurat_obj, reduction_name, n_PCS, batch_col = "Batch") {
    # Embeddings(seurat_obj, reduction = paste0("harmony_pca.", reduction_name))[, 1:n_PCS] %>% compute_lisi(, meta_data = seurat_obj@meta.data, label_colnames = c(batch_col)) -> lis
    # Embeddings(seurat_obj, reduction = paste0("pca_", reduction_name))[, 1:n_PCS] %>% compute_lisi(, meta_data = seurat_obj@meta.data, label_colnames = c(batch_col)) -> lis_no_int
    PC_harmony <- Embeddings(seurat_obj, reduction = paste0("harmony_pca.", reduction_name))[, 1:n_PCS]
    lis <- compute_lisi(PC_harmony, meta_data = seurat_obj@meta.data, label_colnames = c(batch_col))
    PCs_unintegrated <- Embeddings(seurat_obj, reduction = paste0("pca_", reduction_name))[, 1:n_PCS]
    lis_no_int <- compute_lisi(PCs_unintegrated, meta_data = seurat_obj@meta.data, label_colnames = c(batch_col))

    lis_no_int$integrated <- "no"
    lis$integrated <- "yes"
    combined <- rbind(lis_no_int, lis)
    combined$Batch <- combined[[batch_col]]
    plot <- ggplot(combined, mapping = aes(x = integrated, y = Batch)) +
        geom_boxplot() +
        theme_minimal() & labs(title = paste0("Lisi score pre and post integration with a significand t-test statistic: ", round(test_result$statistic, 2)))
    plot
    ggsave(plot, filename = here(fig_path, paste0("Lisi_samples_", reduction_name, ".png")))
    return(combined)
}



.run_milo <- \(Processed_name, n_dims = 50, n_ks = 100, proportion = 0.1, sample_col = "sample_name", reduction_name = Processed_name){
    print("Processing: ")
    print(Processed_name)
    print("reduction name: ")
    print(reduction_name)
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

    UMAP_dim <- paste0("UMAP.HARMONY.", toupper(reduction_name))
    PCA_dim <- paste0("HARMONY_PCA.", toupper(reduction_name))

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
    print("run-milo is finished")
    return(milo.obj)
}


.run_DAtest <- \(milo.obj, condition, contrast.1){
    traj_design <- data.frame(colData(milo.obj))[, c(sample_col, condition)]
    traj_design <- distinct(traj_design)
    dim(traj_design)
    rownames(traj_design) <- traj_design[, "sample_name"]
    traj_design[[condition]] <- factor(traj_design[[condition]])
    da_results <- testNhoods(milo.obj,
        design = reformulate(paste0(" 0 + ", condition)),
        design.df = traj_design, model.contrasts = contrast.1,
        fdr.weighting = "graph-overlap", norm.method = "logMS"
    )

    model <- model.matrix(reformulate(paste0(" 0 + ", condition)), data = traj_design)
    # da_results %>%
    #   arrange(SpatialFDR) %>%
    #   head(n = 10)

    output_figure <- DimPlot(seurat_obj, group.by = "Ann_lvl_1H", cols = cell_colors_2, reduction = reduction_umap) + NoAxes() + plotNhoodGraphDA(milo.obj, da_results, alpha = 0.10) +
        plot_layout(guides = "collect") + plot_annotation(title = "Milo DA results (spatialFDR: 10%)", subtitle = contrast.1)

    da_results <- annotateNhoods(milo.obj, da_results, coldata_col = "Ann_lvl_2")
    ggplot(da_results, aes(Ann_lvl_2_fraction)) +
        geom_histogram(bins = 50)

    da_results$celltype <- ifelse(da_results$Ann_lvl_2_fraction < 0.7, "Mixed", da_results$Ann_lvl_2)
    output_figure_2 <- plotDAbeeswarm(da_results, group.by = "celltype") + plot_annotation(title = "Milo DA results (spatialFDR: 10%)", subtitle = contrast.1)

    ggsave(filename = here(fig_path, "milo", paste0("miloR_umap", contrast.1, "_", k, ".png")), plot = output_figure)
    ggsave(filename = here(fig_path, "milo", paste0("miloR_beeswarm", contrast.1, "_", k, ".png")), plot = output_figure_2)
    return(da_results)
}
.remove_VDJ <- \(genes){
    genes <- genes[!grepl("^TR[ABDG][VDJ][^D]", genes, ignore.case = TRUE)]
    genes <- genes[grepl("^IG[HLK][VDJCAGM]", genes, ignore.case = TRUE) & toupper(genes) == "JCHAIN"  |  !genes %in% c(pseudo_genes_BCR)]
    return(genes)
}