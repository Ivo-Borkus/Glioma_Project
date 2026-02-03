exc_path <- here("Glioma_Project/01_data/Excel")
rmd_path <- here("Glioma_Project/02_knits/figs")
obj_path <- here("Glioma_Project/01_data")
fig_path <- here("Glioma_Project/04_figs")

batch_cols <- c(
    "tumour_full", "tumour_abbreviated", "grade", "Cell_Type", "Cell_Subtype", "sample_name", "IDH_status",
    "FACS_sorting", "tier_correct", "cell_het_paper_id", "storage_condition", "prior_treatment", "gender", "Tumour_Margin",
    "Primary_or_Recurrent", "Previous_Surgery", "Radiotherapy", "Chemotherapy", "MGMT_Methylation", "Quality"
)
continous_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "age_at_surgery_correct", "nuclear_fraction")

sample_col <- "sample_name"
feature_vec_overall <- c(
    "PTPRC", "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "TOP2A", "MKI67", # T and prolif
    "CD1C", "HLA-DQB1", # DC
    "S100AB", "TREM2", "LGALS3", # Monocyte macrophage
    "TMEM119", # Microglial
    "MOG" # Oligos
)
feature_vec_T <- c("PTPRC", "CD3E", "CD3D", "CD8B", "CD4", "LAG3", "CX3CR1", "CCR7", "SELL", "NCAM1", "NKG7", "KLRB1", "FOXP3", "TOP2A", "MKI67", "PDCD1", "FCGR3A")


# CD45 negative cluster main markers
markers <- list(
    "Tumor SEMA3D" = c("SEMA3D", "SOX4", "TNC"),
    "Tumor CDH18" = c("CDH18", "SLIT2", "EPHA5", "KRT19"),
    "Tumor Cycling" = c("TOP2A", "MKI67"),
    "Tumor SGCZ" = c("PTPRD", "SGCZ", "DGKB"),
    "Tumor SEMA6D" = c("SEMA6D", "COL12A1", "MYOC"),
    "Tumor TENM2" = c("TENM2", "COL25A1", "KCND2"),
    "Fibroblast" = c("COL4A1", "ACTA2", "TAGLN"),
    "Tumor FBN2" = c("FBN2", "COL19A1", "FN1", "MEOX2"),
    "Endothelial cells" = c("PECAM1", "ANGPT2")
)

# Basic immune markers to identify cell populations
immune_markers <- list(
    "CD4 Naive/CM" = c("CD4", "ANXA1", "PASK", "SELL", "LEF1", "NOSIP", "CCR7", "TCF7", "ACTN1", "FOXP1", "KLF2", "ITGA6", "CD8A-", "CD8B-", "GZMK-"),
    "CD4 Effector/Mem" = c("CD4", "ZNF683", "KLRB1", "PRDM1", "CX3CR1", "EOMES", "KLRG1", "TNFSF13B", "GZMK", "CCL5", "CCL4", "NKG7", "CD69", "ITGAE", "CD8A-", "CD8B-"),
    "T helper" = c("CD4", "CXCR3", "GATA3", "RORC", "RORA", "IL17F", "IL17A", "CCR6", "CXCR6", "IFNG", "IL4", "IL6ST", "CXCR5", "CXCL13", "PDCD1", "CD8A-", "CD8B-"),
    "CD4 IFN response" = c("CD4", "IFI16", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "CD8A-", "CD8B-"),
    "CD4 Proliferative" = c("CD4", "MKI67", "TOP2A", "STMN1", "UBE2C", "PCLAF", "CENPF", "CDK1", "CD8A-", "CD8B-"),
    "T reg" = c("IL32", "CCR7", "LEF1", "TCF7", "FOXP3", "CTLA4", "IL2RA", "ICOS", "TIGIT", "TOX2", "IKZF2", "GATA3", "CD28", "CD8A-", "CD8B-"),
    "Gamma Delta" = c("TRGC1", "TRGC2", "TRDC", "CD8A-", "CD8B-", "CD4-"),
    "CD8 Naive/CM" = c("CD4-", "ANXA1", "PASK", "SELL", "LEF1", "NOSIP", "CCR7", "TCF7", "ACTN1", "FOXP1", "KLF2", "ITGA6", "CD8A", "CD8B", "GZMK-"),
    "CD8 Mem" = c("CD8A", "CD8B", "ZNF683", "KLRB1", "PRDM1", "CX3CR1", "EOMES", "KLRG1", "TNFSF13B", "CD4-"),
    "CD8 Cytotoxic" = c("CD8A", "CD8B", "GZMK", "GZMH", "CCL5", "CCL4", "NKG7", "CD69", "PRF1", "ITGAE", "CD4-", "CST7", "GZMA", "CCL4L2", "KLRG1", "CTSW", "GZMH", "GZMM", "KLRK1", "HLA-C", "PRF1", "XCL2", "XCL1"),
    "CD8 IFN response" = c("CD8A", "CD8B", "IFI16", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "ISG15", "CD4-"),
    "CD8 Exhausted" = c("CD8A", "CD8B", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX", "TOX2", "LAYN", "CTLA4", "CD4-"),
    "CD8 Proliferative" = c("CD8A", "CD8B", "MKI67", "TOP2A", "STMN1", "UBE2C", "PCLAF", "CENPF", "CDK1", "CD4-"),
    "NK" = c("NCAM1", "FCGR3A", "CX3CR1", "GNLY", "KLRC2", "KLRD1", "KLRC3", "KLRK1", "KLRC1", "GNLY", "NKG7"),
    "Naive B cell" = c("MS4A1", "IGHD", "IGHM", "CCR7", "SELL", "TCL1A", "CD79A", "VPREB3", "FCRL1", "NIBAN3", "CD79B", "HVCN1", "CD72", "FCER2", "CD83", "CD19", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
    "Memory B cell" = c("CD79A", "MS4A1", "CD27", "TNFRSF13B", "ITGAX", "PRDM1", "CD24", "BANK1", "CD74", "HLA-DRA", "IGHA1", "BLK", "SPIB", "P2RX5", "IGHA2", "CD37", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
    "Plasma cells" = c("MZB1", "SDC1", "IGHG1", "JCHAIN", "IGHA1", "IGHG3", "IGLC3", "IGLC1", "IGHGP", "DERL3", "IGHG4", "XBP1", "IRF4", "CD3E-", "CD3G-", "CD3D-", "CD4-", "CD8A-", "CD8B-"),
    "Monocytes" = c("CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1"),
    "Non_classical monocytes" = c("FCGR3A", "CX3CR1", "HLA-DRB1", "HLA-DRA")
)

# general signatures describing functions of immune infiltrates
signatures <- list(
    "Cytotoxic T-cell" = c(
        "CD8A", "CD8B", "PRF1", "CCL5", "GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "IFNG",
        "NKG7", "KLRD1", "KLRK1", "IL2", "IL12RB2", "EOMES", "TBX21"
    ),
    "Exhausted T-cell" = c(
        "PDCD1", "CTLA4", "CXCL13", "TOX", "LAG3", "PRDM1", "TIGIT", "BATF", "HNF1A", "LAYN",
        "CD38", "ENTPD1", "HAVCR2", "FASLG", "PDCD1LG2", "CD160", "EOMES", "IRF4"
    ),
    "Tregs" = c(
        "CD4", "FOXP3", "CTLA4", "IL2RA", "RORA", "ENTPD1", "TNFRSF9", "TNFRSF18", "LGALS3", "ITGAE",
        "S1PR1", "CCR9", "ITGA4", "STAT5", "FOXO1", "IKZF2", "LRRC32", "IL10", "TGFB1", "TNFRSF25"
    ),
    "M2 M" = c(
        "ARG1", "CCL18", "CD163", "CD209", "CD274", "CHIT1", "CSF1", "CSF1R", "IL10", "ITGA4",
        "LGALS9", "MARCO", "MRC1", "RNASE1", "SELENOP", "SPP1", "TGFB1", "TGFB2", "TREM2",
        "MRC2", "VEGFA", "MS4A4A", "CD200R1", "STAB1", "IL4R", "IL13RA1"
    ),
    "M1 M" = c(
        "C1QA", "C1QC", "CCL2", "IL1B", "CCL4", "CCL7", "CCL8", "NFKBIA", "CD40", "CXCL2",
        "CXCL3", "CXCL9", "CXCL10", "CXCL11", "IDO1", "NFKB1", "TNF", "CXCL8", "G0S2",
        "IL6", "INHBA", "S100A8", "S100A9", "IL12A", "IL12B", "IRF5", "STAT1", "NOS2", "FCGR1A"
    )
)

cell_colors <- c(
    "Mg_classic"     = "#5E4FA2",
    "Mg_Inflam"      = "lightpink",
    "Mg_Anti-Inflam" = "#4B99C6",
    "Mg_PLCG2"       = "darkblue",
    "Mg_active"      = "#2A6F97",
    "BAMs"           = "#D73027",
    "TAMs"           = "#FC8D59",
    "Ma-Inflam"      = "#E6550D",
    "Ma-IFN"         = "#FDAE61",
    "Ma-Mt+"         = "#FEE08B",
    "Mono-Classic"   = "#A53E2B",
    "Myeloid_Prolif" = "cyan",
    "DCs_CD1C+"      = "gold",
    "pDCs"           = "yellow",
    "Mast-cells"     = "#008B8B",
    "CD8+_TNFSF9"    = "darkgreen",
    "CD8+_Cytotoxic" = "#2F9A7D",
    "CD8+ Eff-Exh"   = "skyblue",
    "CD4+ cm"        = "#2B60DE",
    "CD4+ Tfh"       = "#74ADD1",
    "Tregs"          = "#2ECC71",
    "T-IFN"          = "#313695",
    "T_Prolif"       = "#FFBF00",
    "MAIT-cells"     = "pink",
    "NK"             = "brown",
    "Hb+"            = "#BDBDBD",
    "Bad-quality"    = "#000000"
)
cell_colors_2 <- c(
    "Microglia" = "darkblue",
    "Macrophage" = "red",
    "Monocyte" = "#A53E2B",
    "Myeloid_Prolif" = "cyan",
    "DCs" = "gold",
    "Mast-cells" = "#008B8B",
    "T-CD8+" = "darkgreen",
    "T-CD4+" = "#2ECC71",
    "T-other" = "pink",
    "T-Prolif" = "#FFBF00",
    "NK" = "brown",
    "Hb+" = "#BDBDBD",
    "Bad-quality" = "#000000"
)

Mgenes <- c(
    "TREM2", "APOE",
    "P2RY12", "SLC12A5",
    "IL1B", "CCL3L1",
    "SELENOP",
    "PLCG2",
    "TNF",
    "LYVE1", "F13A1",
    "CD36", "CCL18", "LGALS3",
    "SLC2A1",
    "IFIT1",
    "MT1H",
    "S100A12",
    "TOP2A", "MKI67",
    "CD1C", "CLEC10A",
    "IL3RA", "LILRA4",
    "HDC", "IL1RL1",
    "HBA1"
)
