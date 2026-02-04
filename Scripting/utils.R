exc_path <- here("Glioma_Project/01_data/Excel")
rmd_path <- here("Glioma_Project/02_knits/figs")
obj_path <- here("Glioma_Project/01_data")
fig_path <- here("Glioma_Project/04_figs")
blank <- element_blank()
x_axis_theme <- theme(axis.text.x = element_text(size = 11, angle = 30, hjust = 1, color = "black"))
y_axis_theme <- theme(axis.text.y = element_text(size = 11, color = "black"))

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


## Extracted from scRepertoire package!

pseudo_genes_BCR <- unique(c(
        "IGHJ1P", "IGHJ2P", "IGHJ3P", "IGLC4", "IGLC5", "IGHEP1", "IGHEP2",
        "IGHV1-12","IGHV1-14", "IGHV1-17", "IGHV1-67", "IGHV1-68",
        "IGHV2-10", "IGHV3-6", "IGHV3-19", "IGHV3-22", "IGHV3-25",
        "IGHV3-29", "IGHV3-32", "IGHV3-36", "IGHV3-37", "IGHV3-41",
        "IGHV3-42", "IGHV3-47", "IGHV3-50", "IGHV3-52", "IGHV3-54",
        "IGHV3-57", "IGHV3-60", "IGHV3-62", "IGHV3-63", "IGHV3-65",
        "IGHV3-71", "IGHV3-75", "IGHV3-76", "IGHV3-79", "IGHV4-55",
        "IGHV4-80", "IGHV5-78", "IGHV7-27", "IGHV7-40", "IGHV7-56",
        "IGHVIII-44", "IGHVIII-82", "IGKV1-22", "IGKV1-32", "IGKV1-35",
        "IGKV1D-22", "IGKV1D-27", "IGKV1D-32", "IGKV1D-35", "IGKVOR-2",
        "IGKVOR-3", "IGKVOR-4", "IGKV2-4", "IGKV2-10", "IGKV2-14", "IGKV2-18",
        "IGKV2-19", "IGKV2-23", "IGKV2-26", "IGKV2-36", "IGKV2-38",
        "IGKV2D-10", "IGKV2D-14", "IGKV2D-18", "IGKV2D-19", "IGKV2D-23",
        "IGKV2D-36", "IGKV2D-38", "IGKV3-25", "IGKV3-31", "IGKV3-34",
        "IGKV7-3", "IGLCOR22-1", "IGLCOR22-2", "IGLJCOR18", "IGLV1-41",
        "IGLV1-62", "IGLV2-5", "IGLV2-28", "IGLV2-34", "IGLV3-2",
        "IGLV3-4", "IGLV3-6", "IGLV3-7", "IGLV3-13", "IGLV3-15",
        "IGLV3-17", "IGLV3-24", "IGLV3-26", "IGLV3-29", "IGLV3-30",
        "IGLV3-31", "IGLV7-35", "IGLV10-67", "IGLVI-20", "IGLVI-38",
        "IGLVI-42", "IGLVI-56", "IGLVI-63", "IGLVI-68", "IGLVI-70",
        "IGLVIV-53", "IGLVIV-59", "IGLVIV-64", "IGLVIV-65", "IGLVV-58",
        "IGLVV-66", "IGHV1OR15-2", "IGHV1OR15-3", "IGHV1OR15-4", "IGHV1OR15-6",
        "IGHV1OR16-1", "IGHV1OR16-2", "IGHV1OR16-3", "IGHV1OR16-4", "IGHV3-30-2",
        "IGHV3-33-2", "IGHV3-69-1", "IGHV3OR15-7", "IGHV3OR16-6", "IGHV3OR16-7",
        "IGHV3OR16-11", "IGHV3OR16-14", "IGHV3OR16-15", "IGHV3OR16-16", "IGHV7-34-1",
        "IGHVII-1-1", "IGHVII-15-1", "IGHVII-20-1", "IGHVII-22-1", "IGHVII-26-2",
        "IGHVII-28-1", "IGHVII-30-1", "IGHVII-30-21", "IGHVII-31-1", "IGHVII-33-1",
        "IGHVII-40-1", "IGHVII-43-1", "IGHVII-44-2", "IGHVII-46-1", "IGHVII-49-1",
        "IGHVII-51-2", "IGHVII-53-1", "IGHVII-60-1", "IGHVII-62-1", "IGHVII-65-1",
        "IGHVII-67-1", "IGHVII-74-1", "IGHVII-78-1", "IGHVIII-2-1", "IGHVIII-5-1",
        "IGHVIII-5-2", "IGHVIII-11-1", "IGHVIII-13-1", "IGHVIII-16-1", "IGHVIII-22-2",
        "IGHVIII-25-1", "IGHVIII-26-1", "IGHVIII-38-1", "IGHVIII-47-1", "IGHVIII-67-2",
        "IGHVIII-67-3", "IGHVIII-67-4", "IGHVIII-76-1", "IGHVIV-44-1", "IGKV1OR1-1",
        "IGKV1OR2-1", "IGKV1OR2-2", "IGKV1OR2-3", "IGKV1OR2-6", "IGKV1OR2-9",
        "IGKV1OR2-11", "IGKV1OR2-118", "IGKV1OR9-1", "IGKV1OR9-2", "IGKV1OR10-1",
        "IGKV1OR15-118", "IGKV1OR22-1", "IGKV1OR22-5", "IGKV1ORY-1", "IGKV2OR2-1",
        "IGKV2OR2-2", "IGKV2OR2-4", "IGKV2OR2-7", "IGKV2OR2-7D", "IGKV2OR2-8",
        "IGKV2OR2-10", "IGKV2OR22-3", "IGKV2OR22-4", "IGKV3OR2-5", "IGKV3OR22-2",
        "IGKV8OR8-1", "IGLVIV-66-1", "IGLVIVOR22-1", "IGLVIVOR22-2", "IGLVVI-22-1",
        "IGLVVI-25-1", "IGLVVII-41-1"
    ))