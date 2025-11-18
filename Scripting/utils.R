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
