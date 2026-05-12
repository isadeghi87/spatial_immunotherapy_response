#!/usr/bin/env Rscript

.spata2_dir <- getOption("spata2_script_dir", default = getwd())

PROJECT_ROOT <- normalizePath(file.path(.spata2_dir, "..", "..", ".."), mustWork = FALSE)
SPATIAL_OUT_ROOT <- file.path(PROJECT_ROOT, "out", "02b_spatial_biomarkers")

FIG_DIR <- file.path(SPATIAL_OUT_ROOT, "figures")
TAB_DIR <- file.path(SPATIAL_OUT_ROOT, "tables")
QC_DIR <- file.path(SPATIAL_OUT_ROOT, "qc")
QC_FIG_DIR <- file.path(QC_DIR, "plots")
QC_TAB_DIR <- file.path(QC_DIR, "tables")
OBJ_DIR <- file.path(SPATIAL_OUT_ROOT, "spata2_objects")
NOTE_DIR <- file.path(SPATIAL_OUT_ROOT, "notes")

CNV_ROOT <- file.path(SPATIAL_OUT_ROOT, "cnv")
CNV_FIG_DIR <- file.path(CNV_ROOT, "figures")
CNV_TAB_DIR <- file.path(CNV_ROOT, "tables")
CNV_RDS_DIR <- file.path(CNV_ROOT, "rds")
CNV_QC_DIR <- file.path(CNV_ROOT, "qc")

GSE238264_DIR <- file.path(PROJECT_ROOT, "data", "spatial_data", "GSE238264_HCC_IO_Visium")
GSE281759_DIR <- file.path(PROJECT_ROOT, "data", "spatial_data", "GSE281759_HCC_MultiRegion")
GTF_PATH <- file.path(PROJECT_ROOT, "data", "references", "gencode.v44.basic.annotation.gtf.gz")

SPATA_SEED <- 1337L

GSE238264_SAMPLES <- c(
  HCC1R = "Responder",
  HCC2R = "Responder",
  HCC3R = "Responder",
  HCC4R = "Responder",
  HCC5NR = "Non_Responder",
  HCC6NR = "Non_Responder",
  HCC7NR = "Non_Responder"
)

infer_gse281759_group <- function(sample_id) {
  short_id <- sub("^GSM[0-9]+_", "", sample_id)
  if (grepl("N$", short_id)) return("Normal")
  if (grepl("P$", short_id)) return("Primary")
  if (grepl("R$", short_id)) return("Recurrence")
  "Tumor"
}

IO_SIGNATURES <- list(
  T_cell_exhaustion = c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1"),
  T_cell_activation = c("IFNG", "GZMB", "PRF1", "TNF", "CD69", "IL2"),
  T_cell_cytotoxicity = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7"),
  Progenitor_CD8 = c("TCF7", "SLAMF6", "LEF1", "CCR7", "IL7R", "SELL"),
  Terminal_exhaustion = c("TOX", "ENTPD1", "HAVCR2", "LAYN", "TIGIT", "CXCL13"),
  Treg_signature = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  M1_macrophage = c("NOS2", "IL1B", "IL6", "TNF", "CD80", "CD86"),
  M2_macrophage = c("CD163", "MRC1", "MSR1", "CD200R1", "IL10"),
  IFN_gamma = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "IDO1"),
  TGF_beta = c("TGFB1", "TGFB2", "TGFB3", "ACTA2", "COL1A1", "SERPINE1"),
  Antigen_presentation = c("HLA.A", "HLA.B", "HLA.C", "B2M", "TAP1", "TAP2", "PSMB9"),
  Immune_checkpoint = c("CD274", "PDCD1LG2", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "VSIR"),
  Angiogenesis = c("VEGFA", "PECAM1", "KDR", "FLT1", "ANGPT2", "FGF2"),
  WNT_beta_catenin = c("CTNNB1", "APC", "AXIN2", "LGR5", "MYC"),
  mregDC = c("CCR7", "LAMP3", "CCL19", "CCL22", "CD80", "CD83", "FSCN1", "MARCKSL1"),
  TLS_signature = c("CXCL13", "CXCR5", "CR2", "CD79A", "MS4A1", "CD19", "IGHG1")
)

CELL_TYPE_MARKERS <- list(
  CD8_T = c("CD8A", "CD8B", "GZMB", "PRF1"),
  CD4_T = c("CD4", "IL7R", "CCR7"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  NK = c("NKG7", "GNLY", "KLRD1"),
  B_cell = c("CD79A", "MS4A1", "CD19"),
  Macrophage = c("CD68", "CD163", "C1QA", "CSF1R"),
  DC = c("CLEC9A", "XCR1", "BATF3"),
  Hepatocyte = c("ALB", "APOB", "HP", "TF"),
  Endothelial = c("PECAM1", "CDH5", "VWF"),
  Fibroblast = c("COL1A1", "COL3A1", "DCN", "LUM"),
  Tumor = c("GPC3", "AFP", "GOLM1", "EPCAM")
)

REFERENCE_CELL_TYPES <- c(
  "CD8_T", "CD4_T", "Treg", "NK", "B_cell", "Macrophage",
  "DC", "Endothelial", "Fibroblast"
)

IMMUNOTYPE_LEVELS <- c("Immune_Hot", "Immune_Excluded", "Immune_Desert", "Immune_Mixed")

IMMUNOTYPE_PALETTE <- c(
  Immune_Hot = "#CC3311",
  Immune_Excluded = "#F28E2B",
  Immune_Desert = "#0077BB",
  Immune_Mixed = "#6C757D"
)

RESPONSE_PALETTE <- c(Responder = "#0077BB", Non_Responder = "#CC3311")

TISSUE_PALETTE <- c(
  Normal = "#4E79A7",
  Tumor = "#E15759",
  Primary = "#E15759",
  Recurrence = "#B07AA1"
)

CELLTYPE_PALETTE <- c(
  CD8_T = "#2E8B57",
  CD4_T = "#66A61E",
  Treg = "#A6D854",
  NK = "#1F78B4",
  B_cell = "#6A3D9A",
  Macrophage = "#D95F02",
  DC = "#E6AB02",
  Hepatocyte = "#C0392B",
  Endothelial = "#377EB8",
  Fibroblast = "#8C564B",
  Tumor = "#7F0000",
  Mixed = "#7F7F7F"
)

KEY_SIGNATURES <- c(
  "T_cell_cytotoxicity",
  "T_cell_exhaustion",
  "TLS_signature",
  "M1_macrophage",
  "TGF_beta",
  "IFN_gamma"
)

CNV_FOCUS_ARMS <- c("1p", "1q", "4q", "5q", "6p", "8p", "8q", "13q", "16p", "16q", "17p", "17q", "20q")