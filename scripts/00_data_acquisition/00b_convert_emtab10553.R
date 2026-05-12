#!/usr/bin/env Rscript
# Convert E-MTAB-10553 (Pfister et al. 2021) Seurat .rds → MTX + CSV for Python
# NASH-HCC + anti-PD1 scRNA-seq data (humanLiver_scRNAseq_seuratObject_Wang2021.rds)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  project_dir <- normalizePath(file.path(dirname(script_path), "../.."))
} else {
  project_dir <- getwd()
}

data_dir <- file.path(project_dir, "data/pre_treatment/E-MTAB-10553_Pfister2021")
out_dir  <- file.path(project_dir, "data/processed/E-MTAB-10553")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(data_dir, "humanLiver_scRNAseq_seuratObject_Wang2021.rds")
if (!file.exists(rds_file)) {
  stop("RDS file not found: ", rds_file)
}

cat("Loading Seurat object...\n")
obj <- readRDS(rds_file)
cat("Class:", class(obj), "\n")
cat("Dimensions:", ncol(obj), "cells x", nrow(obj), "genes\n")

# Extract count matrix
cat("Extracting count matrix...\n")
if (inherits(obj, "Seurat")) {
  counts <- GetAssayData(obj, slot = "counts")
  meta   <- obj@meta.data
} else {
  stop("Unexpected object class: ", class(obj))
}

cat("Writing MTX...\n")
writeMM(counts, file.path(out_dir, "matrix.mtx"))
system(paste("gzip -f", file.path(out_dir, "matrix.mtx")))

cat("Writing barcodes...\n")
write.table(colnames(counts), file.path(out_dir, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
system(paste("gzip -f", file.path(out_dir, "barcodes.tsv")))

cat("Writing genes...\n")
write.table(data.frame(gene_id = rownames(counts), gene_name = rownames(counts)),
            file.path(out_dir, "genes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
system(paste("gzip -f", file.path(out_dir, "genes.tsv")))

cat("Writing metadata...\n")
write.csv(meta, file.path(out_dir, "metadata.csv"), row.names = TRUE)

cat("\nMetadata columns:\n")
for (col in colnames(meta)) {
  if (is.character(meta[[col]]) || is.factor(meta[[col]])) {
    vals <- table(meta[[col]])
    if (length(vals) <= 20) {
      cat("  ", col, ":", paste(names(vals), vals, sep = "=", collapse = ", "), "\n")
    } else {
      cat("  ", col, ":", length(vals), "unique values\n")
    }
  }
}

cat("\n[DONE] E-MTAB-10553 exported to:", out_dir, "\n")
system(paste("ls -lh", out_dir))
