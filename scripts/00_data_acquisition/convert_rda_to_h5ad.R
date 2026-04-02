#!/usr/bin/env Rscript
# Convert GSE206325 SingleCellExperiment .Rda to h5ad for Python
# Magen et al. 2023, Nat Med — neoadjuvant PD-1 blockade in HCC (~1M cells)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  project_dir <- normalizePath(file.path(dirname(script_path), "../.."))
} else {
  project_dir <- getwd()
}
data_dir <- file.path(project_dir, "data/pre_treatment/GSE206325_Liu2022_antiPD1")
out_dir  <- file.path(project_dir, "data/processed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading GSE206325 SCE object...\n")
rda_file <- file.path(data_dir, "GSE206325_data_SingleCellExperiment_object.Rda.gz")
if (!file.exists(rda_file)) {
  stop("File not found: ", rda_file)
}

# First decompress the .gz to get the .Rda
rda_unzipped <- sub("\\.gz$", "", rda_file)
if (!file.exists(rda_unzipped)) {
  cat("Decompressing .Rda.gz...\n")
  system(paste("gunzip -k", shQuote(rda_file)))
}

cat("Loading .Rda file...\n")
load(rda_unzipped)
# Find the SCE object
sce_name <- ls()[sapply(ls(), function(x) is(get(x), "SingleCellExperiment"))]
if (length(sce_name) == 0) {
  cat("No SCE found, checking all objects:\n")
  for (obj_name in ls()) {
    cat("  ", obj_name, ":", class(get(obj_name)), "\n")
  }
  stop("No SingleCellExperiment found")
}
cat("Found SCE:", sce_name[1], "\n")
sce <- get(sce_name[1])

cat("Dimensions:", dim(sce), "\n")
cat("Assays:", assayNames(sce), "\n")
cat("colData columns:", paste(colnames(colData(sce)), collapse=", "), "\n")

# Export count matrix
cat("Writing count matrix (MTX format)...\n")
counts_mat <- assay(sce, assayNames(sce)[1])
writeMM(counts_mat, file.path(out_dir, "GSE206325_matrix.mtx"))

# Export barcodes
write.table(colnames(sce), file.path(out_dir, "GSE206325_barcodes.tsv"),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# Export genes
gene_df <- data.frame(
  gene_id = rownames(sce),
  gene_name = if ("Symbol" %in% colnames(rowData(sce))) rowData(sce)$Symbol else rownames(sce)
)
write.table(gene_df, file.path(out_dir, "GSE206325_genes.tsv"),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# Export metadata
meta <- as.data.frame(colData(sce))
write.csv(meta, file.path(out_dir, "GSE206325_metadata.csv"), row.names=TRUE)

cat("colData summary:\n")
for (col in colnames(meta)) {
  if (is.character(meta[[col]]) || is.factor(meta[[col]])) {
    vals <- table(meta[[col]])
    if (length(vals) <= 20) {
      cat("  ", col, ":", paste(names(vals), vals, sep="=", collapse=", "), "\n")
    } else {
      cat("  ", col, ":", length(vals), "unique values\n")
    }
  } else {
    cat("  ", col, ": numeric, range", range(meta[[col]], na.rm=TRUE), "\n")
  }
}

# Also export TCR data
tcr_file <- file.path(data_dir, "GSE206325_revised_processed_TCR.csv.gz")
if (file.exists(tcr_file)) {
  tcr <- read.csv(gzfile(tcr_file))
  write.csv(tcr, file.path(out_dir, "GSE206325_tcr.csv"), row.names=FALSE)
  cat("TCR data:", nrow(tcr), "clonotypes exported\n")
}

# Export sample annotations
annot_file <- file.path(data_dir, "GSE206325_sample_annots_Liver_Treated_patients.csv.gz")
if (file.exists(annot_file)) {
  annots <- read.csv(gzfile(annot_file))
  write.csv(annots, file.path(out_dir, "GSE206325_sample_annotations.csv"), row.names=FALSE)
  cat("Sample annotations:", nrow(annots), "patients\n")
  cat("Columns:", paste(colnames(annots), collapse=", "), "\n")
  # Print response info
  for (col in colnames(annots)) {
    if (grepl("response|recist|pathol", col, ignore.case=TRUE)) {
      cat("  ", col, ":", paste(table(annots[[col]]), collapse=", "), "\n")
    }
  }
}

# Export cluster annotations
clust_file <- file.path(data_dir, "GSE206325_full_HCC_cluster_annotation.csv.gz")
if (file.exists(clust_file)) {
  clusters <- read.csv(gzfile(clust_file))
  write.csv(clusters, file.path(out_dir, "GSE206325_cluster_annotations.csv"), row.names=FALSE)
  cat("Cluster annotations:", nrow(clusters), "entries\n")
}

# Compress large files
cat("Compressing matrix...\n")
system(paste("gzip -f", file.path(out_dir, "GSE206325_matrix.mtx")))

cat("\n[DONE] GSE206325 exported to:", out_dir, "\n")
cat("Files:\n")
system(paste("ls -lh", out_dir))
