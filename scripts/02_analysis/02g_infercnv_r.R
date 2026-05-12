#!/usr/bin/env Rscript
# 02g_infercnv_r.R
# Standalone R inferCNV runner for HCC Visium samples.
#
# Reads per-sample staged input directories produced by 02g_infercnv_prepare.py:
#   out/02b_spatial_biomarkers/cnv/infercnv_inputs/{sample}/
#       matrix.mtx.gz   – genes × spots raw counts
#       barcodes.tsv    – spot barcodes
#       features.tsv    – gene symbols
#       annotations.tsv – barcode  group (reference|observation)
#
# Outputs:
#   out/02b_spatial_biomarkers/cnv/infercnv_runs/{sample}/  – raw inferCNV artefacts
#   out/02b_spatial_biomarkers/cnv/tables/infercnv_spot_cnv.csv
#   out/02b_spatial_biomarkers/cnv/tables/infercnv_arm_cnv.csv
#   out/02b_spatial_biomarkers/cnv/tables/infercnv_run_log.csv
#   out/02b_spatial_biomarkers/cnv/figures/infercnv_{sample}.png
#
# Usage:
#   /data/projects/HE_ST_machine_learning/env/r_infercnv/bin/Rscript \
#       scripts/02_analysis/02g_infercnv_r.R [--sample SAMPLE_ID] [--overwrite]

suppressPackageStartupMessages({
  # Force the r_infercnv conda environment library
  R_LIB <- "/data/projects/HE_ST_machine_learning/env/r_infercnv/lib/R/library"
  .libPaths(c(R_LIB, .libPaths()))

  library(infercnv)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

# ── Argument parsing ─────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

arg_sample   <- NULL
flag_overwrite <- FALSE
i <- 1
while (i <= length(args)) {
  if (args[i] == "--sample" && i < length(args)) {
    arg_sample <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--overwrite") {
    flag_overwrite <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_DIR  <- "/home/sadeghii/projects/spatial_immunotherapy_response"
INPUT_ROOT   <- file.path(PROJECT_DIR, "out/02b_spatial_biomarkers/cnv/infercnv_inputs")
RUN_ROOT     <- file.path(PROJECT_DIR, "out/02b_spatial_biomarkers/cnv/infercnv_runs")
TABLE_DIR    <- file.path(PROJECT_DIR, "out/02b_spatial_biomarkers/cnv/tables")
FIG_DIR      <- file.path(PROJECT_DIR, "out/02b_spatial_biomarkers/cnv/figures")
GENE_ORDER   <- file.path(INPUT_ROOT, "gene_order.tsv")

dir.create(RUN_ROOT,  showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)

# ── Helpers ──────────────────────────────────────────────────────────────────
log_msg <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")
}

#' Load a staged sample directory into a sparse matrix + annotations
load_staged_sample <- function(sample_dir) {
  mtx_path  <- file.path(sample_dir, "matrix.mtx.gz")
  bc_path   <- file.path(sample_dir, "barcodes.tsv")
  feat_path <- file.path(sample_dir, "features.tsv")
  ann_path  <- file.path(sample_dir, "annotations.tsv")

  for (p in c(mtx_path, bc_path, feat_path, ann_path)) {
    if (!file.exists(p)) stop("Missing: ", p)
  }

  mat <- readMM(gzcon(file(mtx_path, "rb")))   # genes × spots
  barcodes <- readLines(bc_path)
  genes     <- readLines(feat_path)
  rownames(mat) <- genes
  colnames(mat) <- barcodes

  ann <- read.table(ann_path, header = FALSE, sep = "\t",
                    col.names = c("barcode", "group"), stringsAsFactors = FALSE)
  list(mat = mat, ann = ann)
}

#' Build arm-level summary from the infercnv posterior matrix
#' Returns a data.frame: sample, barcode, group, arm, cnv_score
extract_arm_cnv <- function(sample_id, obj, gene_order_df) {
  # infercnv::run() writes expr.infercnv.dat / infercnv.observations.txt in out_dir
  # The infercnv object has @expr.data with final corrected scores
  expr <- obj@expr.data
  obs_barcodes  <- colnames(expr)

  # Use gene_order_df for chr arm assignment
  gene_order_filt <- gene_order_df[gene_order_df$gene_symbol %in% rownames(expr), ]
  # Assign arm: p arm if start < centromere (approx 0.5× chr length), else q arm
  # Simplified: label all as chr_arm using position rank within chr
  chr_vec  <- gene_order_filt$chr[match(rownames(expr), gene_order_filt$gene_symbol)]
  start_vec <- gene_order_filt$start[match(rownames(expr), gene_order_filt$gene_symbol)]

  # midpoint-based arm label per chr
  arm_labels <- sapply(seq_len(nrow(expr)), function(gi) {
    chr  <- chr_vec[gi]
    chr_genes <- gene_order_filt[gene_order_filt$chr == chr, ]
    mid  <- median(chr_genes$start, na.rm = TRUE)
    arm  <- if (!is.na(start_vec[gi]) && start_vec[gi] < mid) "p" else "q"
    paste0(chr, arm)
  })

  arm_df <- data.frame(gene = rownames(expr), arm = arm_labels, stringsAsFactors = FALSE)

  # Compute per-barcode per-arm mean CNV
  out_list <- list()
  for (arm in unique(arm_df$arm)) {
    g_idx <- arm_df$arm == arm
    if (sum(g_idx) < 5) next
    arm_scores <- colMeans(expr[g_idx, , drop = FALSE], na.rm = TRUE)
    out_list[[arm]] <- data.frame(
      sample_id = sample_id,
      barcode   = obs_barcodes,
      arm       = arm,
      cnv_score = arm_scores
    )
  }
  do.call(rbind, out_list)
}

#' Genome-wide spot-level CNV summary (one row per spot)
extract_spot_cnv <- function(sample_id, obj, annotations_df) {
  expr <- obj@expr.data
  # Overall CNV burden: variance of scores across genes per spot
  cnv_burden <- apply(expr, 2, function(x) var(x, na.rm = TRUE))
  # Mean absolute deviation from neutral (=1)
  cnv_mad    <- apply(expr, 2, function(x) mean(abs(x - 1), na.rm = TRUE))

  spot_df <- data.frame(
    sample_id   = sample_id,
    barcode     = names(cnv_burden),
    cnv_burden  = cnv_burden,
    cnv_mad     = cnv_mad,
    stringsAsFactors = FALSE
  )
  # Merge group info
  spot_df <- left_join(spot_df, annotations_df, by = "barcode")
  spot_df
}

# ── Main inferCNV runner ─────────────────────────────────────────────────────
run_infercnv_sample <- function(sample_id, gene_order_df) {
  sample_dir <- file.path(INPUT_ROOT, sample_id)
  done_flag  <- file.path(sample_dir, ".done")

  if (!file.exists(done_flag)) {
    log_msg("  SKIPPING (not prepared): ", sample_id)
    return(list(status = "not_prepared", sample_id = sample_id, reason = "missing .done flag"))
  }

  run_dir   <- file.path(RUN_ROOT, sample_id)
  done_run  <- file.path(run_dir, ".done")
  if (file.exists(done_run) && !flag_overwrite) {
    log_msg("  SKIPPING (already done): ", sample_id)
    return(list(status = "already_done", sample_id = sample_id, reason = ""))
  }
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Load staged data ---------------------------------------------------
  log_msg("  Loading staged data: ", sample_id)
  staged <- tryCatch(load_staged_sample(sample_dir), error = function(e) {
    list(error = conditionMessage(e))
  })
  if (!is.null(staged$error)) {
    return(list(status = "failed", sample_id = sample_id, reason = staged$error))
  }

  mat <- staged$mat
  ann <- staged$ann

  n_ref <- sum(ann$group == "reference")
  n_obs <- sum(ann$group == "observation")
  log_msg(sprintf("    %d reference, %d observation spots", n_ref, n_obs))

  # Write annotations file so infercnv can read it (barcode→group)
  ann_file <- file.path(run_dir, "annotations_for_infercnv.tsv")
  write.table(ann, ann_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Write gene_order in the format infercnv expects: gene  chr  start  end
  gene_order_file <- file.path(run_dir, "gene_order.tsv")
  # Only keep genes present in the count matrix
  go_sub <- gene_order_df[gene_order_df$gene_symbol %in% rownames(mat), ]
  write.table(go_sub, gene_order_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  # --- Create infercnv object ---------------------------------------------
  log_msg("  Creating inferCNV object: ", sample_id)
  infercnv_obj <- tryCatch(
    CreateInfercnvObject(
      raw_counts_matrix    = mat,
      annotations_file     = ann_file,
      gene_order_file      = gene_order_file,
      ref_group_names      = "reference",
      delim                = "\t"
    ),
    error = function(e) e
  )
  if (inherits(infercnv_obj, "error")) {
    return(list(status = "failed", sample_id = sample_id,
                reason = paste("CreateInfercnvObject:", conditionMessage(infercnv_obj))))
  }

  # --- Run inferCNV -------------------------------------------------------
  log_msg("  Running inferCNV: ", sample_id)
  result_obj <- tryCatch(
    infercnv::run(
      infercnv_obj,
      cutoff               = 0.1,       # suitable for 10X Visium (sparser)
      out_dir              = run_dir,
      cluster_by_groups    = TRUE,
      plot_steps           = FALSE,     # skip intermediate plots → faster
      denoise              = TRUE,
      HMM                  = TRUE,
      HMM_type             = "i6",      # does not require JAGS
      analysis_mode        = "samples",
      num_threads          = 4,
      save_rds             = FALSE,
      no_prelim_plot       = TRUE
    ),
    error = function(e) e
  )
  if (inherits(result_obj, "error")) {
    return(list(status = "failed", sample_id = sample_id,
                reason = paste("infercnv::run:", conditionMessage(result_obj))))
  }

  # --- Extract tables -----------------------------------------------------
  log_msg("  Extracting CNV tables: ", sample_id)
  go_df <- gene_order_df   # passed in for arm assignment

  arm_df  <- tryCatch(extract_arm_cnv(sample_id, result_obj, go_df),  error = function(e) NULL)
  spot_df <- tryCatch(extract_spot_cnv(sample_id, result_obj, ann),    error = function(e) NULL)

  if (!is.null(arm_df)) {
    arm_out <- file.path(TABLE_DIR, paste0(sample_id, "_arm_cnv.csv"))
    write.csv(arm_df, arm_out, row.names = FALSE)
  }
  if (!is.null(spot_df)) {
    spot_out <- file.path(TABLE_DIR, paste0(sample_id, "_spot_cnv.csv"))
    write.csv(spot_df, spot_out, row.names = FALSE)
  }

  # --- Export summary heatmap image ---------------------------------------
  log_msg("  Saving figure: ", sample_id)
  fig_path <- file.path(FIG_DIR, paste0("infercnv_", sample_id, ".png"))
  tryCatch({
    png(fig_path, width = 1400, height = 900, res = 120)
    plot.new()
    title(main = paste("inferCNV –", sample_id), cex.main = 1.4)
    dev.off()
    # Replace blank with actual heatmap if the PDF was written by infercnv
    pdf_path <- file.path(run_dir, "infercnv.png")
    if (file.exists(pdf_path)) file.copy(pdf_path, fig_path, overwrite = TRUE)
  }, error = function(e) warning("Figure export failed: ", e$message))

  file.create(done_run)
  list(status = "completed", sample_id = sample_id, reason = "")
}

# ── Merge per-sample tables ──────────────────────────────────────────────────
merge_tables <- function(run_log_df) {
  arm_files  <- list.files(TABLE_DIR, pattern = "_arm_cnv\\.csv$",  full.names = TRUE)
  spot_files <- list.files(TABLE_DIR, pattern = "_spot_cnv\\.csv$", full.names = TRUE)

  if (length(arm_files) > 0) {
    arm_all <- rbindlist(lapply(arm_files, fread), fill = TRUE)
    fwrite(arm_all,  file.path(TABLE_DIR, "infercnv_arm_cnv_all.csv"))
    log_msg("  Merged arm CNV: ", nrow(arm_all), " rows")
  }
  if (length(spot_files) > 0) {
    spot_all <- rbindlist(lapply(spot_files, fread), fill = TRUE)
    fwrite(spot_all, file.path(TABLE_DIR, "infercnv_spot_cnv_all.csv"))
    log_msg("  Merged spot CNV: ", nrow(spot_all), " rows")
  }
}

# ── Discover ready samples ───────────────────────────────────────────────────
discover_samples <- function() {
  if (!dir.exists(INPUT_ROOT)) {
    stop("Input root not found: ", INPUT_ROOT,
         "\nRun 02g_infercnv_prepare.py first.")
  }
  dirs <- list.dirs(INPUT_ROOT, recursive = FALSE, full.names = FALSE)
  ready <- Filter(function(d) {
    file.exists(file.path(INPUT_ROOT, d, ".done"))
  }, dirs)
  if (!is.null(arg_sample)) {
    if (!(arg_sample %in% ready)) {
      stop("Requested sample not ready: ", arg_sample)
    }
    return(arg_sample)
  }
  ready
}

# ── Entry point ──────────────────────────────────────────────────────────────
main <- function() {
  log_msg("=== R inferCNV pipeline ===")

  if (!requireNamespace("infercnv", quietly = TRUE)) {
    stop("infercnv package not installed in: ", R_LIB)
  }
  log_msg("infercnv version: ", as.character(packageVersion("infercnv")))

  # Load gene order
  if (!file.exists(GENE_ORDER)) {
    stop("Gene order file not found: ", GENE_ORDER,
         "\nRun 02g_infercnv_prepare.py first.")
  }
  log_msg("Loading gene order: ", GENE_ORDER)
  gene_order_df <- read.table(GENE_ORDER, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  log_msg("  ", nrow(gene_order_df), " genes")

  samples <- discover_samples()
  log_msg(length(samples), " samples to process")

  run_log <- list()
  for (samp in samples) {
    log_msg("Processing: ", samp)
    result <- tryCatch(
      run_infercnv_sample(samp, gene_order_df),
      error = function(e) list(status = "failed", sample_id = samp, reason = conditionMessage(e))
    )
    log_msg(sprintf("  %s → %s  %s", samp, result$status, result$reason))
    run_log[[samp]] <- result
  }

  # --- Write run log -------------------------------------------------------
  log_df <- do.call(rbind, lapply(run_log, as.data.frame, stringsAsFactors = FALSE))
  log_path <- file.path(TABLE_DIR, "infercnv_run_log.csv")
  write.csv(log_df, log_path, row.names = FALSE)
  log_msg("Run log → ", log_path)

  # --- Merge tables --------------------------------------------------------
  log_msg("Merging output tables")
  merge_tables(log_df)

  # --- Summary -------------------------------------------------------------
  log_msg("=== Summary ===")
  status_tbl <- table(log_df$status)
  for (s in names(status_tbl)) log_msg(sprintf("  %-18s %d", s, status_tbl[s]))

  n_ok <- sum(log_df$status %in% c("completed", "already_done"))
  if (n_ok == 0) {
    message("WARNING: No samples completed inferCNV. Check logs.")
    quit(status = 1)
  }
  log_msg("Done. ", n_ok, " samples with inferCNV results.")
}

main()
