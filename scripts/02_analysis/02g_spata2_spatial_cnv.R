#!/usr/bin/env Rscript

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_path <- normalizePath(sub("^--file=", "", .script_arg[[1]]), mustWork = FALSE)
options(spata2_script_dir = file.path(dirname(.script_path), "spata2"))
source(file.path(getOption("spata2_script_dir"), "utils.R"))

set.seed(SPATA_SEED)
ensure_dirs(CNV_ROOT, CNV_FIG_DIR, CNV_TAB_DIR, CNV_RDS_DIR, CNV_QC_DIR, NOTE_DIR)

cat("\n============================================================\n")
cat("SPATA2 spatial CNV reanalysis\n")
cat("============================================================\n")

manifest <- make_input_manifest()
write_input_manifest(manifest)
stop_on_missing_inputs(manifest)

if (!file.exists(GTF_PATH)) {
  stop(sprintf("Gene annotation file not found: %s", GTF_PATH), call. = FALSE)
}

sample_rows <- manifest %>% filter(status == "present")
sample_objects <- vector("list", nrow(sample_rows))
names(sample_objects) <- sample_rows$sample_id

for (i in seq_len(nrow(sample_rows))) {
  sample_id <- sample_rows$sample_id[[i]]
  dataset <- sample_rows$dataset[[i]]
  group_label <- sample_rows$group[[i]]
  path <- sample_rows$expected_path[[i]]
  cat(sprintf("  [%02d/%02d] %s (%s)\n", i, nrow(sample_rows), sample_id, dataset))
  sample_objects[[sample_id]] <- load_spatial_sample(path = path, sample_id = sample_id, dataset = dataset, group_label = group_label)
}

all_genes <- unique(unlist(lapply(sample_objects, function(x) rownames(x$counts))))
gene_order_df <- gene_order_from_gtf(GTF_PATH, all_genes)
if (nrow(gene_order_df) == 0) {
  stop("No overlapping genes were found between the spatial datasets and the GTF gene order.", call. = FALSE)
}
write_csv(gene_order_df, file.path(CNV_TAB_DIR, "spatial_cnv_gene_order.csv"))

arm_tables <- list()
burden_tables <- list()
run_log <- list()
sample_qc <- list()

for (sample_id in names(sample_objects)) {
  sample_obj <- sample_objects[[sample_id]]
  sample_run_dir <- file.path(CNV_RDS_DIR, sample_id)
  sample_prefix <- file.path(sample_run_dir, sample_id)
  ensure_dirs(sample_run_dir)

  sample_gene_order <- gene_order_df %>% filter(gene_symbol %in% rownames(sample_obj$counts))
  ref_group <- infer_reference_group(sample_obj$spot_df)
  qc_row <- tibble(
    sample_id = sample_id,
    dataset = sample_obj$dataset,
    group = sample_obj$group,
    n_spots = ncol(sample_obj$counts),
    n_reference = sum(ref_group == "reference"),
    n_observation = sum(ref_group == "observation"),
    n_genes_overlap = nrow(sample_gene_order),
    status = "pending",
    reason = NA_character_
  )

  if (nrow(sample_gene_order) < 500) {
    qc_row$status <- "skipped"
    qc_row$reason <- "too_few_genes_after_gtf_overlap"
    sample_qc[[sample_id]] <- qc_row
    run_log[[sample_id]] <- tibble(sample_id = sample_id, status = "skipped", reason = qc_row$reason)
    next
  }
  if (qc_row$n_reference < 20 || qc_row$n_observation < 20) {
    qc_row$status <- "skipped"
    qc_row$reason <- "insufficient_reference_or_observation_spots"
    sample_qc[[sample_id]] <- qc_row
    run_log[[sample_id]] <- tibble(sample_id = sample_id, status = "skipped", reason = qc_row$reason)
    next
  }

  run_status <- tryCatch({
    run_spatial_infercnv(sample_obj, out_prefix = sample_prefix, gene_order_df = sample_gene_order, cutoff = 0.1)
    final_rds_candidates <- c(
      file.path(sample_run_dir, "run.final.infercnv_obj"),
      file.path(sample_run_dir, "run.final.infercnv_obj.rds"),
      file.path(sample_run_dir, "infercnv.observations.txt"),
      file.path(sample_run_dir, paste0(sample_id, "_infercnv_final.rds"))
    )
    final_rds <- final_rds_candidates[file.exists(final_rds_candidates)][1]
    if (is.na(final_rds) || !grepl("infercnv_obj|infercnv_final", basename(final_rds))) {
      final_rds <- file.path(sample_run_dir, "run.final.infercnv_obj")
      if (!file.exists(final_rds)) final_rds <- file.path(sample_run_dir, "run.final.infercnv_obj.rds")
    }
    if (!file.exists(final_rds)) stop("infercnv completed without a final RDS output")

    infer_obj <- readRDS(final_rds)
    saveRDS(infer_obj, file.path(sample_run_dir, paste0(sample_id, "_infercnv_final.rds")))
    expr_mtr <- infer_obj@expr.data
    arm_tbl <- infercnv_expr_to_arm_table(expr_mtr, sample_gene_order, sample_id)
    if (is.null(arm_tbl) || nrow(arm_tbl) == 0) stop("no arm-level CNV table could be constructed")

    burden_tbl <- arm_tbl %>%
      group_by(sample_id, barcode) %>%
      summarise(cnv_burden = mean(abs(cnv), na.rm = TRUE), .groups = "drop") %>%
      left_join(
        sample_obj$spot_df %>%
          select(barcodes, group, dominant_cell_type, spatial_immunotype, annotation_confidence, x, y),
        by = c("barcode" = "barcodes")
      )

    arm_tables[[sample_id]] <- arm_tbl
    burden_tables[[sample_id]] <- burden_tbl
    qc_row$status <- "completed"
    sample_qc[[sample_id]] <- qc_row
    tibble(sample_id = sample_id, status = "completed", reason = NA_character_)
  }, error = function(e) {
    qc_row$status <- "failed"
    qc_row$reason <- conditionMessage(e)
    sample_qc[[sample_id]] <- qc_row
    tibble(sample_id = sample_id, status = "failed", reason = conditionMessage(e))
  })

  run_log[[sample_id]] <- run_status
}

run_log_df <- bind_rows(run_log)
sample_qc_df <- bind_rows(sample_qc)
write_csv(run_log_df, file.path(CNV_TAB_DIR, "spatial_cnv_run_log.csv"))
write_csv(sample_qc_df, file.path(CNV_QC_DIR, "spatial_cnv_sample_qc.csv"))

if (!any(run_log_df$status == "completed")) {
  stop(
    paste0(
      "Spatial CNV code has been staged, but no sample completed infercnv. See ",
      file.path(CNV_TAB_DIR, "spatial_cnv_run_log.csv"), "."
    ),
    call. = FALSE
  )
}

spot_arm_tbl <- bind_rows(arm_tables)
spot_burden_tbl <- bind_rows(burden_tables)

write_csv(spot_arm_tbl, file.path(CNV_TAB_DIR, "spatial_cnv_spot_arm.csv"))
write_csv(spot_burden_tbl, file.path(CNV_TAB_DIR, "spatial_cnv_spot_burden.csv"))

arm_summary <- spot_arm_tbl %>%
  group_by(sample_id, arm) %>%
  summarise(mean_cnv = mean(cnv, na.rm = TRUE), .groups = "drop")
write_csv(arm_summary, file.path(CNV_TAB_DIR, "spatial_cnv_arm_summary.csv"))

arm_heatmap_tbl <- arm_summary %>%
  left_join(
    bind_rows(lapply(sample_objects, function(so) tibble(sample_id = so$sample_id, response = so$group))) %>% distinct(),
    by = "sample_id"
  )
plot_arm_heatmap(arm_heatmap_tbl, "02g_spatial_cnv_arm_heatmap")

burden_summary <- spot_burden_tbl %>%
  group_by(sample_id, group, dominant_cell_type, spatial_immunotype) %>%
  summarise(
    mean_cnv_burden = mean(cnv_burden, na.rm = TRUE),
    median_cnv_burden = median(cnv_burden, na.rm = TRUE),
    n_spots = n(),
    .groups = "drop"
  )
write_csv(burden_summary, file.path(CNV_TAB_DIR, "spatial_cnv_burden_summary.csv"))

response_burden <- spot_burden_tbl %>%
  filter(group %in% c("Responder", "Non_Responder")) %>%
  group_by(sample_id, group) %>%
  summarise(
    mean_cnv_burden = mean(cnv_burden, na.rm = TRUE),
    median_cnv_burden = median(cnv_burden, na.rm = TRUE),
    n_spots = n(),
    .groups = "drop"
  )
write_csv(response_burden, file.path(CNV_TAB_DIR, "gse238264_cnv_burden_by_sample.csv"))

response_burden_stats <- adjust_stat_table(pairwise_group_stats(response_burden, "group", c("Responder", "Non_Responder"), c("mean_cnv_burden", "median_cnv_burden")))
write_csv(response_burden_stats, file.path(CNV_TAB_DIR, "gse238264_cnv_burden_stats.csv"))

focus_arm_summary <- arm_summary %>%
  filter(arm %in% CNV_FOCUS_ARMS) %>%
  left_join(sample_qc_df %>% select(sample_id, group), by = "sample_id")
write_csv(focus_arm_summary, file.path(CNV_TAB_DIR, "spatial_cnv_focus_arm_summary.csv"))

gain_loss_calls <- arm_summary %>%
  mutate(
    cnv_state = case_when(
      mean_cnv >= 0.03 ~ "gain",
      mean_cnv <= -0.03 ~ "loss",
      TRUE ~ "neutral"
    )
  )
write_csv(gain_loss_calls, file.path(CNV_TAB_DIR, "spatial_cnv_gain_loss_calls.csv"))

completed_ids <- run_log_df$sample_id[run_log_df$status == "completed"]
infercnv_list <- list()
spot_df_list <- list()
for (sid in completed_ids) {
  rds_path <- file.path(CNV_RDS_DIR, sid, paste0(sid, "_infercnv_final.rds"))
  if (file.exists(rds_path)) {
    infercnv_list[[sid]] <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    spot_df_list[[sid]] <- sample_objects[[sid]]$spot_df
  }
}

if (length(infercnv_list) > 0) {
  cat("\nGenerating genome-wide CNV heatmaps ...\n")
  plot_genome_cnv_heatmap(
    infercnv_list = infercnv_list,
    spot_df_list = spot_df_list,
    gene_order_df = gene_order_df,
    file_stub = "02g_spatial_cnv_genome_heatmap_celltype",
    group_col = "dominant_cell_type",
    plot_width_px = 9000L,
    plot_height_px = 2600L,
    plot_dpi = 300L
  )

  gse238264_infer <- infercnv_list[names(infercnv_list) %in% names(GSE238264_SAMPLES)]
  gse238264_spots_cnv <- spot_df_list[names(spot_df_list) %in% names(GSE238264_SAMPLES)]
  if (length(gse238264_infer) > 0) {
    plot_genome_cnv_heatmap(
      infercnv_list = gse238264_infer,
      spot_df_list = gse238264_spots_cnv,
      gene_order_df = gene_order_df,
      file_stub = "02g_spatial_cnv_genome_heatmap_response",
      group_col = "group",
      plot_width_px = 9000L,
      plot_height_px = 2200L,
      plot_dpi = 300L
    )
    plot_binned_cnv_heatmap(
      infercnv_list = gse238264_infer,
      spot_df_list = gse238264_spots_cnv,
      gene_order_df = gene_order_df,
      file_stub = "02g_spatial_cnv_binned_heatmap_response",
      max_spots_per_sample = 300L,
      bin_size = 60L
    )
  }

  plot_binned_cnv_heatmap(
    infercnv_list = infercnv_list,
    spot_df_list = spot_df_list,
    gene_order_df = gene_order_df,
    file_stub = "02g_spatial_cnv_binned_heatmap_all_samples",
    max_spots_per_sample = 220L,
    bin_size = 60L
  )
}

cnv_map_samples <- unique(spot_burden_tbl$sample_id)[seq_len(min(4, n_distinct(spot_burden_tbl$sample_id)))]
cnv_map_plots <- lapply(cnv_map_samples, function(sample_id) {
  plot_spatial_numeric(
    spot_burden_tbl %>% filter(sample_id == .env$sample_id),
    feature = "cnv_burden",
    title = sample_id,
    subtitle = "InferCNV spatial CNA burden",
    palette = "magma"
  )
})
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_spatial_maps.png"), wrap_plots(cnv_map_plots, ncol = 2), width = 11, height = 9, dpi = 500)

burden_plot <- ggplot(spot_burden_tbl, aes(x = dominant_cell_type, y = cnv_burden, fill = dominant_cell_type)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = CELLTYPE_PALETTE, drop = FALSE) +
  labs(title = "Spatial CNV burden by dominant cell type", x = NULL, y = "Mean absolute CNV deviation", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_burden_by_celltype.png"), burden_plot, width = 12, height = 6, dpi = 500)

response_plot <- ggplot(response_burden, aes(x = group, y = mean_cnv_burden, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2) +
  scale_fill_manual(values = RESPONSE_PALETTE, drop = FALSE) +
  labs(title = "GSE238264 sample-level CNV burden", x = NULL, y = "Mean CNV burden", fill = NULL) +
  theme_pub()
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_response_burden.png"), response_plot, width = 7, height = 5, dpi = 500)

focus_arm_plot <- focus_arm_summary %>%
  ggplot(aes(x = arm, y = mean_cnv, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(RESPONSE_PALETTE, TISSUE_PALETTE), drop = FALSE) +
  labs(title = "Focused arm-level CNV summary", x = NULL, y = "Mean CNV deviation", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_focus_arms.png"), focus_arm_plot, width = 14, height = 6, dpi = 500)

gain_loss_plot <- gain_loss_calls %>%
  filter(cnv_state != "neutral") %>%
  count(arm, cnv_state, name = "n_samples") %>%
  ggplot(aes(x = arm, y = n_samples, fill = cnv_state)) +
  geom_col(position = position_stack()) +
  scale_fill_manual(values = c(gain = "#B2182B", loss = "#2166AC"), drop = FALSE) +
  labs(title = "Arm-level gain/loss frequency", x = NULL, y = "Samples", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_gain_loss_frequency.png"), gain_loss_plot, width = 14, height = 6, dpi = 500)

write_note_lines(
  file.path(NOTE_DIR, "02g_spata2_spatial_cnv_note.txt"),
  c(
    "Question: Re-run spatial inferCNV across the Visium cohorts and export higher-resolution CNV summaries and figures.",
    "Method: Per-sample inferCNV using SPATA2-derived reference/observation groups, arm-level summaries, and genome-binned ComplexHeatmap renders.",
    sprintf("Completed samples: %d of %d.", sum(run_log_df$status == "completed"), nrow(run_log_df)),
    "Key outputs: run log, CNV QC table, spot arm/burden tables, arm heatmaps, binned genome heatmaps, and response/cell-type burden figures.",
    "QC status: per-sample reference and observation counts exported before model execution."
  )
)

cat("\nSpatial CNV reanalysis complete. Outputs written to:\n")
cat("  Figures:", CNV_FIG_DIR, "\n")
cat("  QC:", CNV_QC_DIR, "\n")
cat("  Tables:", CNV_TAB_DIR, "\n")#!/usr/bin/env Rscript

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_path <- normalizePath(sub("^--file=", "", .script_arg[[1]]), mustWork = FALSE)
options(spata2_script_dir = file.path(dirname(.script_path), "spata2"))
source(file.path(getOption("spata2_script_dir"), "utils.R"))

set.seed(SPATA_SEED)
ensure_dirs(CNV_ROOT, CNV_FIG_DIR, CNV_TAB_DIR, CNV_RDS_DIR, NOTE_DIR)

cat("\n============================================================\n")
cat("SPATA2 spatial CNV reanalysis\n")
cat("============================================================\n")

manifest <- make_input_manifest()
write_input_manifest(manifest)
stop_on_missing_inputs(manifest)

if (!file.exists(GTF_PATH)) {
  stop(sprintf("Gene annotation file not found: %s", GTF_PATH), call. = FALSE)
}

sample_rows <- manifest %>% filter(status == "present")

sample_objects <- vector("list", nrow(sample_rows))
names(sample_objects) <- sample_rows$sample_id

for (i in seq_len(nrow(sample_rows))) {
  sample_id <- sample_rows$sample_id[[i]]
  dataset <- sample_rows$dataset[[i]]
  group_label <- sample_rows$group[[i]]
  path <- sample_rows$expected_path[[i]]
  cat(sprintf("  [%02d/%02d] %s (%s)\n", i, nrow(sample_rows), sample_id, dataset))
  sample_objects[[sample_id]] <- load_spatial_sample(path = path, sample_id = sample_id, dataset = dataset, group_label = group_label)
}

all_genes <- unique(unlist(lapply(sample_objects, function(x) rownames(x$counts))))
gene_order_df <- gene_order_from_gtf(GTF_PATH, all_genes)
if (nrow(gene_order_df) == 0) {
  stop("No overlapping genes were found between the spatial datasets and the GTF gene order.", call. = FALSE)
}

write_csv(gene_order_df, file.path(CNV_TAB_DIR, "spatial_cnv_gene_order.csv"))

arm_tables <- list()
burden_tables <- list()
run_log <- list()

for (sample_id in names(sample_objects)) {
  sample_obj <- sample_objects[[sample_id]]
  out_prefix <- file.path(CNV_RDS_DIR, sample_id)
  sample_gene_order <- gene_order_df %>% filter(gene_symbol %in% rownames(sample_obj$counts))

  if (nrow(sample_gene_order) < 500) {
    run_log[[sample_id]] <- tibble(sample_id = sample_id, status = "skipped", reason = "too_few_genes_after_gtf_overlap")
    next
  }

  run_status <- tryCatch({
    run_spatial_infercnv(sample_obj, out_prefix = out_prefix, gene_order_df = sample_gene_order, cutoff = 0.1)
    final_rds <- file.path(dirname(out_prefix), "run.final.infercnv_obj")
    if (!file.exists(final_rds)) final_rds <- file.path(dirname(out_prefix), "run.final.infercnv_obj.rds")
    if (!file.exists(final_rds)) stop("infercnv completed without a final RDS output")

    infer_obj <- readRDS(final_rds)
    saveRDS(infer_obj, file.path(CNV_RDS_DIR, paste0(sample_id, "_infercnv_final.rds")))
    expr_mtr <- infer_obj@expr.data
    arm_tbl <- infercnv_expr_to_arm_table(expr_mtr, sample_gene_order, sample_id)
    if (is.null(arm_tbl) || nrow(arm_tbl) == 0) stop("no arm-level CNV table could be constructed")

    burden_tbl <- arm_tbl %>%
      group_by(sample_id, barcode) %>%
      summarise(cnv_burden = mean(abs(cnv), na.rm = TRUE), .groups = "drop") %>%
      left_join(sample_obj$spot_df %>% select(barcodes, group, dominant_cell_type, x, y), by = c("barcode" = "barcodes"))

    arm_tables[[sample_id]] <<- arm_tbl
    burden_tables[[sample_id]] <<- burden_tbl
    tibble(sample_id = sample_id, status = "completed", reason = NA_character_)
  }, error = function(e) {
    tibble(sample_id = sample_id, status = "failed", reason = conditionMessage(e))
  })

  run_log[[sample_id]] <- run_status
}

run_log_df <- bind_rows(run_log)
write_csv(run_log_df, file.path(CNV_TAB_DIR, "spatial_cnv_run_log.csv"))

if (!any(run_log_df$status == "completed")) {
  stop(
    paste0(
      "Spatial CNV code has been staged, but no sample completed infercnv. See ",
      file.path(CNV_TAB_DIR, "spatial_cnv_run_log.csv"), "."
    ),
    call. = FALSE
  )
}

spot_arm_tbl <- bind_rows(arm_tables)
spot_burden_tbl <- bind_rows(burden_tables)

write_csv(spot_arm_tbl, file.path(CNV_TAB_DIR, "spatial_cnv_spot_arm.csv"))
write_csv(spot_burden_tbl, file.path(CNV_TAB_DIR, "spatial_cnv_spot_burden.csv"))

arm_summary <- spot_arm_tbl %>%
  group_by(sample_id, arm) %>%
  summarise(mean_cnv = mean(cnv, na.rm = TRUE), .groups = "drop")
write_csv(arm_summary, file.path(CNV_TAB_DIR, "spatial_cnv_arm_summary.csv"))

burden_summary <- spot_burden_tbl %>%
  group_by(sample_id, group, dominant_cell_type) %>%
  summarise(
    mean_cnv_burden = mean(cnv_burden, na.rm = TRUE),
    median_cnv_burden = median(cnv_burden, na.rm = TRUE),
    n_spots = n(),
    .groups = "drop"
  )
write_csv(burden_summary, file.path(CNV_TAB_DIR, "spatial_cnv_burden_summary.csv"))

arm_heatmap_tbl <- arm_summary %>%
  left_join(
    bind_rows(lapply(sample_objects, function(so) tibble(sample_id = so$sample_id, response = so$group))) %>% distinct(),
    by = "sample_id"
  )
if (nrow(arm_heatmap_tbl) > 0) {
  plot_arm_heatmap(arm_heatmap_tbl, "02g_spatial_cnv_arm_heatmap")
}

# --- Gene-level full-genome CNV heatmap using ComplexHeatmap ---
# Load the final infercnv objects we saved during the loop
completed_ids <- run_log_df$sample_id[run_log_df$status == "completed"]
infercnv_list <- list()
spot_df_list_cnv <- list()
for (sid in completed_ids) {
  rds_path <- file.path(CNV_RDS_DIR, paste0(sid, "_infercnv_final.rds"))
  if (file.exists(rds_path)) {
    infercnv_list[[sid]] <- tryCatch(readRDS(rds_path), error = function(e) NULL)
    spot_df_list_cnv[[sid]] <- sample_objects[[sid]]$spot_df
  }
}

if (length(infercnv_list) > 0) {
  cat("\nGenerating genome-wide CNV heatmaps ...\n")
  plot_genome_cnv_heatmap(
    infercnv_list = infercnv_list,
    spot_df_list = spot_df_list_cnv,
    gene_order_df = gene_order_df,
    file_stub = "02g_spatial_cnv_genome_heatmap_celltype",
    group_col = "dominant_cell_type",
    plot_width_px = 7000L, plot_height_px = 2200L, plot_dpi = 300L
  )

  gse238264_infer <- infercnv_list[names(infercnv_list) %in% names(GSE238264_SAMPLES)]
  gse238264_spots_cnv <- spot_df_list_cnv[names(spot_df_list_cnv) %in% names(GSE238264_SAMPLES)]
  if (length(gse238264_infer) > 0) {
    plot_genome_cnv_heatmap(
      infercnv_list = gse238264_infer,
      spot_df_list = gse238264_spots_cnv,
      gene_order_df = gene_order_df,
      file_stub = "02g_spatial_cnv_genome_heatmap_response",
      group_col = "group",
      plot_width_px = 7000L, plot_height_px = 2000L, plot_dpi = 300L
    )
  }
}

cnv_map_samples <- unique(spot_burden_tbl$sample_id)[seq_len(min(4, n_distinct(spot_burden_tbl$sample_id)))]
cnv_map_plots <- lapply(cnv_map_samples, function(sample_id) {
  plot_spatial_numeric(
    spot_burden_tbl %>% filter(sample_id == .env$sample_id),
    feature = "cnv_burden",
    title = sample_id,
    subtitle = "InferCNV spatial CNA burden",
    palette = "magma"
  )
})
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_spatial_maps.png"), wrap_plots(cnv_map_plots, ncol = 2), width = 11, height = 9, dpi = 500)

burden_plot <- ggplot(spot_burden_tbl, aes(x = dominant_cell_type, y = cnv_burden, fill = dominant_cell_type)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = CELLTYPE_PALETTE, drop = FALSE) +
  labs(title = "Spatial CNV burden by dominant cell type", x = NULL, y = "Mean absolute CNV deviation", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(CNV_FIG_DIR, "02g_spatial_cnv_burden_by_celltype.png"), burden_plot, width = 12, height = 6, dpi = 500)

cat("\nSpatial CNV reanalysis complete. Outputs written to:\n")
cat("  Figures:", CNV_FIG_DIR, "\n")
cat("  Tables:", CNV_TAB_DIR, "\n")