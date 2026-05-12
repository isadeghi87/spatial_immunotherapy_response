#!/usr/bin/env Rscript

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_path <- normalizePath(sub("^--file=", "", .script_arg[[1]]), mustWork = FALSE)
options(spata2_script_dir = file.path(dirname(.script_path), "spata2"))
source(file.path(getOption("spata2_script_dir"), "utils.R"))

set.seed(SPATA_SEED)
ensure_dirs(FIG_DIR, TAB_DIR, QC_DIR, QC_FIG_DIR, QC_TAB_DIR, OBJ_DIR, NOTE_DIR)

cat("\n============================================================\n")
cat("SPATA2 spatial transcriptomics reanalysis\n")
cat("============================================================\n")

manifest <- make_input_manifest()
write_input_manifest(manifest)
stop_on_missing_inputs(manifest)

sample_rows <- manifest %>% filter(status == "present")
sample_objects <- vector("list", nrow(sample_rows))
names(sample_objects) <- sample_rows$sample_id

per_sample_summary <- list()
signature_coverage <- list()
gse238264_moran <- list()
gse238264_cooccurrence <- list()
gse238264_neighbor_enrichment <- list()
gse238264_entropy <- list()
gse238264_signature_cor <- list()

for (i in seq_len(nrow(sample_rows))) {
  sample_id <- sample_rows$sample_id[[i]]
  dataset <- sample_rows$dataset[[i]]
  group_label <- sample_rows$group[[i]]
  path <- sample_rows$expected_path[[i]]

  cat(sprintf("  [%02d/%02d] %s (%s)\n", i, nrow(sample_rows), sample_id, dataset))
  sample_obj <- load_spatial_sample(path = path, sample_id = sample_id, dataset = dataset, group_label = group_label)
  sample_objects[[sample_id]] <- sample_obj
  saveRDS(sample_obj$spata, file.path(OBJ_DIR, paste0(sample_id, "_spata2.rds")))

  per_sample_summary[[sample_id]] <- summarise_sample(sample_obj)
  signature_coverage[[sample_id]] <- bind_rows(
    make_signature_coverage_table(sample_obj, IO_SIGNATURES, family = "io_signature"),
    make_signature_coverage_table(sample_obj, CELL_TYPE_MARKERS, family = "cell_type_marker")
  )

  sample_df <- sample_obj$spot_df
  qc_panel <-
    (plot_spatial_numeric(sample_df, "nCount_Spatial", title = paste(sample_id, "library size"), palette = "magma") |
       plot_spatial_numeric(sample_df, "nFeature_Spatial", title = paste(sample_id, "detected genes"), palette = "plasma") |
       plot_spatial_numeric(sample_df, "pct_mt", title = paste(sample_id, "% mitochondrial"), palette = "inferno")) /
    (plot_spatial_discrete(sample_df, "dominant_cell_type", palette = CELLTYPE_PALETTE, title = paste(sample_id, "cell type")) |
       plot_spatial_discrete(sample_df, "spatial_immunotype", palette = IMMUNOTYPE_PALETTE, title = paste(sample_id, "immunotype")) |
       plot_spatial_numeric(sample_df, "immune_gradient", title = paste(sample_id, "immune gradient"), palette = "viridis"))
  safe_save(file.path(QC_FIG_DIR, paste0(sample_id, "_qc_panel.png")), qc_panel, width = 14, height = 9, dpi = 450)

  if (dataset == "GSE238264") {
    nn_idx <- compute_knn_graph(sample_obj$coords)
    sig_cols <- intersect(names(IO_SIGNATURES), colnames(sample_df))
    gse238264_moran[[sample_id]] <- bind_rows(lapply(sig_cols, function(sig) {
      stat <- moran_stat(sample_df[[sig]], nn_idx)
      tibble(sample_id = sample_id, signature = sig, moran_I = stat[["moran_I"]], p_perm = stat[["p_perm"]])
    }))
    gse238264_cooccurrence[[sample_id]] <- cooccurrence_table(as.character(sample_df$dominant_cell_type), nn_idx) %>% mutate(sample_id = sample_id)
    gse238264_neighbor_enrichment[[sample_id]] <- neighbor_enrichment_table(as.character(sample_df$dominant_cell_type), nn_idx) %>% mutate(sample_id = sample_id)
    gse238264_entropy[[sample_id]] <- tibble(
      sample_id = sample_id,
      entropy_cell_type = shannon_entropy(table(sample_df$dominant_cell_type)),
      entropy_immunotype = shannon_entropy(table(sample_df$spatial_immunotype))
    )
    gse238264_signature_cor[[sample_id]] <- bind_rows(lapply(sig_cols, function(sig) {
      cor_obj <- suppressWarnings(stats::cor.test(sample_df[[sig]], sample_df$immune_gradient, method = "spearman"))
      tibble(sample_id = sample_id, signature = sig, spearman_rho = unname(cor_obj$estimate), p_value = cor_obj$p.value)
    }))
  }
}

all_spots <- bind_rows(lapply(sample_objects, `[[`, "spot_df"))
sample_stats <- bind_rows(per_sample_summary)
spot_catalog <- make_spot_qc_catalog(sample_objects)
qc_summary <- summarise_spot_qc_catalog(spot_catalog)
annotation_catalog <- all_spots %>%
  select(sample_id, dataset, group, barcodes, dominant_cell_type, dominant_score, second_score, margin_score, annotation_confidence)
coverage_tbl <- bind_rows(signature_coverage)

write_csv(spot_catalog, file.path(QC_TAB_DIR, "spata2_spot_qc_catalog.csv"))
write_csv(qc_summary, file.path(QC_TAB_DIR, "spata2_spot_qc_summary.csv"))
write_csv(coverage_tbl, file.path(TAB_DIR, "spata2_signature_coverage.csv"))
write_csv(annotation_catalog, file.path(TAB_DIR, "spata2_celltype_annotation_catalog.csv"))
write_csv(sample_stats, file.path(TAB_DIR, "spata2_per_sample_summary.csv"))

gse238264_spots <- all_spots %>% filter(dataset == "GSE238264")
gse281759_spots <- all_spots %>% filter(dataset == "GSE281759")
write_csv(gse238264_spots, file.path(TAB_DIR, "gse238264_celltype_annotations.csv"))
write_csv(gse281759_spots, file.path(TAB_DIR, "gse281759_celltype_annotations.csv"))

gse238264_per_sample <- sample_stats %>%
  filter(dataset == "GSE238264") %>%
  transmute(
    sample,
    response = group,
    n_spots,
    across(starts_with("mean_")),
    across(starts_with("sd_")),
    across(starts_with("frac_")),
    median_nCount,
    median_nFeature,
    median_pct_mt
  )
write_csv(gse238264_per_sample, file.path(TAB_DIR, "gse238264_per_sample_stats.csv"))

gse238264_comp <- gse238264_spots %>%
  count(sample_id, dominant_cell_type, name = "n_spots") %>%
  group_by(sample_id) %>%
  mutate(frac = n_spots / sum(n_spots)) %>%
  ungroup()
write_csv(gse238264_comp, file.path(TAB_DIR, "gse238264_celltype_composition.csv"))

gse238264_immunotype_comp <- gse238264_spots %>%
  count(sample_id, group, spatial_immunotype, name = "n_spots") %>%
  group_by(sample_id) %>%
  mutate(frac = n_spots / sum(n_spots)) %>%
  ungroup()
write_csv(gse238264_immunotype_comp, file.path(TAB_DIR, "gse238264_immunotype_composition.csv"))

sig_cols <- intersect(names(IO_SIGNATURES), colnames(gse238264_spots))
gse238264_sig_stats <- adjust_stat_table(pairwise_group_stats(gse238264_spots, "group", c("Responder", "Non_Responder"), sig_cols)) %>%
  arrange(p_adj, p_value)
write_csv(gse238264_sig_stats, file.path(TAB_DIR, "gse238264_RvsNR_signature_stats.csv"))

gse238264_sample_sig_stats <- adjust_stat_table(pairwise_group_stats(gse238264_per_sample, "response", c("Responder", "Non_Responder"), paste0("mean_", sig_cols))) %>%
  arrange(p_adj, p_value)
write_csv(gse238264_sample_sig_stats, file.path(TAB_DIR, "gse238264_RvsNR_sample_level_signature_stats.csv"))

gse238264_perm <- sample_level_permutation(gse238264_per_sample, "response", paste0("mean_", sig_cols))
write_csv(gse238264_perm, file.path(TAB_DIR, "gse238264_sample_level_permutation.csv"))

gse238264_loo <- leave_one_out_stats(gse238264_per_sample, "response", paste0("mean_", sig_cols))
write_csv(gse238264_loo, file.path(TAB_DIR, "gse238264_leave_one_out.csv"))

gse238264_gradient <- gse238264_spots %>%
  group_by(sample_id, group) %>%
  summarise(
    mean_immune_gradient = mean(immune_gradient, na.rm = TRUE),
    mean_tls = mean(TLS_signature, na.rm = TRUE),
    mean_activation = mean(T_cell_activation, na.rm = TRUE),
    mean_tgfb = mean(TGF_beta, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(gse238264_gradient, file.path(TAB_DIR, "gse238264_immune_gradient.csv"))

write_csv(bind_rows(gse238264_moran), file.path(TAB_DIR, "gse238264_morans_I.csv"))
write_csv(bind_rows(gse238264_cooccurrence), file.path(TAB_DIR, "gse238264_ligand_receptor.csv"))
write_csv(bind_rows(gse238264_neighbor_enrichment), file.path(TAB_DIR, "gse238264_neighbor_enrichment.csv"))
write_csv(adjust_stat_table(bind_rows(gse238264_signature_cor)), file.path(TAB_DIR, "gse238264_signature_gradient_correlations.csv"))
write_csv(bind_rows(gse238264_entropy), file.path(TAB_DIR, "gse238264_spatial_entropy.csv"))

if (length(sig_cols) >= 2) {
  gse238264_cluster_profiles <- gse238264_spots %>%
    mutate(spatial_cluster = kmeans(scale(select(., all_of(sig_cols))), centers = min(4, n_distinct(sample_id)), iter.max = 50)$cluster) %>%
    group_by(sample_id, spatial_cluster) %>%
    summarise(across(all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  write_csv(gse238264_cluster_profiles, file.path(TAB_DIR, "gse238264_spatial_cluster_profiles.csv"))
}

gse281759_sig_cols <- intersect(names(IO_SIGNATURES), colnames(gse281759_spots))
gse281759_tumor_normal <- adjust_stat_table(pairwise_group_stats(gse281759_spots %>% filter(group %in% c("Tumor", "Normal")), "group", c("Tumor", "Normal"), gse281759_sig_cols))
write_csv(gse281759_tumor_normal, file.path(TAB_DIR, "gse281759_tumor_vs_normal_stats.csv"))

gse281759_primary_recurrence <- adjust_stat_table(pairwise_group_stats(gse281759_spots %>% filter(group %in% c("Primary", "Recurrence")), "group", c("Primary", "Recurrence"), gse281759_sig_cols))
write_csv(gse281759_primary_recurrence, file.path(TAB_DIR, "gse281759_primary_vs_recurrence_stats.csv"))

# ---- QC figures ------------------------------------------------------------

qc_bar <- ggplot(qc_summary, aes(x = reorder(sample_id, n_spots), y = n_spots, fill = dataset)) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(title = "Spots per section", x = NULL, y = "Number of spots", fill = NULL) +
  theme_pub()

qc_umi <- ggplot(spot_catalog, aes(x = sample_id, y = log10(nCount_Spatial + 1), fill = dataset)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  labs(title = "Library size per spot", x = NULL, y = "log10(UMI + 1)", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

qc_gene <- ggplot(spot_catalog, aes(x = sample_id, y = nFeature_Spatial, fill = dataset)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  labs(title = "Genes detected per spot", x = NULL, y = "Detected genes", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

qc_conf <- annotation_catalog %>%
  count(sample_id, annotation_confidence, name = "n") %>%
  group_by(sample_id) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = sample_id, y = frac, fill = annotation_confidence)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(Low = "#D55E00", Medium = "#E69F00", High = "#009E73"), drop = FALSE) +
  labs(title = "Cell-type annotation confidence", x = NULL, y = "Fraction of spots", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

safe_save(file.path(QC_FIG_DIR, "02b_qc_overview.png"), (qc_bar | qc_conf) / (qc_umi | qc_gene), width = 18, height = 12, dpi = 450)

# ---- Figures ---------------------------------------------------------------

celltype_map_samples <- c(
  head(names(GSE238264_SAMPLES)[GSE238264_SAMPLES == "Responder"], 2),
  head(names(GSE238264_SAMPLES)[GSE238264_SAMPLES == "Non_Responder"], 2)
)
celltype_map_plots <- lapply(celltype_map_samples, function(sample_id) {
  sample_df <- sample_objects[[sample_id]]$spot_df
  plot_spatial_discrete(sample_df, "dominant_cell_type", CELLTYPE_PALETTE, title = sample_id, subtitle = sample_objects[[sample_id]]$group)
})
safe_save(file.path(FIG_DIR, "02b_ct_spatial_celltype_maps.png"), wrap_plots(celltype_map_plots, ncol = 2), width = 12, height = 10, dpi = 500)

immunotype_map_plots <- lapply(celltype_map_samples, function(sample_id) {
  sample_df <- sample_objects[[sample_id]]$spot_df
  plot_spatial_discrete(sample_df, "spatial_immunotype", IMMUNOTYPE_PALETTE, title = sample_id, subtitle = "Spatial immunotype")
})
safe_save(file.path(FIG_DIR, "02b_visium_spatial_immunotype_maps.png"), wrap_plots(immunotype_map_plots, ncol = 2), width = 12, height = 10, dpi = 500)

umap_df <- gse238264_spots %>% mutate(response = factor(group, levels = c("Responder", "Non_Responder")))
umap_panel <- plot_umap_by_feature(umap_df, "dominant_cell_type", discrete = TRUE, palette = CELLTYPE_PALETTE, title = "Cell type") |
  plot_umap_by_feature(umap_df, "response", discrete = TRUE, palette = RESPONSE_PALETTE, title = "Response") |
  plot_umap_by_feature(umap_df, "spatial_immunotype", discrete = TRUE, palette = IMMUNOTYPE_PALETTE, title = "Immunotype")
safe_save(file.path(FIG_DIR, "02b_ct_umap_celltype_response.png"), umap_panel, width = 16, height = 5, dpi = 500)

marker_cols <- grep("^ct_", colnames(gse238264_spots), value = TRUE)
marker_heatmap_mat <- gse238264_spots %>%
  group_by(dominant_cell_type) %>%
  summarise(across(all_of(marker_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()
rownames(marker_heatmap_mat) <- marker_heatmap_mat$dominant_cell_type
marker_heatmap_mat$dominant_cell_type <- NULL
png(file.path(FIG_DIR, "02b_ct_marker_heatmap.png"), width = 2200, height = 1400, res = 240)
pheatmap(as.matrix(marker_heatmap_mat), border_color = NA, main = "Cell-type marker score heatmap")
dev.off()

violin_df <- gse238264_spots %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse238264_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
violin_plot <- ggplot(violin_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = RESPONSE_PALETTE) +
  labs(title = "GSE238264: responder vs non-responder signatures", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_RvsNR_violins.png"), violin_plot, width = 12, height = 6, dpi = 500)

sample_heatmap_mat <- gse238264_per_sample %>%
  select(sample, starts_with("mean_")) %>%
  as.data.frame()
rownames(sample_heatmap_mat) <- sample_heatmap_mat$sample
sample_heatmap_mat$sample <- NULL
sample_annotation <- data.frame(response = gse238264_per_sample$response, row.names = gse238264_per_sample$sample)
png(file.path(FIG_DIR, "02b_visium_sample_signature_heatmap.png"), width = 2200, height = 1200, res = 240)
pheatmap(as.matrix(sample_heatmap_mat), border_color = NA, annotation_row = sample_annotation, main = "GSE238264 sample-level signature means")
dev.off()

sig_map_samples <- celltype_map_samples
sig_candidates <- KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse238264_spots)]
sig_map_names <- if (length(sig_candidates) > 0) sig_candidates[seq_len(min(4, length(sig_candidates)))] else character(0)
sig_map_plots <- unlist(lapply(sig_map_samples, function(sample_id) {
  sample_df <- sample_objects[[sample_id]]$spot_df
  lapply(sig_map_names, function(sig) {
    plot_spatial_numeric(sample_df, sig, title = sample_id, subtitle = sig, palette = "D")
  })
}), recursive = FALSE)
if (length(sig_map_plots) > 0) {
  safe_save(file.path(FIG_DIR, "02b_visium_spatial_signature_maps.png"), wrap_plots(sig_map_plots, ncol = length(sig_map_names)), width = 18, height = 10, dpi = 500)
}

immunotype_plot <- ggplot(gse238264_immunotype_comp, aes(x = sample_id, y = frac, fill = spatial_immunotype)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = IMMUNOTYPE_PALETTE, drop = FALSE) +
  labs(title = "Spatial immunotype composition", x = NULL, y = "Fraction of spots", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_immunotype_composition.png"), immunotype_plot, width = 10, height = 5, dpi = 500)

moran_plot_df <- bind_rows(gse238264_moran)
moran_plot <- ggplot(moran_plot_df, aes(x = signature, y = moran_I, fill = sample_id)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "GSE238264 spatial autocorrelation", x = NULL, y = "Moran's I", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_morans_I.png"), moran_plot, width = 12, height = 5, dpi = 500)

neighbor_heatmap <- bind_rows(gse238264_neighbor_enrichment) %>%
  group_by(source, target) %>%
  summarise(mean_log2_enrichment = mean(log2_enrichment, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = target, values_from = mean_log2_enrichment, values_fill = 0) %>%
  as.data.frame()
rownames(neighbor_heatmap) <- neighbor_heatmap$source
neighbor_heatmap$source <- NULL
png(file.path(FIG_DIR, "02b_visium_neighbor_enrichment_heatmap.png"), width = 1800, height = 1400, res = 240)
pheatmap(as.matrix(neighbor_heatmap), border_color = NA, main = "Dominant cell-type neighbor enrichment")
dev.off()

tn_plot_df <- gse281759_spots %>%
  filter(group %in% c("Tumor", "Normal")) %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse281759_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
tn_plot <- ggplot(tn_plot_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = TISSUE_PALETTE[c("Tumor", "Normal")]) +
  labs(title = "GSE281759: tumor vs normal", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_gse281759_tumor_vs_normal.png"), tn_plot, width = 12, height = 6, dpi = 500)

pr_plot_df <- gse281759_spots %>%
  filter(group %in% c("Primary", "Recurrence")) %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse281759_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
pr_plot <- ggplot(pr_plot_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = TISSUE_PALETTE[c("Primary", "Recurrence")]) +
  labs(title = "GSE281759: primary vs recurrence", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_gse281759_primary_vs_recurrence.png"), pr_plot, width = 12, height = 6, dpi = 500)

write_note_lines(
  file.path(NOTE_DIR, "02b_spata2_reanalysis_note.txt"),
  c(
    "Question: Reanalyse the Visium spatial cohorts with SPATA2-based QC, annotation, and response/topology summaries.",
    "Comparison: GSE238264 responder vs non-responder; GSE281759 tumor vs normal and primary vs recurrence.",
    sprintf("Inputs: %d Visium sections from %s and %s.", nrow(sample_rows), "GSE238264", "GSE281759"),
    "Key outputs: QC catalog, annotation catalog, signature coverage table, response statistics, neighborhood enrichment, and publication-grade summary figures.",
    "QC status: input manifest generated; per-section QC panels regenerated; sample-level QC summary exported."
  )
)

cat("\nSpatial SPATA2 reanalysis complete. Outputs written to:\n")
cat("  Figures:", FIG_DIR, "\n")
cat("  QC:", QC_DIR, "\n")
cat("  Tables:", TAB_DIR, "\n")

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.script_path <- normalizePath(sub("^--file=", "", .script_arg[[1]]), mustWork = FALSE)
options(spata2_script_dir = file.path(dirname(.script_path), "spata2"))
source(file.path(getOption("spata2_script_dir"), "utils.R"))

set.seed(SPATA_SEED)
ensure_dirs(FIG_DIR, TAB_DIR, QC_DIR, OBJ_DIR, NOTE_DIR)

cat("\n============================================================\n")
cat("SPATA2 spatial transcriptomics reanalysis\n")
cat("============================================================\n")

manifest <- make_input_manifest()
write_input_manifest(manifest)
stop_on_missing_inputs(manifest)

sample_rows <- manifest %>% filter(status == "present")
sample_objects <- vector("list", nrow(sample_rows))
names(sample_objects) <- sample_rows$sample_id

per_sample_summary <- list()
gse238264_moran <- list()
gse238264_cooccurrence <- list()
gse238264_entropy <- list()

for (i in seq_len(nrow(sample_rows))) {
  sample_id <- sample_rows$sample_id[[i]]
  dataset <- sample_rows$dataset[[i]]
  group_label <- sample_rows$group[[i]]
  path <- sample_rows$expected_path[[i]]

  cat(sprintf("  [%02d/%02d] %s (%s)\n", i, nrow(sample_rows), sample_id, dataset))
  sample_obj <- load_spatial_sample(path = path, sample_id = sample_id, dataset = dataset, group_label = group_label)
  sample_objects[[sample_id]] <- sample_obj
  saveRDS(sample_obj$spata, file.path(OBJ_DIR, paste0(sample_id, "_spata2.rds")))

  per_sample_summary[[sample_id]] <- summarise_sample(sample_obj)

  sample_df <- sample_obj$spot_df
  qc_panel <- (
    plot_spatial_numeric(sample_df, "nCount_Spatial", title = paste(sample_id, "library size"), palette = "magma") +
      plot_spatial_numeric(sample_df, "nFeature_Spatial", title = paste(sample_id, "detected genes"), palette = "plasma")
  ) / (
    plot_spatial_numeric(sample_df, "pct_mt", title = paste(sample_id, "% mitochondrial"), palette = "inferno") +
      plot_umap_by_feature(sample_df, "dominant_cell_type", discrete = TRUE, palette = CELLTYPE_PALETTE, title = paste(sample_id, "UMAP"))
  )
  safe_save(file.path(QC_DIR, paste0(sample_id, "_qc_panel.png")), qc_panel, width = 12, height = 10, dpi = 450)

  if (dataset == "GSE238264") {
    nn_idx <- compute_knn_graph(sample_obj$coords)
    moran_tbl <- bind_rows(lapply(intersect(names(IO_SIGNATURES), colnames(sample_df)), function(sig) {
      stat <- moran_stat(sample_df[[sig]], nn_idx)
      tibble(sample = sample_id, signature = sig, moran_I = stat[["moran_I"]], p_perm = stat[["p_perm"]])
    }))
    gse238264_moran[[sample_id]] <- moran_tbl
    gse238264_cooccurrence[[sample_id]] <- cooccurrence_table(as.character(sample_df$dominant_cell_type), nn_idx) %>% mutate(sample = sample_id)
    gse238264_entropy[[sample_id]] <- tibble(
      sample = sample_id,
      entropy_cell_type = shannon_entropy(table(sample_df$dominant_cell_type)),
      entropy_immunotype = shannon_entropy(table(sample_df$spatial_immunotype))
    )
  }
}

all_spots <- bind_rows(lapply(sample_objects, `[[`, "spot_df"))
sample_stats <- bind_rows(per_sample_summary)

gse238264_spots <- all_spots %>% filter(dataset == "GSE238264")
gse281759_spots <- all_spots %>% filter(dataset == "GSE281759")

write_csv(sample_stats, file.path(TAB_DIR, "spata2_per_sample_summary.csv"))
write_csv(gse238264_spots, file.path(TAB_DIR, "gse238264_celltype_annotations.csv"))
write_csv(gse281759_spots, file.path(TAB_DIR, "gse281759_celltype_annotations.csv"))

gse238264_per_sample <- sample_stats %>%
  filter(dataset == "GSE238264") %>%
  transmute(
    sample,
    response = group,
    n_spots,
    across(starts_with("mean_")),
    across(starts_with("sd_")),
    across(starts_with("frac_")),
    median_nCount,
    median_nFeature,
    median_pct_mt
  )
write_csv(gse238264_per_sample, file.path(TAB_DIR, "gse238264_per_sample_stats.csv"))

gse238264_comp <- gse238264_spots %>%
  count(sample_id, dominant_cell_type, name = "n_spots") %>%
  group_by(sample_id) %>%
  mutate(frac = n_spots / sum(n_spots)) %>%
  ungroup()
write_csv(gse238264_comp, file.path(TAB_DIR, "gse238264_celltype_composition.csv"))

sig_cols <- intersect(names(IO_SIGNATURES), colnames(gse238264_spots))
gse238264_sig_stats <- pairwise_group_stats(gse238264_spots, "group", c("Responder", "Non_Responder"), sig_cols) %>%
  arrange(p_value)
write_csv(gse238264_sig_stats, file.path(TAB_DIR, "gse238264_RvsNR_signature_stats.csv"))

gse238264_perm <- sample_level_permutation(gse238264_per_sample, "response", paste0("mean_", sig_cols))
write_csv(gse238264_perm, file.path(TAB_DIR, "gse238264_sample_level_permutation.csv"))

gse238264_loo <- leave_one_out_stats(gse238264_per_sample, "response", paste0("mean_", sig_cols))
write_csv(gse238264_loo, file.path(TAB_DIR, "gse238264_leave_one_out.csv"))

gse238264_gradient <- gse238264_spots %>%
  group_by(sample_id, group) %>%
  summarise(
    mean_immune_gradient = mean(immune_gradient, na.rm = TRUE),
    mean_tls = mean(TLS_signature, na.rm = TRUE),
    mean_activation = mean(T_cell_activation, na.rm = TRUE),
    mean_tgfb = mean(TGF_beta, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(gse238264_gradient, file.path(TAB_DIR, "gse238264_immune_gradient.csv"))

write_csv(bind_rows(gse238264_moran), file.path(TAB_DIR, "gse238264_morans_I.csv"))
write_csv(bind_rows(gse238264_cooccurrence), file.path(TAB_DIR, "gse238264_ligand_receptor.csv"))
write_csv(bind_rows(gse238264_entropy), file.path(TAB_DIR, "gse238264_spatial_entropy.csv"))

gse238264_cluster_profiles <- gse238264_spots %>%
  mutate(spatial_cluster = kmeans(scale(select(., all_of(sig_cols))), centers = min(4, n_distinct(sample_id)), iter.max = 50)$cluster) %>%
  group_by(sample_id, spatial_cluster) %>%
  summarise(across(all_of(sig_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
write_csv(gse238264_cluster_profiles, file.path(TAB_DIR, "gse238264_spatial_cluster_profiles.csv"))

gse281759_tumor_normal <- pairwise_group_stats(gse281759_spots %>% filter(group %in% c("Tumor", "Normal")), "group", c("Tumor", "Normal"), sig_cols)
write_csv(gse281759_tumor_normal, file.path(TAB_DIR, "gse281759_tumor_vs_normal_stats.csv"))

gse281759_primary_recurrence <- pairwise_group_stats(gse281759_spots %>% filter(group %in% c("Primary", "Recurrence")), "group", c("Primary", "Recurrence"), sig_cols)
write_csv(gse281759_primary_recurrence, file.path(TAB_DIR, "gse281759_primary_vs_recurrence_stats.csv"))

# ---- Figures ---------------------------------------------------------------

celltype_map_samples <- names(sample_objects)[seq_len(min(4, length(sample_objects)))]
celltype_map_plots <- lapply(celltype_map_samples, function(sample_id) {
  sample_df <- sample_objects[[sample_id]]$spot_df
  plot_spatial_discrete(sample_df, "dominant_cell_type", CELLTYPE_PALETTE,
                        title = sample_id,
                        subtitle = sample_objects[[sample_id]]$group)
})
safe_save(file.path(FIG_DIR, "02b_ct_spatial_celltype_maps.png"), wrap_plots(celltype_map_plots, ncol = 2), width = 12, height = 10, dpi = 450)

umap_df <- gse238264_spots %>% mutate(response = factor(group, levels = c("Responder", "Non_Responder")))
umap_panel <- plot_umap_by_feature(umap_df, "dominant_cell_type", discrete = TRUE, palette = CELLTYPE_PALETTE, title = "Cell type") +
  plot_umap_by_feature(umap_df, "response", discrete = TRUE, palette = RESPONSE_PALETTE, title = "Response")
safe_save(file.path(FIG_DIR, "02b_ct_umap_celltype_response.png"), umap_panel, width = 12, height = 5, dpi = 450)

marker_heatmap_mat <- gse238264_spots %>%
  group_by(dominant_cell_type) %>%
  summarise(across(starts_with("ct_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  as.data.frame()
rownames(marker_heatmap_mat) <- marker_heatmap_mat$dominant_cell_type
marker_heatmap_mat$dominant_cell_type <- NULL
png(file.path(FIG_DIR, "02b_ct_marker_heatmap.png"), width = 1800, height = 1200, res = 220)
pheatmap(as.matrix(marker_heatmap_mat), border_color = NA, main = "Cell-type module score heatmap")
dev.off()

violin_df <- gse238264_spots %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse238264_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
violin_plot <- ggplot(violin_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = RESPONSE_PALETTE) +
  labs(title = "GSE238264: R vs NR signature distributions", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_RvsNR_violins.png"), violin_plot, width = 12, height = 6, dpi = 450)

sample_heatmap_mat <- gse238264_per_sample %>%
  select(sample, starts_with("mean_")) %>%
  as.data.frame()
rownames(sample_heatmap_mat) <- sample_heatmap_mat$sample
sample_heatmap_mat$sample <- NULL
png(file.path(FIG_DIR, "02b_visium_sample_signature_heatmap.png"), width = 1800, height = 1200, res = 220)
pheatmap(as.matrix(sample_heatmap_mat), border_color = NA, main = "GSE238264 sample-level signature means")
dev.off()

sig_map_samples <- names(sample_objects)[names(sample_objects) %in% names(GSE238264_SAMPLES)][seq_len(min(4, sum(names(sample_objects) %in% names(GSE238264_SAMPLES))))]
sig_map_plots <- unlist(lapply(sig_map_samples, function(sample_id) {
  sample_df <- sample_objects[[sample_id]]$spot_df
  lapply(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(sample_df)][1:min(4, length(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(sample_df)]))], function(sig) {
    plot_spatial_numeric(sample_df, sig, title = sample_id, subtitle = sig, palette = "D")
  })
}), recursive = FALSE)
if (length(sig_map_plots) > 0) {
  safe_save(file.path(FIG_DIR, "02b_visium_spatial_signature_maps.png"), wrap_plots(sig_map_plots, ncol = 4), width = 16, height = 10, dpi = 450)
}

immunotype_df <- gse238264_spots %>%
  count(sample_id, group, spatial_immunotype, name = "n") %>%
  group_by(sample_id) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
immunotype_plot <- ggplot(immunotype_df, aes(x = sample_id, y = frac, fill = spatial_immunotype)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = IMMUNOTYPE_PALETTE, drop = FALSE) +
  labs(title = "Spatial immunotype composition", x = NULL, y = "Fraction of spots", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_immunotype_composition.png"), immunotype_plot, width = 10, height = 5, dpi = 450)

moran_plot_df <- bind_rows(gse238264_moran)
moran_plot <- ggplot(moran_plot_df, aes(x = signature, y = moran_I, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "GSE238264 spatial autocorrelation", x = NULL, y = "Moran's I", fill = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_morans_I.png"), moran_plot, width = 12, height = 5, dpi = 450)

tn_plot_df <- gse281759_spots %>%
  filter(group %in% c("Tumor", "Normal")) %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse281759_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
tn_plot <- ggplot(tn_plot_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = TISSUE_PALETTE[c("Tumor", "Normal")]) +
  labs(title = "GSE281759: Tumor vs Normal", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_gse281759_tumor_vs_normal.png"), tn_plot, width = 12, height = 6, dpi = 450)

pr_plot_df <- gse281759_spots %>%
  filter(group %in% c("Primary", "Recurrence")) %>%
  select(sample_id, group, all_of(KEY_SIGNATURES[KEY_SIGNATURES %in% colnames(gse281759_spots)])) %>%
  pivot_longer(cols = -c(sample_id, group), names_to = "signature", values_to = "score")
pr_plot <- ggplot(pr_plot_df, aes(x = signature, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  scale_fill_manual(values = TISSUE_PALETTE[c("Primary", "Recurrence")]) +
  labs(title = "GSE281759: Primary vs Recurrence", x = NULL, y = "Module score", fill = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
safe_save(file.path(FIG_DIR, "02b_visium_gse281759_primary_vs_recurrence.png"), pr_plot, width = 12, height = 6, dpi = 450)

cat("\nSpatial SPATA2 reanalysis complete. Outputs written to:\n")
cat("  Figures:", FIG_DIR, "\n")
cat("  Tables:", TAB_DIR, "\n")