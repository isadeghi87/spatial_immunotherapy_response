#!/usr/bin/env Rscript

.spata2_dir <- getOption("spata2_script_dir", default = getwd())
source(file.path(.spata2_dir, "config.R"))

suppressPackageStartupMessages({
  library(SPATA2)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(readr)
  library(stringr)
  library(purrr)
  library(scales)
  library(RANN)
  library(uwot)
  library(infercnv)
})

ensure_dirs <- function(...) {
  for (dir_path in c(...)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(colour = "grey35", hjust = 0, margin = margin(b = 6)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = rel(0.85)),
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text = element_text(size = rel(0.8)),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

theme_spatial <- function(base_size = 10) {
  theme_void(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(colour = "grey35", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text = element_text(size = rel(0.8))
    )
}

safe_save <- function(filename, plot, width = 9, height = 7, dpi = 400, bg = "white") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = bg)
}

adaptive_point_size <- function(n_obs, min_size = 0.25, max_size = 1.6) {
  if (!is.finite(n_obs) || n_obs <= 0) return(max_size)
  size <- 24 / sqrt(n_obs)
  max(min_size, min(max_size, size))
}

clip_numeric <- function(x, lower_q = 0.01, upper_q = 0.99) {
  if (!is.numeric(x) || all(!is.finite(x))) return(x)
  lo <- stats::quantile(x, lower_q, na.rm = TRUE)
  hi <- stats::quantile(x, upper_q, na.rm = TRUE)
  pmin(pmax(x, lo), hi)
}

shannon_entropy <- function(x) {
  probs <- as.numeric(x) / sum(x)
  probs <- probs[is.finite(probs) & probs > 0]
  -sum(probs * log2(probs))
}

make_input_manifest <- function() {
  manifest <- bind_rows(
    tibble(
      dataset = "GSE238264",
      sample_id = names(GSE238264_SAMPLES),
      expected_path = file.path(GSE238264_DIR, paste0(names(GSE238264_SAMPLES), ".h5ad")),
      group = unname(GSE238264_SAMPLES)
    ),
    tibble(
      dataset = "GSE281759",
      sample_id = sub("\\.h5ad$", "", basename(list.files(GSE281759_DIR, pattern = "\\.h5ad$", full.names = FALSE))),
      expected_path = list.files(GSE281759_DIR, pattern = "\\.h5ad$", full.names = TRUE),
      group = NA_character_
    )
  ) %>%
    mutate(
      symlink_target = vapply(expected_path, function(path) Sys.readlink(path), character(1)),
      exists = file.exists(expected_path),
      status = ifelse(exists, "present", ifelse(nzchar(symlink_target), "broken_symlink", "missing"))
    )

  if (nrow(manifest) == 0) {
    manifest <- tibble(
      dataset = character(), sample_id = character(), expected_path = character(),
      group = character(), symlink_target = character(), exists = logical(), status = character()
    )
  }

  manifest$group[manifest$dataset == "GSE281759"] <- vapply(
    manifest$sample_id[manifest$dataset == "GSE281759"], infer_gse281759_group, character(1)
  )
  manifest
}

write_input_manifest <- function(manifest) {
  ensure_dirs(TAB_DIR, NOTE_DIR)
  write_csv(manifest, file.path(TAB_DIR, "spata2_input_manifest.csv"))
  if (any(manifest$status != "present")) {
    lines <- c(
      "SPATA2 reanalysis input audit",
      sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      "Missing or broken inputs:",
      apply(manifest %>% filter(status != "present"), 1, function(row) {
        sprintf("- %s [%s]: %s (target=%s)", row[["sample_id"]], row[["status"]], row[["expected_path"]], row[["symlink_target"]])
      })
    )
    writeLines(lines, file.path(NOTE_DIR, "spata2_input_audit.txt"))
  }
}

stop_on_missing_inputs <- function(manifest) {
  missing_manifest <- manifest %>% filter(status != "present")
  if (nrow(missing_manifest) == 0) return(invisible(NULL))
  stop(
    paste0(
      "SPATA2 spatial reanalysis cannot run because ", nrow(missing_manifest),
      " input file(s) are missing or broken symlinks. See ",
      file.path(TAB_DIR, "spata2_input_manifest.csv"),
      " and ", file.path(NOTE_DIR, "spata2_input_audit.txt"), "."
    ),
    call. = FALSE
  )
}

resolve_counts_assay <- function(sce) {
  assay_names <- assayNames(sce)
  if ("counts" %in% assay_names) return("counts")
  if ("X" %in% assay_names) return("X")
  assay_names[[1]]
}

extract_coords_df <- function(sce) {
  rd_names <- reducedDimNames(sce)
  if ("spatial" %in% rd_names) {
    coords <- reducedDim(sce, "spatial")
    return(data.frame(
      barcodes = colnames(sce),
      x = coords[, 1],
      y = coords[, 2],
      row.names = colnames(sce),
      check.names = FALSE
    ))
  }

  meta <- as.data.frame(colData(sce))
  meta$barcodes <- colnames(sce)
  coord_pairs <- list(
    c("pxl_col", "pxl_row"),
    c("imagecol", "imagerow"),
    c("x", "y"),
    c("spatial_1", "spatial_2")
  )
  for (pair in coord_pairs) {
    if (all(pair %in% colnames(meta))) {
      return(data.frame(
        barcodes = meta$barcodes,
        x = meta[[pair[[1]]]],
        y = meta[[pair[[2]]]],
        row.names = meta$barcodes,
        check.names = FALSE
      ))
    }
  }
  stop("Could not identify spatial coordinates in reducedDims or colData.", call. = FALSE)
}

log_normalize_counts <- function(count_mtr, scale_factor = 1e4) {
  lib_size <- Matrix::colSums(count_mtr)
  lib_size[lib_size == 0] <- 1
  norm_mtr <- t(t(count_mtr) / lib_size) * scale_factor
  norm_mtr@x <- log1p(norm_mtr@x)
  norm_mtr
}

module_score_matrix <- function(expr_mtr, gene_sets, prefix = NULL) {
  gene_names <- toupper(rownames(expr_mtr))
  rownames(expr_mtr) <- gene_names
  score_list <- lapply(names(gene_sets), function(set_name) {
    genes <- intersect(toupper(gene_sets[[set_name]]), gene_names)
    if (length(genes) < 2) return(rep(NA_real_, ncol(expr_mtr)))
    Matrix::colMeans(expr_mtr[genes, , drop = FALSE])
  })
  score_df <- as.data.frame(do.call(cbind, score_list))
  colnames(score_df) <- if (is.null(prefix)) names(gene_sets) else paste0(prefix, names(gene_sets))
  rownames(score_df) <- colnames(expr_mtr)
  score_df
}

annotate_celltype_confidence <- function(cell_score_df) {
  score_mtr <- as.matrix(cell_score_df)
  if (nrow(score_mtr) == 0 || ncol(score_mtr) == 0) {
    return(tibble(
      dominant_cell_type = factor(character(), levels = c(names(CELL_TYPE_MARKERS), "Mixed")),
      dominant_score = numeric(),
      second_score = numeric(),
      margin_score = numeric(),
      annotation_confidence = factor(character(), levels = c("Low", "Medium", "High"))
    ))
  }

  label_idx <- max.col(score_mtr, ties.method = "first")
  max_vals <- score_mtr[cbind(seq_len(nrow(score_mtr)), label_idx)]
  second_vals <- vapply(seq_len(nrow(score_mtr)), function(i) {
    vals <- sort(score_mtr[i, ], decreasing = TRUE)
    if (length(vals) < 2) return(NA_real_)
    vals[[2]]
  }, numeric(1))
  margin_vals <- max_vals - second_vals
  labels <- colnames(score_mtr)[label_idx]

  low_signal_cut <- stats::quantile(max_vals, 0.25, na.rm = TRUE)
  low_margin_cut <- stats::quantile(margin_vals, 0.25, na.rm = TRUE)
  low_signal <- is.na(max_vals) | max_vals < low_signal_cut
  low_margin <- is.na(margin_vals) | margin_vals < low_margin_cut
  labels[low_signal & low_margin] <- "Mixed"

  high_margin_cut <- stats::quantile(margin_vals, 0.75, na.rm = TRUE)
  mid_margin_cut <- stats::quantile(margin_vals, 0.40, na.rm = TRUE)
  confidence <- dplyr::case_when(
    labels == "Mixed" ~ "Low",
    margin_vals >= high_margin_cut ~ "High",
    margin_vals >= mid_margin_cut ~ "Medium",
    TRUE ~ "Low"
  )

  tibble(
    dominant_cell_type = factor(labels, levels = c(names(CELL_TYPE_MARKERS), "Mixed")),
    dominant_score = as.numeric(max_vals),
    second_score = as.numeric(second_vals),
    margin_score = as.numeric(margin_vals),
    annotation_confidence = factor(confidence, levels = c("Low", "Medium", "High"))
  )
}

annotate_dominant_celltype <- function(cell_score_df) {
  annotate_celltype_confidence(cell_score_df)$dominant_cell_type
}

classify_immunotype <- function(signature_df) {
  if (!all(c("T_cell_activation", "TGF_beta") %in% colnames(signature_df))) {
    return(factor(rep(NA_character_, nrow(signature_df)), levels = IMMUNOTYPE_LEVELS))
  }
  act <- signature_df$T_cell_activation
  tgf <- signature_df$TGF_beta
  p66_act <- stats::quantile(act, 0.66, na.rm = TRUE)
  p33_act <- stats::quantile(act, 0.33, na.rm = TRUE)
  p66_tgf <- stats::quantile(tgf, 0.66, na.rm = TRUE)
  labels <- ifelse(act > p66_act & tgf < p66_tgf, "Immune_Hot",
                   ifelse(act > p66_act & tgf >= p66_tgf, "Immune_Excluded",
                          ifelse(act < p33_act, "Immune_Desert", "Immune_Mixed")))
  factor(labels, levels = IMMUNOTYPE_LEVELS)
}

compute_qc_metrics <- function(count_mtr) {
  gene_names <- rownames(count_mtr)
  mt_mask <- grepl("^MT-", gene_names, ignore.case = TRUE)
  n_count <- Matrix::colSums(count_mtr)
  n_feature <- Matrix::colSums(count_mtr > 0)
  pct_mt <- if (any(mt_mask)) Matrix::colSums(count_mtr[mt_mask, , drop = FALSE]) / pmax(n_count, 1) * 100 else rep(0, ncol(count_mtr))
  tibble(
    barcodes = colnames(count_mtr),
    nCount_Spatial = as.numeric(n_count),
    nFeature_Spatial = as.numeric(n_feature),
    pct_mt = as.numeric(pct_mt)
  )
}

compute_embedding <- function(expr_mtr, n_neighbors = 15L, seed = SPATA_SEED) {
  set.seed(seed)
  gene_var <- apply(as.matrix(expr_mtr), 1, stats::var)
  keep_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(1500, length(gene_var)))]
  emb_input <- t(as.matrix(expr_mtr[keep_genes, , drop = FALSE]))
  if (nrow(emb_input) < 3) {
    return(tibble(UMAP1 = rep(0, nrow(emb_input)), UMAP2 = rep(0, nrow(emb_input))))
  }
  pca <- stats::prcomp(emb_input, center = TRUE, scale. = TRUE)
  pcs <- pca$x[, seq_len(min(20, ncol(pca$x))), drop = FALSE]
  umap <- uwot::umap(pcs, n_neighbors = min(n_neighbors, max(2, nrow(pcs) - 1)), min_dist = 0.25, metric = "cosine", ret_model = FALSE, verbose = FALSE)
  tibble(UMAP1 = umap[, 1], UMAP2 = umap[, 2])
}

build_spata_object <- function(sample_id, count_mtr, coords_df, metadata_df) {
  spatial_method <- if (ncol(count_mtr) <= 7000) "VisiumSmall" else "VisiumLarge"
  sp_obj <- SPATA2::initiateSpataObject(
    sample_name = sample_id,
    count_mtr = count_mtr,
    modality = "gene",
    coords_df = coords_df,
    spatial_method = spatial_method,
    verbose = FALSE
  )
  sp_obj <- SPATA2::normalizeCounts(sp_obj, method = "LogNormalize")
  sp_obj
}

load_spatial_sample <- function(path, sample_id, dataset, group_label) {
  sce <- readH5AD(path)
  assay_name <- resolve_counts_assay(sce)
  count_mtr <- assay(sce, assay_name)
  if (!inherits(count_mtr, "dgCMatrix")) count_mtr <- as(count_mtr, "dgCMatrix")
  gene_names <- rownames(count_mtr)
  if (is.null(gene_names) || all(gene_names == "")) gene_names <- rownames(sce)
  gene_names <- make.unique(toupper(gene_names))
  rownames(count_mtr) <- gene_names

  coords_df <- extract_coords_df(sce)
  common <- intersect(colnames(count_mtr), coords_df$barcodes)
  count_mtr <- count_mtr[, common, drop = FALSE]
  coords_df <- coords_df[common, , drop = FALSE]
  meta <- as.data.frame(colData(sce))
  meta$barcodes <- colnames(sce)
  meta <- meta[match(common, meta$barcodes), , drop = FALSE]
  meta$sample_id <- sample_id
  meta$dataset <- dataset
  meta$group <- group_label
  meta$barcodes <- common

  qc_df <- compute_qc_metrics(count_mtr)
  norm_mtr <- log_normalize_counts(count_mtr)
  sig_df <- module_score_matrix(norm_mtr, IO_SIGNATURES)
  cell_score_df <- module_score_matrix(norm_mtr, CELL_TYPE_MARKERS)
  annot_df <- annotate_celltype_confidence(cell_score_df[common, , drop = FALSE])
  cell_df <- cell_score_df
  colnames(cell_df) <- sub("^", "ct_", colnames(cell_df))
  embedding_df <- compute_embedding(norm_mtr)

  spot_df <- meta %>%
    dplyr::select(barcodes, sample_id, dataset, group) %>%
    bind_cols(coords_df[, c("x", "y")]) %>%
    left_join(qc_df, by = "barcodes") %>%
    bind_cols(sig_df[common, , drop = FALSE]) %>%
    bind_cols(cell_df[common, , drop = FALSE]) %>%
    bind_cols(annot_df) %>%
    bind_cols(embedding_df) %>%
    mutate(
      spatial_immunotype = classify_immunotype(sig_df[common, , drop = FALSE]),
      immune_gradient = T_cell_activation - TGF_beta
    )

  spata_obj <- build_spata_object(sample_id, count_mtr, coords_df, spot_df)

  list(
    sample_id = sample_id,
    dataset = dataset,
    group = group_label,
    counts = count_mtr,
    norm = norm_mtr,
    coords = coords_df,
    spot_df = spot_df,
    spata = spata_obj
  )
}

plot_spatial_numeric <- function(df, feature, title = feature, subtitle = NULL, palette = "C") {
  plot_df <- df %>% mutate(plot_value = clip_numeric(.data[[feature]]))
  point_size <- adaptive_point_size(nrow(plot_df), min_size = 0.45, max_size = 1.5)
  ggplot(plot_df, aes(x = x, y = y, colour = plot_value)) +
    geom_point(size = point_size, alpha = 0.9, stroke = 0, shape = 16) +
    coord_fixed() +
    scale_y_reverse() +
    scale_colour_viridis_c(option = palette, na.value = "grey80") +
    labs(title = title, subtitle = subtitle, colour = NULL) +
    theme_spatial()
}

plot_spatial_discrete <- function(df, feature, palette, title = feature, subtitle = NULL) {
  point_size <- adaptive_point_size(nrow(df), min_size = 0.45, max_size = 1.5)
  ggplot(df, aes(x = x, y = y, colour = .data[[feature]])) +
    geom_point(size = point_size, alpha = 0.95, stroke = 0, shape = 16) +
    coord_fixed() +
    scale_y_reverse() +
    scale_colour_manual(values = palette, na.value = "grey80", drop = FALSE) +
    labs(title = title, subtitle = subtitle, colour = NULL) +
    theme_spatial()
}

plot_umap_by_feature <- function(df, feature, discrete = FALSE, palette = NULL, title = feature) {
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2))
  point_size <- adaptive_point_size(nrow(df), min_size = 0.25, max_size = 0.8)
  if (discrete) {
    p <- p + geom_point(aes(colour = .data[[feature]]), size = point_size, alpha = 0.85) +
      scale_colour_manual(values = palette, drop = FALSE, na.value = "grey80")
  } else {
    p <- p + geom_point(aes(colour = clip_numeric(.data[[feature]])), size = point_size, alpha = 0.85) +
      scale_colour_viridis_c(option = "magma", na.value = "grey80")
  }
  p + labs(title = title, colour = NULL) + theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
}

compute_knn_graph <- function(coords_df, k = 6L) {
  knn <- RANN::nn2(data = coords_df[, c("x", "y")], query = coords_df[, c("x", "y")], k = k + 1)
  knn$nn.idx[, -1, drop = FALSE]
}

moran_stat <- function(x, nn_idx) {
  ok <- is.finite(x)
  if (sum(ok) < 3) return(c(moran_I = NA_real_, p_perm = NA_real_))
  x <- x[ok]
  nn_idx <- nn_idx[ok, , drop = FALSE]
  remap <- seq_len(length(ok))[ok]
  lookup <- setNames(seq_along(remap), remap)
  valid <- matrix(NA_integer_, nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  for (i in seq_len(nrow(nn_idx))) {
    valid[i, ] <- suppressWarnings(as.integer(lookup[as.character(nn_idx[i, ])]))
  }
  centered <- x - mean(x)
  w_sum <- sum(!is.na(valid))
  if (w_sum == 0) return(c(moran_I = NA_real_, p_perm = NA_real_))
  numerator <- 0
  for (i in seq_len(nrow(valid))) {
    nbrs <- valid[i, !is.na(valid[i, ])]
    numerator <- numerator + sum(centered[i] * centered[nbrs])
  }
  denom <- sum(centered^2)
  moran_i <- length(x) / w_sum * numerator / denom
  perm_vals <- replicate(100, {
    xp <- sample(centered)
    num <- 0
    for (i in seq_len(nrow(valid))) {
      nbrs <- valid[i, !is.na(valid[i, ])]
      num <- num + sum(xp[i] * xp[nbrs])
    }
    length(xp) / w_sum * num / sum(xp^2)
  })
  p_perm <- mean(abs(perm_vals) >= abs(moran_i))
  c(moran_I = moran_i, p_perm = p_perm)
}

cooccurrence_table <- function(labels, nn_idx) {
  pair_list <- vector("list", nrow(nn_idx))
  for (i in seq_len(nrow(nn_idx))) {
    nbrs <- labels[nn_idx[i, ]]
    pair_list[[i]] <- tibble(source = labels[[i]], target = nbrs)
  }
  bind_rows(pair_list) %>% count(source, target, name = "n")
}

summarise_sample <- function(sample_obj) {
  sample_df <- sample_obj$spot_df
  row <- tibble(
    sample = sample_obj$sample_id,
    dataset = sample_obj$dataset,
    group = sample_obj$group,
    n_spots = nrow(sample_df),
    median_nCount = median(sample_df$nCount_Spatial, na.rm = TRUE),
    median_nFeature = median(sample_df$nFeature_Spatial, na.rm = TRUE),
    median_pct_mt = median(sample_df$pct_mt, na.rm = TRUE)
  )
  for (sig in names(IO_SIGNATURES)) {
    if (sig %in% colnames(sample_df)) {
      row[[paste0("mean_", sig)]] <- mean(sample_df[[sig]], na.rm = TRUE)
      row[[paste0("sd_", sig)]] <- stats::sd(sample_df[[sig]], na.rm = TRUE)
    }
  }
  for (label in IMMUNOTYPE_LEVELS) {
    row[[paste0("frac_", label)]] <- mean(sample_df$spatial_immunotype == label, na.rm = TRUE)
  }
  row
}

pairwise_group_stats <- function(df, group_col, groups, value_cols, id_cols = character()) {
  bind_rows(lapply(value_cols, function(col_name) {
    sub_df <- df %>% filter(.data[[group_col]] %in% groups)
    g1 <- sub_df %>% filter(.data[[group_col]] == groups[[1]]) %>% pull(.data[[col_name]])
    g2 <- sub_df %>% filter(.data[[group_col]] == groups[[2]]) %>% pull(.data[[col_name]])
    if (sum(is.finite(g1)) < 3 || sum(is.finite(g2)) < 3) {
      return(tibble(feature = col_name, group_1 = groups[[1]], group_2 = groups[[2]], p_value = NA_real_, effect = NA_real_))
    }
    test <- suppressWarnings(stats::wilcox.test(g1, g2, exact = FALSE))
    tibble(
      feature = col_name,
      group_1 = groups[[1]],
      group_2 = groups[[2]],
      p_value = test$p.value,
      effect = median(g1, na.rm = TRUE) - median(g2, na.rm = TRUE)
    )
  }))
}

sample_level_permutation <- function(sample_stats, group_col, value_cols, n_perm = 1000L) {
  set.seed(SPATA_SEED)
  bind_rows(lapply(value_cols, function(col_name) {
    sub_df <- sample_stats %>% filter(is.finite(.data[[col_name]]))
    if (n_distinct(sub_df[[group_col]]) != 2) return(NULL)
    groups <- unique(sub_df[[group_col]])
    values <- sub_df[[col_name]]
    group_values <- sub_df[[group_col]]
    obs_diff <- mean(values[group_values == groups[[1]]]) - mean(values[group_values == groups[[2]]])
    perm_diffs <- replicate(n_perm, {
      perm_group <- sample(group_values)
      mean(values[perm_group == groups[[1]]]) - mean(values[perm_group == groups[[2]]])
    })
    tibble(
      feature = col_name,
      group_1 = groups[[1]],
      group_2 = groups[[2]],
      observed_diff = obs_diff,
      permutation_p = mean(abs(perm_diffs) >= abs(obs_diff))
    )
  }))
}

leave_one_out_stats <- function(sample_stats, group_col, value_cols) {
  bind_rows(lapply(seq_len(nrow(sample_stats)), function(i) {
    loo_df <- sample_stats[-i, , drop = FALSE]
    perms <- sample_level_permutation(loo_df, group_col = group_col, value_cols = value_cols, n_perm = 500L)
    if (nrow(perms) == 0) return(NULL)
    perms$left_out_sample <- sample_stats$sample[[i]]
    perms
  }))
}

gene_order_from_gtf <- function(gtf_path, genes) {
  con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, open = "rt") else file(gtf_path, open = "rt")
  on.exit(close(con), add = TRUE)
  raw_lines <- readLines(con, warn = FALSE)
  raw_lines <- raw_lines[!startsWith(raw_lines, "#")]
  split_lines <- strsplit(raw_lines, "\t")
  keep <- vapply(split_lines, function(x) length(x) >= 9 && x[[3]] == "gene", logical(1))
  split_lines <- split_lines[keep]
  gene_df <- tibble(
    seqname = vapply(split_lines, `[[`, character(1), 1),
    start = as.integer(vapply(split_lines, `[[`, character(1), 4)),
    end = as.integer(vapply(split_lines, `[[`, character(1), 5)),
    attribute = vapply(split_lines, `[[`, character(1), 9)
  ) %>%
    mutate(
      gene_name = str_match(attribute, 'gene_name "([^"]+)"')[, 2],
      gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
      gene_symbol = toupper(coalesce(gene_name, gene_id))
    ) %>%
    filter(seqname %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  genes_upper <- toupper(genes)
  gene_df <- gene_df %>% filter(gene_symbol %in% genes_upper)
  gene_df[match(genes_upper[genes_upper %in% gene_df$gene_symbol], gene_df$gene_symbol), c("gene_symbol", "seqname", "start", "end")]
}

infer_reference_group <- function(spot_df) {
  ref_group <- ifelse(as.character(spot_df$dominant_cell_type) %in% REFERENCE_CELL_TYPES, "reference", "observation")
  factor(ref_group, levels = c("reference", "observation"))
}

run_spatial_infercnv <- function(sample_obj, out_prefix, gene_order_df, cutoff = 0.1) {
  spot_df <- sample_obj$spot_df
  ref_group <- infer_reference_group(spot_df)
  if (sum(ref_group == "reference") < 20 || sum(ref_group == "observation") < 20) {
    stop(sprintf("Sample %s does not have enough reference/observation spots for infercnv.", sample_obj$sample_id), call. = FALSE)
  }

  annotations_df <- data.frame(barcode = colnames(sample_obj$counts), group = ref_group, stringsAsFactors = FALSE)
  annotation_path <- paste0(out_prefix, "_annotations.tsv")
  gene_order_path <- paste0(out_prefix, "_gene_order.tsv")
  write.table(annotations_df, annotation_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(gene_order_df, gene_order_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infer_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = sample_obj$counts,
    gene_order_file = gene_order_path,
    annotations_file = annotation_path,
    ref_group_names = "reference",
    min_max_counts_per_cell = c(100, Inf)
  )

  infercnv::run(
    infercnv_obj = infer_obj,
    cutoff = cutoff,
    out_dir = dirname(out_prefix),
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    analysis_mode = "cells",
    no_plot = FALSE,
    plot_steps = FALSE,
    resume_mode = FALSE,
    png_res = 600,
    output_format = "png",
    num_threads = 4,
    save_final_rds = TRUE
  )
}

infercnv_expr_to_arm_table <- function(expr_mtr, gene_order_df, sample_id) {
  gene_tbl <- tibble(
    gene_symbol = toupper(rownames(expr_mtr)),
    row_id = rownames(expr_mtr)
  ) %>% inner_join(gene_order_df, by = "gene_symbol")
  if (nrow(gene_tbl) == 0) return(NULL)

  arm_tbl <- gene_tbl %>%
    group_by(seqname) %>%
    mutate(chr_mid = max(end, na.rm = TRUE) / 2) %>%
    ungroup() %>%
    mutate(
      chr = sub("^chr", "", seqname),
      arm = paste0(chr, ifelse(start <= chr_mid, "p", "q"))
    )

  arm_scores <- lapply(unique(arm_tbl$arm), function(arm_name) {
    arm_genes <- arm_tbl$row_id[arm_tbl$arm == arm_name]
    tibble(
      barcode = colnames(expr_mtr),
      arm = arm_name,
      cnv = Matrix::colMeans(expr_mtr[arm_genes, , drop = FALSE] - 1)
    )
  })
  bind_rows(arm_scores) %>% mutate(sample_id = sample_id)
}

plot_arm_heatmap <- function(arm_summary, file_stub) {
  if (nrow(arm_summary) == 0) return(invisible(NULL))

  # Order arms by chromosome then p/q arm
  arm_levels <- unique(arm_summary$arm)
  chr_num <- suppressWarnings(as.integer(sub("[pq]$", "", arm_levels)))
  arm_type <- ifelse(grepl("p$", arm_levels), 0L, 1L)
  arm_levels <- arm_levels[order(chr_num, arm_type)]

  sample_levels <- sort(unique(arm_summary$sample_id))
  heatmap_df <- arm_summary %>%
    mutate(arm = factor(arm, levels = arm_levels),
           sample_id = factor(sample_id, levels = sample_levels)) %>%
    tidyr::pivot_wider(names_from = arm, values_from = mean_cnv, values_fill = 0) %>%
    as.data.frame()
  rownames(heatmap_df) <- as.character(heatmap_df$sample_id)
  heatmap_mat <- as.matrix(heatmap_df[, arm_levels[arm_levels %in% colnames(heatmap_df)], drop = FALSE])

  # Chromosome groupings for column splits
  arm_chr <- sub("[pq]$", "", colnames(heatmap_mat))
  chr_factor <- factor(arm_chr, levels = unique(arm_chr))
  n_chrs <- length(unique(arm_chr))

  col_fun <- circlize::colorRamp2(
    c(-0.15, 0, 0.15),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  # Row annotation: response group if column exists
  row_annot <- NULL
  if ("response" %in% colnames(arm_summary)) {
    resp_map <- arm_summary %>%
      select(sample_id, response) %>%
      distinct() %>%
      tibble::deframe()
    resp_vec <- resp_map[rownames(heatmap_mat)]
    row_annot <- rowAnnotation(
      Response = resp_vec,
      col = list(Response = c(Responder = "#0077BB", Non_Responder = "#CC3311",
                               Tumor = "#E41A1C", Normal = "#4DAF4A",
                               Primary = "#FF7F00", Recurrence = "#984EA3")),
      annotation_name_gp = gpar(fontsize = 9),
      simple_anno_size = unit(4, "mm"),
      show_legend = TRUE
    )
  }

  top_annot <- HeatmapAnnotation(
    Chromosome = anno_block(
      gp = gpar(fill = rep(c("#E8E8E8", "#C8C8C8"), length.out = n_chrs)),
      labels = levels(chr_factor),
      labels_gp = gpar(cex = 0.65, fontface = "bold"),
      height = unit(4, "mm")
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 0)
  )

  ht <- Heatmap(
    heatmap_mat,
    name = "CNV\ndeviation",
    col = col_fun,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = FALSE,
    cluster_rows = nrow(heatmap_mat) > 1,
    cluster_columns = FALSE,
    column_split = chr_factor,
    column_gap = unit(0.5, "mm"),
    top_annotation = top_annot,
    right_annotation = row_annot,
    use_raster = TRUE,
    raster_device = "png",
    raster_quality = 3,
    border = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    column_title = "Chromosomal arm",
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )

  out_path_png <- file.path(CNV_FIG_DIR, paste0(file_stub, ".png"))
  out_path_pdf <- file.path(CNV_FIG_DIR, paste0(file_stub, ".pdf"))
  # Wide format: ~200px per arm, minimum 2000px wide
  n_arms <- ncol(heatmap_mat)
  pw <- max(2000L, n_arms * 45L)
  ph <- max(800L, nrow(heatmap_mat) * 60L + 400L)
  png(out_path_png, width = pw, height = ph, res = 220)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  pdf(out_path_pdf, width = pw / 220, height = ph / 220)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat(sprintf("  Arm heatmap: %s\n", out_path_png))
  invisible(ht)
}

# High-resolution full-genome CNV heatmap (grouped by dominant cell type)
# Maps each infercnv expr.data matrix across chromosome positions using
# ComplexHeatmap with chromosome band annotations — mirrors reference pipeline
# 05_cnv_genome_heatmap.R
plot_genome_cnv_heatmap <- function(infercnv_list, spot_df_list, gene_order_df, file_stub,
                                     group_col = "dominant_cell_type",
                                     max_cells_per_group = 500L,
                                     plot_width_px = 7000L,
                                     plot_height_px = 2000L,
                                     plot_dpi = 300L) {
  # infercnv_list: named list of infercnv objects (per sample)
  # spot_df_list: named list of spot data frames with barcodes + group_col
  # gene_order_df: from gene_order_from_gtf()

  # 1. Determine common genes across all samples using gene_order_df
  common_genes <- gene_order_df$gene_symbol
  for (nm in names(infercnv_list)) {
    common_genes <- intersect(common_genes, rownames(infercnv_list[[nm]]@expr.data))
  }
  if (length(common_genes) == 0) {
    warning("plot_genome_cnv_heatmap: no common genes across infercnv objects.")
    return(invisible(NULL))
  }

  # Use gene_order_df for positional ordering
  gene_df <- gene_order_df %>%
    filter(gene_symbol %in% common_genes) %>%
    mutate(
      chr_num = suppressWarnings(as.integer(sub("^chr", "", seqname))),
      chr_label = sub("^chr", "", seqname)
    ) %>%
    filter(!is.na(chr_num), chr_num >= 1L, chr_num <= 22L) %>%
    arrange(chr_num, start)
  ordered_genes <- gene_df$gene_symbol

  if (length(ordered_genes) < 100) {
    warning("plot_genome_cnv_heatmap: too few ordered genes after GTF filtering.")
    return(invisible(NULL))
  }

  # 2. Aggregate expression by group_col across all samples
  all_groups <- sort(unique(unlist(lapply(spot_df_list, function(df) {
    if (group_col %in% colnames(df)) as.character(df[[group_col]]) else character(0)
  }))))

  group_agg <- list()
  for (nm in names(infercnv_list)) {
    infer_obj <- infercnv_list[[nm]]
    spot_df <- spot_df_list[[nm]]
    if (is.null(spot_df) || !(group_col %in% colnames(spot_df))) next

    expr <- infer_obj@expr.data
    expr <- expr[ordered_genes[ordered_genes %in% rownames(expr)], , drop = FALSE]
    if (nrow(expr) == 0) next

    for (grp in all_groups) {
      bc <- as.character(spot_df$barcodes[as.character(spot_df[[group_col]]) == grp])
      bc <- intersect(bc, colnames(expr))
      if (length(bc) == 0) next

      # Downsample if needed
      if (length(bc) > max_cells_per_group) {
        set.seed(SPATA_SEED)
        bc <- sample(bc, max_cells_per_group)
      }

      sub_expr <- expr[, bc, drop = FALSE]
      mean_vec <- Matrix::rowMeans(sub_expr, na.rm = TRUE)

      if (is.null(group_agg[[grp]])) {
        group_agg[[grp]] <- list(sum = mean_vec, n = 1L)
      } else {
        # Average across samples
        group_agg[[grp]]$sum <- group_agg[[grp]]$sum + mean_vec
        group_agg[[grp]]$n <- group_agg[[grp]]$n + 1L
      }
    }
  }

  if (length(group_agg) == 0) {
    warning("plot_genome_cnv_heatmap: no groups could be aggregated.")
    return(invisible(NULL))
  }

  # Build matrix: rows = groups, cols = genes ordered by genome position
  present_genes <- rownames(group_agg[[1]]$sum)
  gene_order_final <- ordered_genes[ordered_genes %in% present_genes]
  gene_df_final <- gene_df %>% filter(gene_symbol %in% gene_order_final)
  # Re-order to match gene_order_final
  gene_df_final <- gene_df_final[match(gene_order_final, gene_df_final$gene_symbol), ]

  mat <- do.call(rbind, lapply(names(group_agg), function(g) {
    v <- group_agg[[g]]$sum / group_agg[[g]]$n
    v[gene_order_final]
  }))
  rownames(mat) <- names(group_agg)
  colnames(mat) <- gene_order_final

  # Center at 1 (inferCNV expresses CNV as values around 1)
  mat_centered <- mat - 1

  # 3. Chromosome annotation
  chr_labels <- as.character(1:22)
  chr_factor <- factor(gene_df_final$chr_num, levels = as.integer(chr_labels))
  chr_labels_present <- chr_labels[chr_labels %in% as.character(levels(chr_factor)[!is.na(levels(chr_factor))])]

  col_fun <- circlize::colorRamp2(
    c(-0.15, 0, 0.15),
    c("#2166AC", "#F7F7F7", "#B2182B")
  )

  top_annot <- HeatmapAnnotation(
    Chromosome = anno_block(
      gp = gpar(fill = rep(c("#f1a340", "#998ec3"), length.out = 22)),
      labels = as.character(1:22),
      labels_gp = gpar(cex = 0.65, fontface = "bold"),
      height = unit(5, "mm")
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 0)
  )

  ht <- Heatmap(
    mat_centered,
    name = "CNV",
    col = col_fun,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    row_names_side = "left",
    show_column_names = FALSE,
    cluster_rows = nrow(mat_centered) > 1,
    cluster_columns = FALSE,
    column_split = chr_factor,
    column_gap = unit(0.8, "mm"),
    top_annotation = top_annot,
    use_raster = TRUE,
    raster_device = "png",
    raster_quality = 3,
    border = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    column_title = NULL
  )

  out_path_png <- file.path(CNV_FIG_DIR, paste0(file_stub, ".png"))
  out_path_pdf <- file.path(CNV_FIG_DIR, paste0(file_stub, ".pdf"))
  png(out_path_png, width = plot_width_px, height = plot_height_px, res = plot_dpi)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  pdf(out_path_pdf, width = plot_width_px / plot_dpi, height = plot_height_px / plot_dpi)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat(sprintf("  Genome CNV heatmap: %s\n", out_path_png))
  invisible(ht)
}

plot_cnv_spatial_map <- function(df, feature, file_stub, title) {
  p <- plot_spatial_numeric(df, feature = feature, title = title, palette = "magma")
  safe_save(file.path(CNV_FIG_DIR, paste0(file_stub, ".png")), p, width = 7, height = 6, dpi = 500)
}

make_signature_coverage_table <- function(sample_obj, gene_sets, family) {
  gene_names <- toupper(rownames(sample_obj$counts))
  bind_rows(lapply(names(gene_sets), function(set_name) {
    genes <- unique(toupper(gene_sets[[set_name]]))
    detected <- intersect(genes, gene_names)
    tibble(
      sample_id = sample_obj$sample_id,
      dataset = sample_obj$dataset,
      group = sample_obj$group,
      signature_family = family,
      signature_name = set_name,
      n_genes_defined = length(genes),
      n_genes_detected = length(detected),
      frac_genes_detected = length(detected) / max(length(genes), 1),
      detected_genes = paste(detected, collapse = ";")
    )
  }))
}

make_spot_qc_catalog <- function(sample_objects) {
  bind_rows(lapply(sample_objects, function(sample_obj) {
    sample_obj$spot_df %>%
      transmute(
        sample_id,
        dataset,
        group,
        barcodes,
        x,
        y,
        nCount_Spatial,
        nFeature_Spatial,
        pct_mt,
        dominant_cell_type,
        dominant_score,
        margin_score,
        annotation_confidence,
        spatial_immunotype,
        immune_gradient
      )
  })) %>%
    mutate(umi_per_feature = nCount_Spatial / pmax(nFeature_Spatial, 1))
}

summarise_spot_qc_catalog <- function(spot_catalog) {
  spot_catalog %>%
    group_by(sample_id, dataset, group) %>%
    summarise(
      n_spots = n(),
      median_nCount = median(nCount_Spatial, na.rm = TRUE),
      median_nFeature = median(nFeature_Spatial, na.rm = TRUE),
      median_pct_mt = median(pct_mt, na.rm = TRUE),
      median_umi_per_feature = median(umi_per_feature, na.rm = TRUE),
      frac_high_confidence = mean(annotation_confidence == "High", na.rm = TRUE),
      frac_low_confidence = mean(annotation_confidence == "Low", na.rm = TRUE),
      n_cell_types = n_distinct(dominant_cell_type),
      n_immunotypes = n_distinct(spatial_immunotype),
      .groups = "drop"
    )
}

neighbor_enrichment_table <- function(labels, nn_idx) {
  labels <- as.character(labels)
  pair_tbl <- cooccurrence_table(labels, nn_idx)
  if (nrow(pair_tbl) == 0) return(pair_tbl)
  label_freq <- prop.table(table(labels))
  total_edges <- sum(pair_tbl$n)
  pair_tbl %>%
    mutate(
      expected = as.numeric(label_freq[source]) * as.numeric(label_freq[target]) * total_edges,
      observed_fraction = n / total_edges,
      log2_enrichment = log2((n + 1) / (expected + 1))
    )
}

adjust_stat_table <- function(stat_tbl, p_col = "p_value") {
  if (nrow(stat_tbl) == 0 || !(p_col %in% colnames(stat_tbl))) return(stat_tbl)
  stat_tbl %>% mutate(p_adj = stats::p.adjust(.data[[p_col]], method = "BH"))
}

write_note_lines <- function(path, lines) {
  ensure_dirs(dirname(path))
  writeLines(lines, con = path)
}

make_genome_bin_table <- function(gene_order_df, bin_size = 75L) {
  gene_order_df %>%
    mutate(
      chr_num = suppressWarnings(as.integer(sub("^chr", "", seqname))),
      chr_label = sub("^chr", "", seqname)
    ) %>%
    filter(!is.na(chr_num), chr_num >= 1L, chr_num <= 22L) %>%
    arrange(chr_num, start, end) %>%
    group_by(chr_label) %>%
    mutate(
      chr_bin_id = ((row_number() - 1L) %/% bin_size) + 1L,
      genome_bin = paste0(chr_label, "_", sprintf("%03d", chr_bin_id))
    ) %>%
    ungroup() %>%
    select(gene_symbol, seqname, start, end, chr_num, chr_label, genome_bin)
}

bin_infercnv_matrix <- function(expr_mtr, bin_tbl) {
  gene_names <- toupper(rownames(expr_mtr))
  rownames(expr_mtr) <- gene_names
  map_tbl <- bin_tbl %>% filter(gene_symbol %in% gene_names)
  if (nrow(map_tbl) == 0) return(NULL)
  map_tbl <- map_tbl[match(intersect(map_tbl$gene_symbol, rownames(expr_mtr)), map_tbl$gene_symbol), , drop = FALSE]
  expr_sub <- expr_mtr[map_tbl$gene_symbol, , drop = FALSE]
  bin_levels <- unique(map_tbl$genome_bin)
  split_idx <- split(seq_len(nrow(map_tbl)), map_tbl$genome_bin)
  bin_mat <- do.call(cbind, lapply(bin_levels, function(bin_id) {
    idx <- split_idx[[bin_id]]
    Matrix::colMeans(expr_sub[idx, , drop = FALSE] - 1, na.rm = TRUE)
  }))
  if (is.null(dim(bin_mat))) {
    bin_mat <- matrix(bin_mat, nrow = ncol(expr_sub), dimnames = list(colnames(expr_sub), bin_levels))
  } else {
    rownames(bin_mat) <- colnames(expr_sub)
    colnames(bin_mat) <- bin_levels
  }
  list(mat = bin_mat, bin_meta = distinct(map_tbl, genome_bin, chr_label, chr_num))
}

plot_binned_cnv_heatmap <- function(infercnv_list,
                                    spot_df_list,
                                    gene_order_df,
                                    file_stub,
                                    max_spots_per_sample = 250L,
                                    bin_size = 75L,
                                    row_annotation_cols = c("group", "dominant_cell_type", "spatial_immunotype"),
                                    plot_dpi = 300L) {
  bin_tbl <- make_genome_bin_table(gene_order_df, bin_size = bin_size)
  if (nrow(bin_tbl) == 0) return(invisible(NULL))

  heatmap_parts <- list()
  row_meta_parts <- list()
  reference_bin_meta <- NULL

  for (sample_id in names(infercnv_list)) {
    infer_obj <- infercnv_list[[sample_id]]
    spot_df <- spot_df_list[[sample_id]]
    if (is.null(infer_obj) || is.null(spot_df)) next
    expr_mtr <- infer_obj@expr.data
    bc <- intersect(colnames(expr_mtr), as.character(spot_df$barcodes))
    if (length(bc) < 10) next
    set.seed(SPATA_SEED)
    if (length(bc) > max_spots_per_sample) bc <- sample(bc, max_spots_per_sample)
    expr_mtr <- expr_mtr[, bc, drop = FALSE]
    binned <- bin_infercnv_matrix(expr_mtr, bin_tbl)
    if (is.null(binned)) next
    if (is.null(reference_bin_meta)) reference_bin_meta <- binned$bin_meta

    meta_df <- spot_df %>%
      filter(barcodes %in% bc) %>%
      mutate(barcodes = as.character(barcodes)) %>%
      distinct(barcodes, .keep_all = TRUE) %>%
      select(any_of(c("barcodes", "sample_id", row_annotation_cols)))
    meta_df <- meta_df[match(rownames(binned$mat), meta_df$barcodes), , drop = FALSE]
    meta_df$sample_id <- sample_id
    meta_df$cnv_burden_binned <- rowMeans(abs(binned$mat), na.rm = TRUE)

    ord_cols <- row_annotation_cols[row_annotation_cols %in% colnames(meta_df)]
    ord <- if (length(ord_cols) > 0) {
      do.call(order, c(meta_df[ord_cols], list(-meta_df$cnv_burden_binned)))
    } else {
      order(-meta_df$cnv_burden_binned)
    }

    heatmap_parts[[sample_id]] <- binned$mat[ord, , drop = FALSE]
    row_meta_parts[[sample_id]] <- meta_df[ord, , drop = FALSE]
  }

  if (length(heatmap_parts) == 0) return(invisible(NULL))

  heatmap_mat <- do.call(rbind, heatmap_parts)
  row_meta <- bind_rows(row_meta_parts)
  stopifnot(identical(rownames(heatmap_mat), row_meta$barcodes))

  sample_levels <- unique(row_meta$sample_id)
  sample_palette <- setNames(
    colorRampPalette(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"))(length(sample_levels)),
    sample_levels
  )

  annot_fields <- list(Sample = factor(row_meta$sample_id, levels = sample_levels))
  annot_cols <- list(Sample = sample_palette)
  if ("group" %in% colnames(row_meta)) {
    group_levels <- unique(row_meta$group)
    annot_fields$Group <- factor(row_meta$group, levels = group_levels)
    annot_cols$Group <- c(RESPONSE_PALETTE, TISSUE_PALETTE)[names(c(RESPONSE_PALETTE, TISSUE_PALETTE)) %in% group_levels]
  }
  if ("dominant_cell_type" %in% colnames(row_meta)) {
    annot_fields$Celltype <- factor(row_meta$dominant_cell_type, levels = names(CELLTYPE_PALETTE))
    annot_cols$Celltype <- CELLTYPE_PALETTE
  }
  if ("spatial_immunotype" %in% colnames(row_meta)) {
    annot_fields$Immunotype <- factor(row_meta$spatial_immunotype, levels = IMMUNOTYPE_LEVELS)
    annot_cols$Immunotype <- IMMUNOTYPE_PALETTE
  }
  row_ha <- do.call(ComplexHeatmap::rowAnnotation, c(annot_fields, list(col = annot_cols, show_annotation_name = TRUE)))

  chr_factor <- factor(
    reference_bin_meta$chr_label[match(colnames(heatmap_mat), reference_bin_meta$genome_bin)],
    levels = unique(reference_bin_meta$chr_label[order(reference_bin_meta$chr_num)])
  )
  top_annot <- HeatmapAnnotation(
    Chromosome = anno_block(
      gp = gpar(fill = rep(c("#f1a340", "#998ec3"), length.out = nlevels(chr_factor))),
      labels = levels(chr_factor),
      labels_gp = gpar(cex = 0.65, fontface = "bold"),
      height = unit(4, "mm")
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 0)
  )
  col_fun <- circlize::colorRamp2(c(-0.2, 0, 0.2), c("#2166AC", "#F7F7F7", "#B2182B"))

  ht <- Heatmap(
    heatmap_mat,
    name = "CNV",
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(row_meta$sample_id, levels = sample_levels),
    column_split = chr_factor,
    left_annotation = row_ha,
    top_annotation = top_annot,
    use_raster = TRUE,
    raster_device = "png",
    raster_quality = 3,
    border = FALSE,
    column_gap = unit(0.7, "mm")
  )

  plot_width_px <- max(8000L, ncol(heatmap_mat) * 26L)
  plot_height_px <- max(2400L, nrow(heatmap_mat) * 3L + 800L)
  png_path <- file.path(CNV_FIG_DIR, paste0(file_stub, ".png"))
  pdf_path <- file.path(CNV_FIG_DIR, paste0(file_stub, ".pdf"))
  png(png_path, width = plot_width_px, height = plot_height_px, res = plot_dpi)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  pdf(pdf_path, width = plot_width_px / plot_dpi, height = plot_height_px / plot_dpi)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  write_csv(row_meta, file.path(CNV_TAB_DIR, paste0(file_stub, "_rows.csv")))
  write_csv(reference_bin_meta, file.path(CNV_TAB_DIR, paste0(file_stub, "_bins.csv")))
  invisible(ht)
}