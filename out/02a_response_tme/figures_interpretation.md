# Figures interpretation — `out/02a_response_tme`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Tumour-microenvironment composition and signature contrasts between Responders (R, n=7) and Non-Responders (NR, n=17) in the GSE206325 anti-PD-1 discovery cohort. Patient-level Mann–Whitney U with BH-FDR; effect sizes reported as Cliff's δ (positive = higher in R).

## Figure-by-figure interpretation

### `figures/`
#### `figures/02_abundance_barplot.png`

Legacy duplicate of the per-cell-type abundance bars (kept for traceability).

![02_abundance_barplot](figures/02_abundance_barplot.png)
#### `figures/02_activation_exhaustion_scatter.png`

Per-patient scatter of activation vs exhaustion signature means — Responders cluster in a 'high-activation-with-productive-dysfunction' region.

![02_activation_exhaustion_scatter](figures/02_activation_exhaustion_scatter.png)
#### `figures/02_celltype_response_umap.png`

Legacy duplicate of the response-coloured cell-type UMAP.

![02_celltype_response_umap](figures/02_celltype_response_umap.png)
#### `figures/02_composition_stacked_bar.png`

Legacy duplicate of the stacked composition bar plot.

![02_composition_stacked_bar](figures/02_composition_stacked_bar.png)
#### `figures/02_signature_radar.png`

Legacy duplicate of the IO signature radar plot.

![02_signature_radar](figures/02_signature_radar.png)
#### `figures/02_signature_violins.png`

Legacy duplicate of the per-patient signature violins.

![02_signature_violins](figures/02_signature_violins.png)
#### `figures/02_tissue_signature_heatmap.png`

Legacy duplicate of the tissue-level signature heatmap.

![02_tissue_signature_heatmap](figures/02_tissue_signature_heatmap.png)
#### `figures/02_volcano_deg.png`

Legacy duplicate of the pseudo-bulk DE volcano.

![02_volcano_deg](figures/02_volcano_deg.png)
#### `figures/02a_abundance_barplot.png`

Per-cell-type abundance bars for R vs NR with BH-adjusted Mann–Whitney p-values.

![02a_abundance_barplot](figures/02a_abundance_barplot.png)
#### `figures/02a_activation_exhaustion_scatter.png`

Auxiliary diagnostic figure produced by the corresponding analysis script. See `scripts/` and `results_index.yaml` for provenance.

![02a_activation_exhaustion_scatter](figures/02a_activation_exhaustion_scatter.png)
#### `figures/02a_celltype_annotation_enhanced.png`

Enhanced annotation panel combining marker dot plot, UMAP and per-patient prevalence.

![02a_celltype_annotation_enhanced](figures/02a_celltype_annotation_enhanced.png)
#### `figures/02a_celltype_response_umap.png`

Cell-type UMAP coloured by R/NR; visually highlights expansion of Dysfunctional-effector CD8 in R.

![02a_celltype_response_umap](figures/02a_celltype_response_umap.png)
#### `figures/02a_composition_stacked_bar.png`

Per-patient stacked composition bars ordered by response group.

![02a_composition_stacked_bar](figures/02a_composition_stacked_bar.png)
#### `figures/02a_marker_dotplot.png`

Lineage-marker dot plot used for the annotation pass on GSE206325.

![02a_marker_dotplot](figures/02a_marker_dotplot.png)
#### `figures/02a_marker_feature_hexbin.png`

Hexbin marker feature plots — handles high cell density in the immune compartment.

![02a_marker_feature_hexbin](figures/02a_marker_feature_hexbin.png)
#### `figures/02a_signature_radar.png`

Radar plot of 10 IO signatures for R vs NR.

![02a_signature_radar](figures/02a_signature_radar.png)
#### `figures/02a_signature_violins.png`

Per-patient violins of cytotoxicity, exhaustion, IFN-γ, Tfh-like and TLS signatures.

![02a_signature_violins](figures/02a_signature_violins.png)
#### `figures/02a_tissue_signature_heatmap.png`

Tissue-level z-scored heatmap of 20 hallmark + IO signatures across R/NR samples.

![02a_tissue_signature_heatmap](figures/02a_tissue_signature_heatmap.png)
#### `figures/02a_volcano_deg.png`

Pseudo-bulk DE volcano (DESeq2 on patient pseudobulks; |log2FC|≥1, padj<0.05).

![02a_volcano_deg](figures/02a_volcano_deg.png)
#### `figures/02f_tf_regulon_activity.png`

L2 SCENIC-style TF regulon activity heatmap per patient; supports Extended Data Fig. 2.

![02f_tf_regulon_activity](figures/02f_tf_regulon_activity.png)
#### `figures/eda_deg_counts_barplot.png`

Exploratory DEG-count bars per cell type from the pseudo-bulk DE pass.

![eda_deg_counts_barplot](figures/eda_deg_counts_barplot.png)
#### `figures/eda_top_degs_heatmap.png`

Top-DEG heatmap from pseudo-bulk DE — exploratory companion to the volcano plot.

![eda_top_degs_heatmap](figures/eda_top_degs_heatmap.png)

### `figures/_archive/`
#### `figures/_archive/02_abundance_barplot.png`

Legacy duplicate of the per-cell-type abundance bars (kept for traceability).

![02_abundance_barplot](figures/_archive/02_abundance_barplot.png)
#### `figures/_archive/02_activation_exhaustion_scatter.png`

Per-patient scatter of activation vs exhaustion signature means — Responders cluster in a 'high-activation-with-productive-dysfunction' region.

![02_activation_exhaustion_scatter](figures/_archive/02_activation_exhaustion_scatter.png)
#### `figures/_archive/02_celltype_response_umap.png`

Legacy duplicate of the response-coloured cell-type UMAP.

![02_celltype_response_umap](figures/_archive/02_celltype_response_umap.png)
#### `figures/_archive/02_composition_stacked_bar.png`

Legacy duplicate of the stacked composition bar plot.

![02_composition_stacked_bar](figures/_archive/02_composition_stacked_bar.png)
#### `figures/_archive/02_signature_radar.png`

Legacy duplicate of the IO signature radar plot.

![02_signature_radar](figures/_archive/02_signature_radar.png)
#### `figures/_archive/02_signature_violins.png`

Legacy duplicate of the per-patient signature violins.

![02_signature_violins](figures/_archive/02_signature_violins.png)
#### `figures/_archive/02_tissue_signature_heatmap.png`

Legacy duplicate of the tissue-level signature heatmap.

![02_tissue_signature_heatmap](figures/_archive/02_tissue_signature_heatmap.png)
#### `figures/_archive/02_volcano_deg.png`

Legacy duplicate of the pseudo-bulk DE volcano.

![02_volcano_deg](figures/_archive/02_volcano_deg.png)
