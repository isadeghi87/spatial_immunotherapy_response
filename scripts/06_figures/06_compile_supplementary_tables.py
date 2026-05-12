#!/usr/bin/env python3
"""
06_compile_supplementary_tables.py — Consolidate all CSV tables into a
single multi-sheet Excel workbook for journal submission.

Produces: out/06_manuscript_figures/Supplementary_Tables.xlsx
"""
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT / "out" / "06_manuscript_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Table manifest: (sheet_name, csv_path, description) ───────────────
TABLES = [
    # Preprocessing & QC
    ("S1_QC_Summary",
     "out/01_preprocessing/tables/qc_summary.csv",
     "Quality-control metrics per dataset"),
    ("S2_Response_Counts",
     "out/01_preprocessing/tables/response_counts_by_dataset.csv",
     "Responder / Non-Responder patient counts per cohort"),

    # TME characterisation
    ("S3_CellType_Abundance",
     "out/02a_response_tme/tables/celltype_abundance_R_vs_NR.csv",
     "Cell-type abundance comparison (R vs NR, Fisher exact)"),
    ("S4_DEGs_Per_CellType",
     "out/02a_response_tme/tables/deg_R_vs_NR_per_celltype.csv",
     "Differentially expressed genes per cell type (R vs NR)"),
    ("S5_IO_Signatures_RvsNR",
     "out/02a_response_tme/tables/signature_R_vs_NR.csv",
     "IO signature scores comparison (R vs NR)"),
    ("S6_Tissue_Signatures",
     "out/02a_response_tme/tables/tissue_signature_R_vs_NR.csv",
     "Tissue-stratified signature comparisons"),

    # Spatial biomarkers
    ("S7_Visium_CellType_Comp",
     "out/02b_spatial_biomarkers/tables/gse238264_celltype_composition.csv",
     "GSE238264 Visium cell-type deconvolution per section"),
    ("S8_Morans_I",
     "out/02b_spatial_biomarkers/tables/gse238264_morans_I.csv",
     "Moran's I spatial autocorrelation statistics"),
    ("S9_SVGs",
     "out/02b_spatial_biomarkers/tables/gse238264_top_SVGs_per_section.csv",
     "Top spatially variable genes per Visium section"),
    ("S10_LigandReceptor_Spatial",
     "out/02b_spatial_biomarkers/tables/gse238264_ligand_receptor.csv",
     "Ligand-receptor spatial co-localisation scores"),
    ("S11_Spatial_Entropy",
     "out/02b_spatial_biomarkers/tables/gse238264_spatial_entropy.csv",
     "Spatial Shannon entropy per section"),
    ("S12_Permutation_Tests",
     "out/02b_spatial_biomarkers/tables/gse238264_sample_level_permutation.csv",
     "Sample-level signature permutation test results"),
    ("S13_LOO_Sensitivity",
     "out/02b_spatial_biomarkers/tables/gse238264_leave_one_out.csv",
     "Leave-one-out sensitivity analysis"),
    ("S14_Pseudobulk_Comp",
     "out/02b_spatial_biomarkers/tables/pseudobulk_group_comparison.csv",
     "Pseudo-bulk group comparisons"),

    # Cappuyns signatures
    ("S15_Cappuyns_Stats",
     "out/02i_cappuyns_signatures/tables/cappuyns_signature_stats.csv",
     "Cappuyns signature differential statistics (R vs NR, BH-corrected)"),
    ("S16_Cappuyns_Patient",
     "out/02i_cappuyns_signatures/tables/patient_level_cappuyns_scores.csv",
     "Patient-level Cappuyns signature scores"),
    ("S17_Molecular_Subsets",
     "out/02i_cappuyns_signatures/tables/patient_molecular_subsets.csv",
     "Patient molecular subset classification (Cappuyns)"),

    # Concordance
    ("S18_Concordance_Stats",
     "out/02j_concordance/tables/concordance_stats.csv",
     "Chi-squared, Cohen's kappa for internal vs Cappuyns concordance"),
    ("S19_Concordance_Crosstab",
     "out/02j_concordance/tables/concordance_crosstab.csv",
     "Cross-tabulation: internal subtypes vs Cappuyns subsets"),

    # Cappuyns spatial
    ("S20_Cappuyns_Spatial",
     "out/02k_cappuyns_spatial/tables/spatial_cappuyns_section_scores.csv",
     "Cappuyns signature spatial scores per Visium section"),

    # Validation
    ("S21_LASSO_Coefficients",
     "out/03_validation/tables/lasso_coefficients.csv",
     "LASSO logistic regression coefficients (LOOCV)"),
    ("S22_RF_Importance",
     "out/03_validation/tables/rf_feature_importance.csv",
     "Random Forest feature importances"),
    ("S23_Patient_Features",
     "out/03_validation/tables/patient_features.csv",
     "Patient-level feature matrix used for modelling"),
    ("S24_Bootstrap_Effects",
     "out/03_validation/tables/bootstrap_effects.csv",
     "Bootstrap confidence intervals for signature effects"),
    ("S25_Bulk_Signatures",
     "out/03_validation/tables/bulk_signature_scores.csv",
     "Bulk RNA-seq signature scores (GSE202069)"),
    ("S26_Signature_Corr",
     "out/03_validation/tables/signature_correlation.csv",
     "Signature-signature Spearman correlation network"),
    ("S27_Pathway_Enrichment",
     "out/03_validation/tables/pathway_enrichment.csv",
     "Pathway enrichment on R vs NR DEGs"),
    ("S28_TCR_Clonality",
     "out/03_validation/tables/tcr_clonality.csv",
     "TCR clonality metrics per patient"),

    # Complementary
    ("S29_CD8_State_Fractions",
     "out/04_complementary/tables/cd8_state_fractions.csv",
     "CD8 T-cell state fractions per patient"),
    ("S30_CD8_State_Tests",
     "out/04_complementary/tables/cd8_state_tests.csv",
     "Statistical tests for CD8 state differences (R vs NR)"),
    ("S31_Hallmark_Pathways",
     "out/04_complementary/tables/hallmark_pathway_R_vs_NR.csv",
     "Hallmark pathway enrichment (R vs NR)"),
    ("S32_Immune_Ratios",
     "out/04_complementary/tables/immune_ratios.csv",
     "Immune cell ratios per patient (CD8/Treg, M1/M2, etc.)"),
    ("S33_Immune_Ratio_AUCs",
     "out/04_complementary/tables/immune_ratio_aucs.csv",
     "AUC for immune ratio biomarkers"),
    ("S34_Patient_Subtypes",
     "out/04_complementary/tables/patient_immune_subtypes.csv",
     "Patient immune subtype assignments"),
    ("S35_Subtype_Response",
     "out/04_complementary/tables/subtype_response_rates.csv",
     "Response rates per immune subtype"),
    ("S36_LR_Interactions",
     "out/04_complementary/tables/lr_interaction_scores.csv",
     "Ligand-receptor interaction scores"),
    ("S37_Cross_Dataset_Sigs",
     "out/04_complementary/tables/cross_dataset_signatures.csv",
     "Cross-dataset signature reproducibility"),
]

# ── Build workbook ─────────────────────────────────────────────────────
def main():
    out_path = OUT_DIR / "Supplementary_Tables.xlsx"
    print("=" * 70)
    print("Compiling Supplementary Tables")
    print("=" * 70)

    # Create a table-of-contents DataFrame
    toc_rows = []

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, csv_rel, desc in TABLES:
            csv_path = PROJECT / csv_rel
            if not csv_path.exists():
                print(f"  SKIP: {csv_rel} (not found)")
                continue
            try:
                df = pd.read_csv(csv_path)
            except Exception as e:
                print(f"  SKIP: {csv_rel} ({e})")
                continue

            # Excel sheet names max 31 chars
            safe_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=safe_name, index=False)
            toc_rows.append({
                "Sheet": safe_name,
                "Description": desc,
                "Source CSV": csv_rel,
                "Rows": len(df),
                "Columns": len(df.columns),
            })
            print(f"  {safe_name}: {len(df)} rows × {len(df.columns)} cols")

        # Write table-of-contents as first sheet
        if toc_rows:
            toc_df = pd.DataFrame(toc_rows)
            toc_df.to_excel(writer, sheet_name="Table_of_Contents", index=False)
            # Move TOC to first position
            wb = writer.book
            ws_toc = wb["Table_of_Contents"]
            wb.move_sheet(ws_toc, offset=-(len(wb.sheetnames) - 1))

    print(f"\n  Saved: {out_path}")
    print(f"  Sheets: {len(toc_rows)} data tables + Table of Contents")
    print("=" * 70)


if __name__ == "__main__":
    main()
