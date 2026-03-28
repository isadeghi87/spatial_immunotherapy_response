# Spatial Immunotherapy Response — Analysis Prompts

## Prompt 1: IO-HCC Data Assembly
1. Download all publicly available HCC immunotherapy scRNA-seq datasets (see README for full list)
2. Harmonize response labels: RECIST CR+PR = Responder; SD+PD = Non-Responder
3. Record treatment type, line, duration, and clinical features
4. Integrate using scVI with study as batch, preserving response-associated biological variation
5. Create unified cell type annotations focusing on immune subsets

## Prompt 2: Response-Associated TME Characterization
1. Differential cell type abundance: Responders vs Non-Responders (Dirichlet regression)
2. Within each cell type, identify DEGs between R and NR (pseudobulk DESeq2)
3. Focus analyses on:
   - CD8 T cells: exhaustion depth (TOX, TIGIT, LAG3, HAVCR2), stem-like markers (TCF7, SLAMF6)
   - Macrophages: M1/M2 ratio, SPP1+ vs C1QA+ TAM balance
   - Tregs: frequency, CTLA4/TIGIT/ICOS expression
   - DCs: cDC1 frequency (BATF3+CLEC9A+) as cross-presentation capacity
4. Build response signature per cell type

## Prompt 3: Spatial Biomarker Discovery
1. Apply spatial immunotyping to all available spatial HCC data
2. Extract spatial features: EIES, TLS, CD8-tumor distance, border regularity
3. Correlate spatial features with response (where data available)
4. Build LASSO model: spatial features → response prediction
5. Compare with bulk-derived biomarkers (TIDE, TMB, PD-L1)

## Prompt 4: Deconvolution of Bulk IO Trials
1. Deconvolve IMbrave150/CheckMate-459 bulk RNA-seq using CIBERSORTx with HCC scRNA-seq reference
2. Infer spatial immunotype from deconvolved proportions
3. Test association: inferred immunotype vs. objective response rate
4. Survival analysis: PFS and OS by inferred immunotype
5. Establish that spatial architecture biomarkers outperform abundance-only biomarkers
