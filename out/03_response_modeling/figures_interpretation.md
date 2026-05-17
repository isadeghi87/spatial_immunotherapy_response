# Figures interpretation — `out/03_response_modeling`

> Companion to `docs/manuscript.md` (May 2026 revision). All R-vs-NR
> statistics, AUCs and effect sizes referenced below are computed at the
> patient level on the GSE206325 anti-PD-1 discovery cohort (n=24; 7 R,
> 17 NR) unless otherwise stated. External validation uses the Cappuyns
> et al. 2023 atezolizumab + bevacizumab cohort (n=32 evaluable; 19 R,
> 13 NR). Statistics use two-sided Mann–Whitney U with BH-FDR within
> family, Cliff's δ for effect size, and 2,000-iteration patient-level
> bootstrap 95 % CIs.

## Scope
Logistic-regression / LASSO hyper-parameter sensitivity for the multi-layer predictor.

## Figure-by-figure interpretation

### `figures/`
#### `figures/05m_lasso_C_sensitivity.png`

LASSO penalty (C) sensitivity scan for the response classifier — selects C=0.5 for the 15-feature ensemble.

![05m_lasso_C_sensitivity](figures/05m_lasso_C_sensitivity.png)
