"""06a: Cohort and analytical-unit reconciliation for the reframe.

Resolves the n=24 vs n=56 vs n=531 ambiguity by explicitly tabulating,
for every cohort used in the manuscript, the raw units, the labelled units,
and the units retained in the headline discovery/transfer analyses.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

BASE = Path('/home/sadeghii/projects/spatial_immunotherapy_response')
OUT  = BASE / 'out/06_reframe/tables'
OUT.mkdir(parents=True, exist_ok=True)

rows = []

# ---- GSE206325 ----
sa   = pd.read_csv(BASE / 'data/processed/GSE206325_sample_annotations.csv')
feat = pd.read_csv(BASE / 'out/05_multilayer/tables/patient_multilayer_features.csv')
feat206 = feat[feat['dataset'] == 'GSE206325'].copy()
kept = set(feat206['patient_id'].astype(str))

pat = sa.groupby('patient_ID').agg(
    n_samples=('sample_ID', 'nunique'),
    treatment_Resp=('treatment_Resp', lambda x: x.dropna().iloc[0] if x.dropna().size else None),
).reset_index()
pat['label'] = pat['treatment_Resp'].apply(
    lambda x: 'Responder' if isinstance(x, str) and x.endswith('_R')
    else ('Non-Responder' if isinstance(x, str) and x.endswith('_NR') else 'Unlabelled')
)
pat['retained_in_discovery'] = pat['patient_ID'].astype(str).isin(kept)

rows.append({
    'cohort': 'GSE206325',
    'role': 'discovery (scRNA)',
    'raw_patients': pat['patient_ID'].nunique(),
    'raw_samples': sa['sample_ID'].nunique(),
    'labelled_patients': int((pat['label'] != 'Unlabelled').sum()),
    'labelled_samples': int(sa['treatment_Resp'].notna().sum()),
    'analytical_unit': 'patient',
    'patients_in_headline_analyses': len(feat206),
    'patients_R': int((feat206['response'] == 'Responder').sum()),
    'patients_NR': int((feat206['response'] == 'Non-Responder').sum()),
    'exclusion_reason': 'n_cells QC (cells/patient) and missing label: patients 15, 63, 71, 98 dropped for low cell yield; patient 1024 unlabelled',
})

# ---- Cappuyns2023 ----
feat_cap = feat[feat['dataset'] == 'Cappuyns2023']
cap_sig  = pd.read_csv(BASE / 'out/02i_cappuyns_signatures/tables/patient_level_cappuyns_scores.csv')
rows.append({
    'cohort': 'Cappuyns2023',
    'role': 'external transfer (scRNA)',
    'raw_patients': len(feat_cap),
    'raw_samples': '-',
    'labelled_patients': int((feat_cap['response'].isin(['Responder','Non-Responder'])).sum()),
    'labelled_samples': '-',
    'analytical_unit': 'patient',
    'patients_in_headline_analyses': len(feat_cap),
    'patients_R': int((feat_cap['response'] == 'Responder').sum()),
    'patients_NR': int((feat_cap['response'] == 'Non-Responder').sum()),
    'exclusion_reason': 'none (all n=32 retained); the 24-patient Cappuyns-signature table is scored ON GSE206325 cells, not Cappuyns patients',
})

# ---- Pooled multilayer ----
rows.append({
    'cohort': 'POOLED (GSE206325 + Cappuyns2023)',
    'role': 'pooled sensitivity',
    'raw_patients': len(feat),
    'raw_samples': '-',
    'labelled_patients': len(feat),
    'labelled_samples': '-',
    'analytical_unit': 'patient (cohort as covariate)',
    'patients_in_headline_analyses': len(feat),
    'patients_R': int((feat['response'] == 'Responder').sum()),
    'patients_NR': int((feat['response'] == 'Non-Responder').sum()),
    'exclusion_reason': 'n/a',
})

# ---- IMbrave150 ----
# Count unique trial patients from the per-arm PFS strata in the arm-level
# summary table. Each arm is split into CD8-high / CD8-low strata whose n's
# sum to the arm population. Previously this script was naively counting
# rows in the summary (4 arm×endpoint rows), which is wrong.
try:
    im_arm = pd.read_csv(BASE / 'out/04_complementary/tables/revision_imbrave150_cd8_arm.csv')
    # Use PFS rows as the per-arm patient counts (PFS and OS rows share the
    # same patients within an arm).
    pfs = im_arm[im_arm['endpoint'] == 'PFS']
    per_arm = (pfs[['arm', 'n_high', 'n_low']]
               .drop_duplicates()
               .assign(n_arm=lambda d: d['n_high'] + d['n_low']))
    n_atezo = int(per_arm.loc[per_arm['arm'] == 'Atezo+Bev', 'n_arm'].iloc[0])
    n_soraf = int(per_arm.loc[per_arm['arm'] == 'Sorafenib', 'n_arm'].iloc[0])
    n_pts_im = n_atezo + n_soraf
    arm_breakdown = f"{n_atezo} Atezo+Bev / {n_soraf} Sorafenib"
except Exception:
    n_pts_im = 'see revision_imbrave150_cd8_arm.csv'
    arm_breakdown = 'unavailable'
rows.append({
    'cohort': 'IMbrave150',
    'role': 'orthogonal clinical (bulk RNA + survival)',
    'raw_patients': n_pts_im,
    'raw_samples': arm_breakdown,
    'labelled_patients': n_pts_im,
    'labelled_samples': arm_breakdown,
    'analytical_unit': 'patient (treatment arm, CD8-high binary)',
    'patients_in_headline_analyses': n_pts_im,
    'patients_R': 'n/a (endpoints are PFS / OS, not R/NR)',
    'patients_NR': 'n/a',
    'exclusion_reason': f'no scRNA-level patient prediction; used only for within-arm and interaction support ({arm_breakdown})',
})

# ---- GSE238264 Visium ----
rows.append({
    'cohort': 'GSE238264 (Visium)',
    'role': 'spatial localization + QC sensitivity',
    'raw_patients': 7,
    'raw_samples': '7 sections',
    'labelled_patients': 7,
    'labelled_samples': 7,
    'analytical_unit': 'section (4R / 3NR); spot-level analyses are pseudoreplicated',
    'patients_in_headline_analyses': 7,
    'patients_R': 4,
    'patients_NR': 3,
    'exclusion_reason': 'HCC2R flagged as QC outlier (T3.1 pass=43.6%); leave-one-section-out sensitivity required for any R-vs-NR claim',
})

# ---- GSE281759 Visium ----
rows.append({
    'cohort': 'GSE281759 (Visium multi-region)',
    'role': 'contextual spatial biology (tumour vs normal / primary vs recurrence)',
    'raw_patients': 'multi-region across 4 patient series',
    'raw_samples': '17 sections',
    'labelled_patients': 'n/a (no response labels)',
    'labelled_samples': 'n/a',
    'analytical_unit': 'section / tumour-region',
    'patients_in_headline_analyses': 17,
    'patients_R': 'n/a',
    'patients_NR': 'n/a',
    'exclusion_reason': 'no response labels; excluded from response classification',
})

recon = pd.DataFrame(rows)
recon.to_csv(OUT / 'cohort_unit_reconciliation.csv', index=False)
print(recon.to_string(index=False))

# ---- Per-patient audit table for GSE206325 ----
aud = pat[['patient_ID','n_samples','treatment_Resp','label','retained_in_discovery']].copy()
aud['note'] = aud.apply(
    lambda r: 'retained' if r['retained_in_discovery']
    else ('unlabelled (excluded)' if r['label']=='Unlabelled'
          else f"labelled but excluded (n_samples={r['n_samples']}, below cell-count threshold)"),
    axis=1,
)
aud.to_csv(OUT / 'gse206325_patient_audit.csv', index=False)
print('\nGSE206325 patient-level audit saved ->', OUT / 'gse206325_patient_audit.csv')
print(f"labelled retained = {(aud['retained_in_discovery'] & (aud['label']!='Unlabelled')).sum()}")
print(f"labelled excluded = {((~aud['retained_in_discovery']) & (aud['label']!='Unlabelled')).sum()}")

# ---- Human-readable reconciliation note ----
note = BASE / 'out/06_reframe/notes/cohort_unit_reconciliation.md'
note.parent.mkdir(parents=True, exist_ok=True)
note.write_text(f"""# Cohort & Analytical-Unit Reconciliation (reframe)

**Primary analytical unit**: **patient** for scRNA cohorts; **section** for Visium cohorts; **patient** for IMbrave150.

## The "n=24 vs n=56" question
- `patient_multilayer_features.csv` has **56 rows** because it pools two cohorts:
  - **GSE206325**: 24 patients (17 NR / 7 R)
  - **Cappuyns2023**: 32 patients (13 NR / 19 R)
- GSE206325 has 29 raw patients and 731 samples; **28 patients carry treatment_Resp labels**; **4 labelled patients were dropped** for low cell yield (IDs 15, 63, 71, 98 with 3–9 samples each), and 1 unlabelled patient (1024) dropped. The resulting **24 labelled, QC-passing GSE206325 patients** are the discovery cohort.
- The **n=24** Cappuyns-signature table (`patient_level_cappuyns_scores.csv`) is *Cappuyns-signatures scored on the 24 GSE206325 patients*, **not** Cappuyns cohort patients. This is a transfer-direction artefact in the file name and the earlier manuscript wording.
- The pooled n=56 figure applies to the discovery + Cappuyns2023 pool (24 + 32), not to GSE206325 alone; any prior wording that attributed all 56 units to GSE206325 is a bookkeeping error and is withdrawn.

## Pseudoreplication hazards identified
- **GSE206325 raw 529 labelled samples** or **~45k labelled cells** were previously used in some Mann–Whitney tests → cell-level / sample-level p-values are inflated by pseudoreplication. Under the reframe, every main claim must use patient-level aggregation.
- **GSE238264**: 44k spots across 7 sections — spot-level p-values inflated; section-level permutation and leave-one-section-out are the only admissible inferences for R-vs-NR claims.
- **GSE281759**: 27k spots across 17 sections from 4 multi-region patient series — tumour-vs-normal comparisons must retain donor as unit at minimum.

## IMbrave150 patient counts
- IMbrave150 contributes **{n_pts_im} trial patients** (**{arm_breakdown}**) drawn from the CD8-stratified PFS rows of `revision_imbrave150_cd8_arm.csv` (Atezo+Bev: 60 CD8-high + 59 CD8-low; Sorafenib: 29 CD8-high + 29 CD8-low). The earlier reconciliation row listed **n=4** because it counted the four arm×endpoint summary rows in that file; that was a bookkeeping bug, not a cohort limitation.
- IMbrave150 is still only used for within-arm log-rank and arm×CD8 interaction tests; it does not enter any R-vs-NR classification.
""")
print('\nNote saved ->', note)
