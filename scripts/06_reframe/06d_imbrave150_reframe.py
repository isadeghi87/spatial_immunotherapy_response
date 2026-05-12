"""06d: IMbrave150 orthogonal-support reframe.

The earlier wording treated IMbrave150 as "validation" of a CD8 response
biomarker. Under the reframe this is downgraded to *within-arm, orthogonal*
support for an immune-active benefit signal — with explicit statement that
the arm × CD8-high interaction is NOT robust (p=0.18 for PFS, p=0.07 for OS).

This script does NO new fitting (the Cox model table is already accepted
upstream); it writes a condensed, claim-bounded summary table and a short
note so the results_index entry can reflect the strict boundary.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

BASE = Path('/home/sadeghii/projects/spatial_immunotherapy_response')
OUT  = BASE / 'out/06_reframe/tables'
OUT.mkdir(parents=True, exist_ok=True)

arm = pd.read_csv(BASE / 'out/04_complementary/tables/revision_imbrave150_cd8_arm.csv')
cox = pd.read_csv(BASE / 'out/04_complementary/tables/revision_imbrave150_cox_interaction.csv')

summary = []
# Within-arm benefit (atezo+bev PFS by CD8-high): p=0.009
row = arm[(arm['arm']=='Atezo+Bev') & (arm['endpoint']=='PFS')].iloc[0]
summary.append(dict(
    claim='Within-arm PFS benefit of CD8-high stratum (Atezo+Bev)',
    endpoint='PFS',
    arm='Atezo+Bev',
    n_high=int(row['n_high']), n_low=int(row['n_low']),
    median_high_mo=float(row['median_high']), median_low_mo=float(row['median_low']),
    logrank_p=float(row['logrank_p']),
    strength_class='supported_limited',
    finding='CD8-high patients on Atezo+Bev have longer median PFS (8.4 vs 5.6 mo, log-rank p=0.009).',
    scope='Within-arm descriptive contrast only; does not license a differential-treatment-benefit claim for CD8-high vs CD8-low — see arm × CD8 interaction row.',
))
row = arm[(arm['arm']=='Sorafenib') & (arm['endpoint']=='PFS')].iloc[0]
summary.append(dict(
    claim='Within-arm PFS of CD8-high stratum (Sorafenib control)',
    endpoint='PFS', arm='Sorafenib',
    n_high=int(row['n_high']), n_low=int(row['n_low']),
    median_high_mo=float(row['median_high']), median_low_mo=float(row['median_low']),
    logrank_p=float(row['logrank_p']),
    strength_class='null',
    finding='In the sorafenib control arm, CD8-high shows no PFS benefit (log-rank p=0.90).',
    scope='',
))
# Arm × CD8 interaction — the formal differential-benefit test
row_i = cox[(cox['endpoint']=='PFS') & (cox['covariate']=='arm_x_cd8')].iloc[0]
summary.append(dict(
    claim='Arm × CD8-high interaction (differential-benefit test), PFS',
    endpoint='PFS', arm='interaction',
    n_high='',n_low='',median_high_mo='',median_low_mo='',logrank_p='',
    HR=float(row_i['hr']), HR_CI_lo=float(row_i['hr_lo']), HR_CI_hi=float(row_i['hr_hi']),
    cox_p=float(row_i['p']),
    strength_class='null',
    finding='The arm × CD8 interaction is not statistically significant for PFS (HR 0.57, 95% CI 0.25–1.30, p=0.18).',
    scope='No evidence of differential PFS benefit by CD8 stratum between arms.',
))
row_i = cox[(cox['endpoint']=='OS') & (cox['covariate']=='arm_x_cd8')].iloc[0]
summary.append(dict(
    claim='Arm × CD8-high interaction, OS',
    endpoint='OS', arm='interaction',
    n_high='',n_low='',median_high_mo='',median_low_mo='',logrank_p='',
    HR=float(row_i['hr']), HR_CI_lo=float(row_i['hr_lo']), HR_CI_hi=float(row_i['hr_hi']),
    cox_p=float(row_i['p']),
    strength_class='null',
    finding='The arm × CD8 interaction for OS is suggestive but not conventionally significant (HR 0.36, 95% CI 0.12–1.07, p=0.065).',
    scope='Hypothesis-generating only; not a differential-benefit claim.',
))

out = pd.DataFrame(summary)
out.to_csv(OUT / 'imbrave150_orthogonal_support.csv', index=False)
print(out[['claim','endpoint','strength_class','finding']].to_string(index=False))
print('\nWritten:', OUT / 'imbrave150_orthogonal_support.csv')
