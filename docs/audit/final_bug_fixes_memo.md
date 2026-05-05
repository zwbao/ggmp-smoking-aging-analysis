# Final Bug Fixes Memo — Peer-Review Revision Package

**Date:** 2026-05-04
**Auditor goal:** close out Bugs 3 and 4 from `cross_pipeline_audit.md`
(family/genus Python OLS smoke_binary contamination, and ASCVD smoker
mapping in the multi-mediator panel). Bugs 1 and 2 (transform mismatch,
quartile encoding) were closed previously.

Refit scripts:
- `revision-analysis/refit_family_genus_smoking_ols.py`
- `revision-analysis/refit_mediation_panel.py`

---

## Bug 3 — Family / genus smoking OLS contamination

### Root cause
`run_revision_analysis.py:run_taxonomic_robustness()` line 400 used
`meta["smoke_binary"].notna()` to filter the smoking analytic frame.
Because `clean_metadata` builds `smoke_binary` via
`np.where(..., np.where(..., np.nan))` — which coerces the third branch
to the literal string `'nan'` rather than NaN — the `.notna()` filter
is True for **all** rows. The 591 former_smoker / not_everyday samples
(570 with non-missing covariates) leak into the smoking frame with
`smoke_everyday=0`, contaminating the never-smoker baseline.

### Fix
One-line change (in the refit script):
```python
# was:
smoking_meta = meta[meta["smoke_binary"].notna() & ...]
# now:
smoking_meta = meta[
    meta["smk_status"].isin(["everyday", "never_smoker"])
    & meta["age"].notna() & meta["bmi"].notna() & meta["bristol"].notna()
]
```
n: **6,496 (buggy) → 5,926 (corrected)**.

### Outputs (all under `revision-analysis/outputs/`)
- `family_smoking_results_corrected.tsv` (replaces `family_smoking_results.tsv` for downstream use)
- `family_shared_results_corrected.tsv` (smoking-side recomputed; age-side joined from existing CLEAN `family_age_results.tsv`)
- `genus_smoking_results_corrected.tsv`
- `genus_shared_results_corrected.tsv`
- `supp_family_overlap_scatter.{png,pdf}` (overwrote, identical style)
- `supp_genus_overlap_scatter.{png,pdf}` (overwrote)
- Buggy backups preserved as `family_smoking_results_buggy_smoke_binary.tsv`,
  `family_shared_results_buggy_smoke_binary.tsv`,
  `genus_smoking_results_buggy_smoke_binary.tsv`,
  `genus_shared_results_buggy_smoke_binary.tsv`,
  `supp_family_overlap_scatter_buggy_smoke_binary.{png,pdf}`,
  `supp_genus_overlap_scatter_buggy_smoke_binary.{png,pdf}`.

### Before vs after — counts

| Metric                                       | Buggy | Corrected |
|----------------------------------------------|------:|----------:|
| n_samples (smoking model)                    | 6,496 | 5,926     |
| Family — features tested                     |    83 |    83     |
| Family — FDR-significant for smoking         |     8 |     7     |
| Family — shared & both-FDR-sig               |     4 |     3     |
| Family — shared & direction-concordant       |     1 |     0     |
| Genus — features tested                      |   154 |   154     |
| Genus — FDR-significant for smoking          |    15 |    15     |
| Genus — shared & both-FDR-sig                |     5 |     4     |
| Genus — shared & direction-concordant        |     1 |     0     |

The headline "**1 family / 1 genus concordant (Turicibacteraceae /
Turicibacter)**" **NO LONGER HOLDS** under the corrected filter.
Turicibacteraceae and Turicibacter both lose smoking-side FDR
significance (q_smoking ≈ 0.09 and 0.08, respectively), so neither is
in the both-FDR-sig set anymore, and the corrected panel has zero
direction-concordant taxa at the family or genus level.

### Named-taxa status table (corrected)

| Taxon                  | Level  | Smoking q | Age q     | direction_concordant | both_fdr_sig |
|------------------------|--------|----------:|----------:|---------------------:|-------------:|
| Coriobacteriaceae      | family | 4.12e-11  | 0.344     | False                | False        |
| Actinomycetaceae       | family | 8.54e-06  | 5.16e-05  | False                | True         |
| Erysipelotrichaceae    | family | 0.139     | 1.28e-04  | False                | False        |
| Turicibacteraceae      | family | 0.0927    | 0.0118    | True                 | **False**    |
| Atopobium              | genus  | 6.74e-26  | 3.77e-04  | False                | True         |
| Actinomyces            | genus  | 3.36e-06  | 8.03e-05  | False                | True         |
| Dorea                  | genus  | 5.27e-03  | 1.35e-04  | False                | True         |
| Turicibacter           | genus  | 0.0809    | 0.0133    | True                 | **False**    |

Coriobacteriaceae and Erysipelotrichaceae lose their family-level
both_fdr_sig status under the corrected filter (Coriobacteriaceae's age
side was already not FDR-sig; Erysipelotrichaceae's smoking q drops
from 0.0054 to 0.139).

### What still holds for the response letter
- The **MaAsLin2 cross-check (`*_smoking_results_maaslin2.tsv`,
  `maaslin2_summary.tsv`) is CLEAN** — `run_maaslin2_reanalysis.R`
  uses R's `ifelse` correctly, the smoking subset is the right
  n=5,926, and those numbers are what the response letter should
  primarily cite. The MaAsLin2 cross-framework concordance results
  reported under R1-2 are unaffected by Bug 3.
- The Python OLS supplementary TSVs are now CORRECT and the figure
  panels still convey the same general structure; only the
  Turicibacter direction-concordant story shifts.

---

## Bug 4 — ASCVD smoker mapping in the multi-mediator panel

### Root cause
`run_mediation_panel_originalspec.py` line 243 used
```python
df["smk_status"].map({"never_smoker": 0, "current_smoker": 1})
```
But GGMP `smk_status` only contains
`{everyday, former_smoker, not_everyday, never_smoker}` — no
`current_smoker`. So `everyday`, `former_smoker`, and `not_everyday`
all silently became NaN, the ASCVD computation returned NaN for those,
and the mediation regressions effectively ran on **never_smokers only**
(`n_complete = 3,426`).

### Fix
The published `220914/cvrisk.R` (line 19) filters
`smk_status %in% c("never_smoker", "current_smoker")` after applying
`code.R` line 24's mapping
`smk_status3 := str_replace_all(smk_status, c("not_everyday" = "current_smoker", "everyday" = "current_smoker"))`.
That definition leaves `former_smoker` unmapped, so it is implicitly
EXCLUDED from ASCVD. We adopt this exact mapping (which also matches
the PCE rationale that "current smoker = smoking now"):

| `smk_status`     | smoker → |
|------------------|---------:|
| never_smoker     |        0 |
| everyday         |        1 |
| not_everyday     |        1 |
| former_smoker    |   **NaN — excluded** |

This deviates from the brief's suggestion (`former_smoker → 0`); we
follow the published study because the brief explicitly said "unless
cvrisk.R reveals a different mapping that we must match (in which case
use that)".

### Outputs
- `mediation_panel_results_corrected.tsv`
- `mediation_panel_summary_corrected.md`
- `supp_mediation_panel.{png,pdf}` (overwritten)
- Buggy backups: `supp_mediation_panel_smoker_buggy.{png,pdf}` (the
  buggy TSV remains as `mediation_panel_results.tsv`).

### Before vs after — full mediator panel

| Mediator       | n_complete (was→is) | q_a (was→is)             | q_b (was→is)              | indirect (was→is)   | prop. mediated (was→is) |
|----------------|---------------------|--------------------------|---------------------------|---------------------|--------------------------|
| anthrop_SBP    | 3,426 → 4,642       | 6.18e-14 → 1.58e-20      | 2.45e-217 → 4.05e-253     | 2.073 → 2.351       | 49.5% → **39.7%**        |
| anthrop_DBP    | 3,426 → 4,642       | 2.46e-05 → 2.12e-09      | 2.52e-33 → 5.28e-55       | 0.485 → 0.743       | 11.6% → 12.5%            |
| anthrop_BMI    | 3,420 → 4,635       | 0.861 → 0.682            | 0.168 → 0.187             | 0.0022 → -0.0041    | NA → NA (sign disagree)  |
| anthrop_waist  | 3,420 → 4,635       | 0.0282 → 0.0202          | 9.44e-22 → 7.51e-40       | 0.201 → 0.247       | 4.79% → 4.16%            |
| biochem_FBG    | 3,426 → 4,642       | 2.54e-03 → 3.10e-03      | 5.71e-30 → 7.51e-40       | 0.336 → 0.325       | 8.02% → 5.49%            |
| biochem_TCHO   | 3,426 → 4,642       | 5.39e-07 → 3.41e-07      | 7.55e-05 → 1.92e-05       | 0.191 → 0.176       | 4.55% → 2.98%            |
| biochem_TG     | 3,426 → 4,642       | 0.650 → 0.121            | 1.34e-06 → 1.98e-32       | 0.0284 → 0.145      | 0.68% → 2.45%            |
| biochem_UA     | 3,426 → 4,642       | 0.702 → 0.0139           | 9.84e-53 → 1.84e-90       | -0.0653 → 0.401     | NA → 6.76%               |
| biochem_HDL    | 3,426 → 4,642       | 0.0282 → 0.0496          | 6.53e-15 → 1.24e-40       | -0.165 → -0.210     | NA → NA (sign disagree)  |
| biochem_LDL    | 3,426 → 4,642       | 4.05e-03 → 1.39e-02      | 2.96e-05 → 2.35e-02       | 0.117 → 0.0473      | 2.80% → 0.80%            |

(95% CIs for indirect_ab and prop_mediated for the corrected panel are
in `mediation_panel_results_corrected.tsv`. Notable corrected CIs:
SBP indirect 95% CI [1.85, 2.87]; SBP prop. mediated 95% CI [32.4%,
48.7%]; DBP indirect 95% CI [0.50, 1.00]; UA indirect 95% CI [0.091,
0.699].)

### Headline counts

| Metric                                      | Buggy | Corrected |
|---------------------------------------------|------:|----------:|
| n_complete (per mediator, typical)          | 3,426 |    4,642  |
| Mediators jointly FDR-sig (both q_a, q_b)   |     7 |        8  |
| SBP rank by \|indirect\|                    |    #1 |       #1  |
| SBP proportion mediated                     | 49.5% |    39.7%  |
| BMI status                                  | path-a NS | path-a NS (consistent) |

**New mediators that gained joint-FDR significance under the
correction:** biochem_UA (q_a 0.702 → 0.014; was not jointly sig
before).

**Mediators that lost joint-FDR significance:** none — biochem_TG was
NOT jointly sig before (q_a=0.65) and remains NOT (q_a=0.12).

**Direction flips:** biochem_UA's indirect flipped sign from -0.065
(previously path-a NS, sign noise) to +0.401 (now jointly sig,
positive direction).

### Qualitative narrative — does R2-6 survive?
Yes. SBP retains its #1 rank by absolute indirect effect among jointly
FDR-sig mediators, and its proportion-mediated remains the dominant
pathway (39.7% > everything else). The R2-6 conclusion that "SBP is
the dominant cardiometabolic mediator of pro-aging-microbiome →
ASCVD" is qualitatively preserved. **However the headline number
49.5% must be updated to 39.7% (95% CI 32.4–48.7%)**, and the
"7 of 10 jointly FDR-sig" claim must become "8 of 10".

The published Fig. 6b reported 28.5% for SBP. Under the buggy panel
this looked like 49.5% (delta +21 pp, MISMATCHED); under the corrected
panel it is 39.7% (delta +11 pp, still MISMATCHED but closer). The
gap from 28.5% likely reflects a different exposure / covariate spec
in the published primary analysis (e.g., not the multi-OTU module
score) — this is a documented limitation that the response letter
already acknowledges.

---

## Recommendation for response letter

### R1-2 (family / genus robustness paragraph)

**Numerical update needed (not a narrative rewrite).** Replace any
language that says "the Python OLS robustness check at the family /
genus level recapitulated 1 direction-concordant family
(Turicibacteraceae) and 1 direction-concordant genus (Turicibacter)
that were jointly FDR-significant for both age (within never-smokers)
and smoking" with language that:

1. Cites **MaAsLin2** as the primary cross-framework check (those
   numbers are clean and unaffected); and
2. Notes that under the corrected Python OLS smoking model
   (n=5,926 everyday-vs-never), **0 family and 0 genus features are
   simultaneously FDR-significant on both panels AND
   direction-concordant**. The smoking signal at the family level is
   now Actinomycetaceae (q=8.5e-06, opposite-sign with age),
   Dethiosulfovibrionaceae, Pseudonocardiaceae; at the genus level,
   Atopobium, Actinomyces, Dorea, Pyramidobacter (all opposite-sign
   with age). These are biologically-consistent: smoking-induced and
   age-associated taxa often shift in opposite directions in this
   cohort.

**Suggested replacement text** (drop into R1-2):

> We re-ran the family- and genus-level robustness analysis with the
> sample-filter encoding bug corrected (smk_status ∈ {everyday,
> never_smoker}, n = 5,926 vs the previously reported n = 6,496;
> outputs `family_smoking_results_corrected.tsv`,
> `genus_smoking_results_corrected.tsv`). Family-level: 7 of 83
> features were FDR-significant for smoking (was 8); 3 of those were
> also FDR-significant for age within never-smokers (was 4); none
> showed direction-concordant signs (was 1, Turicibacteraceae, which
> drops to q_smoking = 0.09 in the corrected fit). Genus-level: 15 of
> 154 features were FDR-significant for smoking (unchanged); 4 were
> jointly FDR-significant with age (was 5); none direction-concordant
> (was 1, Turicibacter, q_smoking = 0.08 corrected). The MaAsLin2
> cross-framework reanalysis (`run_maaslin2_reanalysis.R`,
> n = 5,926) was unaffected by this bug and remains the primary
> robustness reference.

### R2-6 (multi-mediator panel)

**Numerical update needed (narrative survives).** Update the headline
numbers but keep the conclusion.

**Suggested replacement text** (drop into R2-6):

> Among the ten candidate cardiometabolic mediators systematically
> tested under the original `gai_med.R` specification (no covariate
> adjustment beyond the exposure and mediator), 8 of 10 reached
> FDR-significance on both the path-a and path-b regressions
> (n_complete = 4,642). Systolic blood pressure (SBP) ranked #1 by
> absolute indirect effect with a proportion mediated of 39.7%
> (95% CI 32.4–48.7%), supporting its inclusion as the primary
> mediator in Figure 6b. The other 7 jointly FDR-sig mediators were
> DBP (12.5%), uric acid (6.8%), fasting glucose (5.5%), waist
> circumference (4.2%), total cholesterol (3.0%), LDL (0.8%), and
> HDL (sign-flip; proportion undefined under ACME/(ACME+ADE)). BMI
> remained non-significant on path a (q = 0.68) under this spec,
> consistent with the previous panel. The corrected ASCVD smoker
> indicator follows the published `cvrisk.R` definition (everyday +
> not_everyday → current = 1; never_smoker → 0; former_smoker
> excluded), which fixes the previous mapping that silently dropped
> all current-smokers from ASCVD because GGMP `smk_status` does not
> contain the literal value "current_smoker".

**Note:** the previously reported 49.5% for SBP and "7 of 10 jointly
FDR-sig" came from the buggy panel, which ran on 3,426 never-smokers
only. These numbers should NOT be cited.

### R2-2 (male-only quartile sensitivity)
Unchanged by Bug 3 / Bug 4. The R2-2 paragraph already uses the
Approach-B numbers from `male_only_module_effects_by_age_corrected_4level.tsv`
which were fixed in the previous round.

---

## Files written / overwritten

### Refit scripts
- `revision-analysis/refit_family_genus_smoking_ols.py` (new)
- `revision-analysis/refit_mediation_panel.py` (new)

### TSVs (new)
- `outputs/family_smoking_results_corrected.tsv`
- `outputs/family_shared_results_corrected.tsv`
- `outputs/genus_smoking_results_corrected.tsv`
- `outputs/genus_shared_results_corrected.tsv`
- `outputs/mediation_panel_results_corrected.tsv`

### Markdown summary (new)
- `outputs/mediation_panel_summary_corrected.md`
- `outputs/final_bug_fixes_memo.md` (this file)

### Figures (overwritten — buggy backups preserved)
- `outputs/supp_family_overlap_scatter.{png,pdf}` (buggy → `*_buggy_smoke_binary.{png,pdf}`)
- `outputs/supp_genus_overlap_scatter.{png,pdf}` (buggy → `*_buggy_smoke_binary.{png,pdf}`)
- `outputs/supp_mediation_panel.{png,pdf}` (buggy → `*_smoker_buggy.{png,pdf}`)

### Buggy TSV backups
- `outputs/family_smoking_results_buggy_smoke_binary.tsv`
- `outputs/family_shared_results_buggy_smoke_binary.tsv`
- `outputs/genus_smoking_results_buggy_smoke_binary.tsv`
- `outputs/genus_shared_results_buggy_smoke_binary.tsv`
- (`outputs/mediation_panel_results.tsv` is the buggy version; the
  corrected one lives at `mediation_panel_results_corrected.tsv`.)
