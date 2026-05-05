# Cross-Pipeline Audit — Revision Analysis Package

**Audit date:** 2026-05-04
**Auditor goal:** find ALL transform-mismatch and categorical-encoding /
sample-filter bugs in the supplementary analyses of the peer-review revision,
beyond the two already documented in
`outputs/male_only_otu_concordance_fix_memo.md`.

The two known bugs:

- **Bug 1 (transform mismatch — FIXED):** male-only OTU-level smoking
  coefficients in `male_only_overlap_otu_concordance.tsv` were computed in
  Python OLS arcsin-sqrt space, then plotted against full-sample MaAsLin2
  LOG coefficients. Fix: `run_maaslin2_male_only_otu.R` re-runs the male
  model in the same MaAsLin2 LOG framework; new TSV
  `male_only_overlap_otu_concordance_maaslin2.tsv`.
- **Bug 2 (categorical encoding bug):** `clean_metadata` line 60 builds
  `smoke_binary` via
  `np.where(eq("everyday"), "everyday", np.where(eq("never_smoker"), "never_smoker", np.nan))`.
  Because the inner `np.where` returns an `object` array of strings, the
  third branch becomes the **literal string `'nan'`**, not `NaN`.
  Downstream `meta["smoke_binary"].notna()` filters then let those
  former_smoker / not_everyday rows through, with
  `smoke_everyday = (smoke_binary == "everyday").astype(float) = 0.0`. The
  baseline is therefore "everyday vs (never + former + not_everyday)",
  not the intended "everyday vs never".

I confirmed empirically:
```
Unique values in buggy smoke_binary: ['nan' 'never_smoker' 'everyday']
.notna() count: 6676 (entire frame!)
```

This audit traces every consumer of `clean_metadata` and every
cross-framework coefficient comparison.

---

## File 1 — `run_revision_analysis.py`

### Finding 1.1 — Quartile-level male module forest (left panel of `supp_male_only_sensitivity.png`)

- **Severity: CRITICAL** — changes the reported quartile coef + FDR in the
  R2-2 male-only paragraph of the response letter, and the left panel of
  Supp. Fig. `supp_male_only_sensitivity`.
- **Location:** `male_only_sensitivity()` lines 273-315.
- **Mechanism:** `meta["smoke_binary"].notna()` filter at line 281 lets all
  6,676 rows through; the male subset of n=2,801 includes 402
  former_smoker + 121 not_everyday males (523 total) coded as
  `smoke_everyday=0`, contaminating the never-smoker baseline.
- **Fix:** done. New script
  `refit_male_quartile_module.py` re-fits the per-quartile model two ways:
  - Approach A (binary, drop former+notdaily) → n_male = 2,278
  - Approach B (4-level smk_status factor, all 2,801 males kept)
- **Outputs:**
  - `male_only_module_effects_by_age_corrected_everyday_vs_never.tsv`
  - `male_only_module_effects_by_age_corrected_4level.tsv`
  - re-rendered `supp_male_only_sensitivity.{png,pdf}`
- **Q3 status after fix:** stays non-significant under both approaches
  (Approach A: p=0.39, FDR=0.39; Approach B: p=0.24, FDR=0.24). The
  dampening of the smoking effect at age quartile 3 in males is a real
  pattern, not an artifact of the bug.

Numerical comparison:

| Quartile | Buggy coef | Buggy FDR | Approach A coef | Approach A FDR | Approach B coef | Approach B FDR |
|----------|-----------:|----------:|----------------:|---------------:|----------------:|---------------:|
| Q1       | 0.0995     | 2.28e-5   | 0.0950          | 1.09e-4        | 0.0985          | 3.94e-5        |
| Q2       | 0.0998     | 6.30e-5   | 0.1123          | 1.09e-4        | 0.1143          | 3.94e-5        |
| Q3       | 0.0192     | 0.4123    | 0.0239          | 0.3863         | 0.0313          | 0.2425         |
| Q4       | 0.0504     | 0.0210    | 0.0667          | 0.0128         | 0.0647          | 0.0101         |

The top-line story (Q1, Q2, Q4 significant; Q3 not significant) is
preserved across all three specifications. Magnitudes shift modestly.
Approach B is recommended for publication because it preserves the
4-level smk_status spec used in the published full-sample model.

### Finding 1.2 — OTU-level male-only OLS in `male_only_sensitivity()` (right panel of `supp_male_only_sensitivity.png`)

- **Severity: CRITICAL** — but already fixed by Bug-1 MaAsLin2 LOG re-run.
- **Note:** the original `run_revision_analysis.py` OTU loop at lines
  317-338 inherits BOTH bugs:
  - the same n=2,801 contaminated male subset (Bug 2), AND
  - the arcsin-sqrt vs MaAsLin2-LOG transform mismatch (Bug 1).
- The MaAsLin2 R script `run_maaslin2_male_only_otu.R` happens to also use
  the contaminated n=2,801 frame (it preserves the analytic frame of the
  buggy Python script and uses the 4-level smk_status factor, which is
  Approach B). Because the 4-level factor correctly partitions the 523
  contaminating males into former / not_everyday, the
  `smk_status everyday` coefficient is the right comparison to never. So
  the right panel of the figure is correct under Approach B logic.
- **Recommendation:** no further change needed for the right panel.

### Finding 1.3 — Family-level OLS smoking model in `run_taxonomic_robustness()`

- **Severity: HIGH** (changes a supplementary number; in
  `family_smoking_results.tsv`, n_samples = 6,496 is bug-contaminated).
  The MaAsLin2 R cross-check
  (`family_smoking_results_maaslin2.tsv`) is CLEAN (run from
  `run_maaslin2_reanalysis.R` which uses `ifelse` correctly), so the
  numerical conclusions in the response letter — derived from MaAsLin2 —
  are unaffected. But the Python TSV is wrong if reviewers read it.
- **Location:** `run_taxonomic_robustness()` lines 399-414.
- **Mechanism:** identical to 1.1 — `smoke_binary.notna()` at line 400
  lets the 591 former+not_everyday rows (with non-missing other
  covariates: 570 actually surviving) through with `smoke_everyday=0`.
- **Fix:** simplest fix is to change line 400 from
  `meta["smoke_binary"].notna()` to
  `meta["smk_status"].isin(["everyday", "never_smoker"])`.
- **Empirical n:** buggy n=6,496 → corrected n=5,926.
- **Recommendation for response letter:** quote the MaAsLin2 numbers
  (already clean) and either re-run Python OLS with the correct filter,
  or drop the Python OLS family/genus TSVs from the supplementary
  package since MaAsLin2 supersedes them.

### Finding 1.4 — Genus-level OLS smoking model in `run_taxonomic_robustness()`

- **Severity: HIGH** (same as 1.3, propagates to
  `genus_smoking_results.tsv`).
- **Mechanism / fix:** identical to 1.3.
- **MaAsLin2 cross-check:** clean (`genus_smoking_results_maaslin2.tsv`).

### Finding 1.5 — Age model in `run_taxonomic_robustness()` (within never_smoker)

- **Severity: NONE.**
- **Reason:** uses `meta["smk_status"].eq("never_smoker")` directly (line
  417), bypassing the buggy `smoke_binary`. Clean.

### Finding 1.6 — `clean_metadata()` itself (line 60) — root cause of Bug 2

- **Severity:** root cause.
- **Recommended permanent fix:** replace the np.where idiom with a
  boolean-mask approach using actual NaN. Example:
  ```python
  sb = pd.Series(pd.NA, index=meta.index, dtype="object")
  sb[meta["smk_status"].eq("everyday")] = "everyday"
  sb[meta["smk_status"].eq("never_smoker")] = "never_smoker"
  meta["smoke_binary"] = sb
  ```
  Or the equivalent `meta["smk_status"].where(...)` chain. The
  `refit_male_quartile_module.py:clean_metadata_corrected` shows the
  correct pattern.

### Finding 1.7 — Rarefaction / sequencing depth analysis

- **Severity: NONE.** `compute_rarefaction()` does not use `smoke_binary`
  or `smoke_everyday`; only `meta["ID"]`. Clean.

### Finding 1.8 — Module construction in `male_only_sensitivity()` (lines 287-300)

- **Severity: NONE for the score itself**, but CONTEXT for Finding 1.1.
- The arcsin-sqrt → orient-by-sign(age_coef) → z-score → mean-across-OTUs
  pipeline is fine, but the z-score normalisation is computed on the n=2,801
  male subset; under Approach A (n=2,278) the z-score uses a slightly
  different mean / SD. `refit_male_quartile_module.py` recomputes the
  z-score on Approach A's n=2,278 subset for internal consistency. Approach
  B uses the original n=2,801 z-score (which is what the published
  workflow does).
- **Recommendation:** document the choice in the methods paragraph.
  Approach B is what's plotted.

---

## File 2 — `run_maaslin2_reanalysis.R`

### Finding 2.1 — Family/genus MaAsLin2 family/genus reanalysis

- **Severity: NONE.**
- `meta$smoke_binary` is built with R's `ifelse(meta$smk_status == "everyday", "everyday", ifelse(meta$smk_status == "never_smoker", "never_smoker", NA_character_))`
  (lines 80-81). R's `NA_character_` is real NA, not the string "NA"; and
  `is.na(meta$smoke_binary)` correctly drops them. The smoking model
  (lines 161-165) filters by `!is.na(smoke_binary)`. Smoking subset
  reproduces the expected n=5,926 of correct everyday + never_smoker
  males-and-females.
- The age model uses `smk_status == "never_smoker"` directly. Clean.
- **No transform mismatch:** this script's coefficients live in the same
  MaAsLin2 LOG-LM space as the published main-paper analysis; cross-frame
  comparisons in `maaslin2_summary.tsv` and `*_maaslin2.tsv` files are
  apples-to-apples.

---

## File 3 — `run_maaslin2_male_only_otu.R`

### Finding 3.1 — Male-only OTU MaAsLin2 LOG re-run

- **Severity: NONE for transform; AMBIGUOUS for sample frame.**
- **Transform side:** clean. Coefficients are fitted in TSS+LOG+LM space
  to match the published full-sample model exactly. The Bug-1 fix is
  correct.
- **Sample frame:** the script intentionally preserves the same n=2,801
  male frame as the buggy Python pipeline (lines 102 and the explanatory
  comment), but uses the 4-level `smk_status` factor with reference =
  never_smoker. This is correct under Approach-B logic — the `everyday`
  coefficient is correctly contrasted against never_smoker, with
  former_smoker / not_everyday absorbed into separate factor levels.
- **Recommendation:** in the figure caption / methods, explicitly state
  that the male subset preserved by this script is the 2,801-male
  complete-case frame and that smk_status is treated as a 4-level factor
  with reference = never_smoker — i.e., flag that this is Approach B,
  not the binary-everyday-vs-never frame. This avoids future confusion
  about the n.

---

## File 4 — `run_mediation_panel_originalspec.py`

### Finding 4.1 — Multi-mediator panel under original gai_med.R spec

- **Severity: NONE.**
- **Sample frame:** mediation does NOT use `smoke_binary` or
  `smoke_everyday`. The exposure is `pro_aging_score` (continuous,
  z-scored). The module score is built over ALL 6,676 samples (no
  smoking filter), which is correct because the score is defined for
  everyone, not just smokers. The mediation regressions then run on the
  complete-case subset for each mediator + ASCVD.
- **Local `clean_metadata()` (lines 79-90) does NOT create
  `smoke_binary` or `smoke_everyday` — it leaves `smk_status` as a
  string column.** This script is independent of the buggy
  `clean_metadata` in `run_revision_analysis.py`. Clean.
- **Smoking is used inside the ASCVD computation (line 243):**
  `df["smk_status"].map({"never_smoker": 0, "current_smoker": 1})`. NB:
  this maps `"current_smoker"` to 1, but the metadata uses
  `"everyday"`, `"former_smoker"`, `"not_everyday"`, `"never_smoker"` —
  there is NO `"current_smoker"` value in this dataset. Therefore every
  non-`never_smoker` sample (including `everyday`) maps to NaN and is
  EXCLUDED from the ASCVD calculation entirely.
  - **Severity: HIGH (potentially CRITICAL — needs lead-author
    confirmation).** This means the ASCVD risk used as the mediator
    outcome is computed only on never_smoker rows (and the resulting
    NA-cells are dropped at the regression step). Approximately 4,679
    never_smokers minus those <40 or >79 years old or missing covariates
    will be the effective sample.
  - Likely intent was `df["smk_status"].map({"never_smoker": 0, "everyday": 1, "former_smoker": 1, "not_everyday": 1})`
    or just `(df["smk_status"] == "everyday").astype(int)` — the
    PCE smoker indicator is "current smoker" which most directly matches
    `"everyday"` (and possibly `"not_everyday"`).
  - **Recommendation: AMBIGUOUS — lead author needs to specify whether
    the published Figure 6b ASCVD calculation defined the smoker
    indicator the same way (i.e., only never_smokers contribute) or
    differently. If the panel currently runs on a much smaller cohort
    than expected, this would explain any unexpected n in
    `mediation_panel_results.tsv:n_complete`.
- **Transform mismatch:** none — all OLS regressions inside this script
  use the same Python OLS routine; coefficients are not compared to
  external R / MaAsLin2 outputs. Clean.

### Finding 4.2 — Pro-aging score construction in `build_proaging_score()`

- **Severity: NONE.**
- The score uses arcsin-sqrt → orient-by-age-sign → z-score → mean. This
  is identical to `male_only_sensitivity` and is the convention used in
  the response letter for the "module score" exposure. The score is
  computed over ALL 6,676 samples (no smoking filter), which is the
  intended population-level scaling.

---

## File 5 — `run_lightgbm_optimized.py`

### Finding 5.1 — Smoking classifier dataset

- **Severity: NONE.**
- `load_raw()` lines 110-122 filters by
  `smk.isin(["never_smoker", "everyday"])` directly on
  `meta["smk_status"]`. This bypasses the buggy `clean_metadata`
  entirely (this script does not import or use it). The 591 former /
  not_everyday rows are correctly dropped by `.isin(...)`. Clean.
- Labels are `(smk == "everyday").astype(int)`, so the binary
  classification is correctly defined as everyday-vs-never. Clean.
- **No transform mismatch:** the LightGBM model is end-to-end inside
  this script. Feature importances are LGBM-internal, not compared to
  MaAsLin2 / OLS coefficients. Clean.

---

## Summary table

| File                             | Function / region              | Severity     | Status                  |
|----------------------------------|--------------------------------|-------------:|-------------------------|
| run_revision_analysis.py         | clean_metadata (line 60)       | root cause   | NEEDS PERMANENT FIX     |
| run_revision_analysis.py         | male_only_sensitivity quartile | **CRITICAL** | FIXED via refit script  |
| run_revision_analysis.py         | male_only_sensitivity OTU      | CRITICAL     | FIXED via MaAsLin2 R    |
| run_revision_analysis.py         | family OLS smoking             | **HIGH**     | NEEDS RE-RUN OR DROP    |
| run_revision_analysis.py         | genus OLS smoking              | **HIGH**     | NEEDS RE-RUN OR DROP    |
| run_revision_analysis.py         | family/genus OLS age (never)   | NONE         | clean                   |
| run_revision_analysis.py         | rarefaction                    | NONE         | clean                   |
| run_maaslin2_reanalysis.R        | family/genus MaAsLin2 LOG      | NONE         | clean (use these numbers) |
| run_maaslin2_male_only_otu.R     | male-only MaAsLin2 LOG re-run  | NONE         | clean (Approach B)      |
| run_mediation_panel_originalspec.py | mediation panel             | NONE         | clean                   |
| run_mediation_panel_originalspec.py | ASCVD smoker mapping        | **HIGH/AMBIGUOUS** | NEEDS LEAD-AUTHOR CHECK |
| run_lightgbm_optimized.py        | smoking classifier             | NONE         | clean                   |

## Recommendations to lead author

1. **Update R2-2 male-only paragraph** with the corrected per-quartile
   numbers from
   `male_only_module_effects_by_age_corrected_4level.tsv`
   (Approach B, n=2,801 with 4-level factor): Q1 coef=0.098, FDR=3.94e-5;
   Q2 coef=0.114, FDR=3.94e-5; Q3 coef=0.031, FDR=0.24 (NS); Q4 coef=0.065,
   FDR=0.010. Or use Approach A (n=2,278, drop former+not_everyday) —
   numbers are quantitatively similar.
2. **Replace the supplementary figure** `supp_male_only_sensitivity.{png,pdf}`
   with the version produced by `refit_male_quartile_module.py`
   (already done; buggy backup preserved as
   `supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}`).
3. **Family / genus Python OLS results in
   `family_smoking_results.tsv` / `genus_smoking_results.tsv` are bug-
   contaminated** but the MaAsLin2 cross-check is clean. Either re-run
   the Python OLS with the correct filter (one-line fix in
   `run_taxonomic_robustness`) or drop the Python TSVs from the
   supplementary package, since MaAsLin2 is the published framework.
4. **Permanently fix `clean_metadata`** in `run_revision_analysis.py`
   line 60 to use a `pd.NA` / boolean-mask construction instead of
   `np.where(...np.nan)`. See
   `refit_male_quartile_module.py:clean_metadata_corrected` for the
   reference fix.
5. **Verify the ASCVD smoker mapping** in
   `run_mediation_panel_originalspec.py` line 243. The current code
   uses `{"never_smoker": 0, "current_smoker": 1}`, but the GGMP
   metadata does not contain `"current_smoker"`. Confirm whether the
   published Figure 6b used the same mapping (in which case nothing
   changes) or an `everyday`-mapped definition.
