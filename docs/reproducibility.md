# Reproducibility guide

This document maps each script in [`scripts/`](../scripts) to the supplementary figures and tables in Bao et al. (NTR-2026-092), and lists the inputs, outputs, expected runtime, and random seeds used. All commands assume the repository root as the current working directory.

> **Bug fixes (May 2026).** A four-bug audit pass on the revision-analysis package produced four new fix scripts (`02b`, `06`, `07`, `08`) plus a figure-only re-render (`05`) that supersede sub-pipelines of `01` and `03` for the manuscript-ready outputs. The fix scripts emit `*_corrected*.tsv` files alongside the originals; the supplementary tables, figures, and `outputs-reference/` files in this repository all reflect the **post-fix** state. See [`audit/`](audit/) for full memos and the README "Audit history" section for a one-paragraph summary of each bug.

## Prerequisites

1. Place the data files under `./data/` as described in [`data/README.md`](../data/README.md).
2. Install Python and R dependencies as described in the project [`README.md`](../README.md).
3. Outputs default to `./outputs/` and run logs to `./logs/`. Both are created on first run.

## 1. `scripts/01_revision_analysis.py`

**Produces:**

- Supp Fig 6 + Supp Table 9 (rarefaction / sequencing depth)
- Supp Tables 13, 16 (family/genus age-side OLS — clean; never had a bug)
- `outputs/otu_overlap_from_original_results.tsv` consumed by all downstream scripts.

> The buggy male-only sensitivity (Supp Fig 7 + Supp Tables 10-11) and family/genus smoking-side OLS (Supp Fig 8-9 + Supp Tables 12, 14, 15, 17) sub-pipelines of this script have been superseded by `06_refit_male_quartile_module.py`, `02b_maaslin2_male_only_otu.R` + `05_replot_male_only_figure.py`, and `07_refit_family_genus_smoking_ols.py` (see Audit history in the project README).

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `data/smk_status_sig_res.tsv`, `data/age_sig_res.tsv`,
  `data/smk_amount_categ_sig_res.tsv`, `data/smk_y_categ_sig_res.tsv`

**Output files (under `outputs/`):**

- `sequencing_depth_summary.tsv`, `rarefaction_summary.tsv`, `supp_rarefaction_curve.{png,pdf}` — rarefaction supplement.
- `male_only_module_effects_by_age.tsv`, `male_only_overlap_otu_concordance.tsv`, `supp_male_only_sensitivity.{png,pdf}` — male-only sensitivity.
- `family_smoking_results.tsv`, `family_age_results.tsv`, `family_shared_results.tsv`, `supp_family_overlap_scatter.{png,pdf}` — family-level OLS robustness.
- `genus_smoking_results.tsv`, `genus_age_results.tsv`, `genus_shared_results.tsv`, `supp_genus_overlap_scatter.{png,pdf}` — genus-level OLS robustness.
- `otu_overlap_from_original_results.tsv` — OTU overlap consumed downstream by `03_mediation_panel.py`.
- `analysis_summary.md` — narrative summary.

**Run:**

```bash
python scripts/01_revision_analysis.py
```

**Seed:** rarefaction sub-sampling uses `seed=20260319` (overridable via `--seed`).

**Expected wall-clock runtime:** ~5 minutes on a modern laptop.

## 2. `scripts/02_maaslin2_reanalysis.R`

**Produces:**

- Supp Tables 18-19 (MaAsLin2 family / genus cross-check)

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`

**Output files (under `outputs/`):**

- `family_smoking_results_maaslin2.tsv`, `family_age_results_maaslin2.tsv`, `family_shared_results_maaslin2.tsv`
- `genus_smoking_results_maaslin2.tsv`, `genus_age_results_maaslin2.tsv`, `genus_shared_results_maaslin2.tsv`
- `maaslin2_summary.tsv` — overall summary across levels.
- `outputs/maaslin2_runs/` — MaAsLin2 native run directories.

**Run:**

```bash
Rscript scripts/02_maaslin2_reanalysis.R
```

**Seed:** none required; MaAsLin2 LM fits are deterministic given the inputs.

**Expected wall-clock runtime:** ~10 minutes.

## 3. `scripts/03_mediation_panel.py` (superseded — kept for reproducibility of the original buggy run)

> **This script's outputs are superseded by `08_refit_mediation_panel.py` (Bug 4 fix).** It is retained only so that anyone can reproduce the as-submitted-then-corrected buggy state of `mediation_panel_results.tsv`. Do not use the outputs of this script for any downstream comparison; use `outputs/mediation_panel_results_corrected.tsv` from script 08 instead.

**Produces (buggy):**

- `outputs/mediation_panel_results.tsv`
- `outputs/mediation_panel_summary.md`
- `outputs/supp_mediation_panel.{png,pdf}` (overwritten by script 08 in a normal run)

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/otu_overlap_from_original_results.tsv` (produced by `01_revision_analysis.py`; **run script 01 first**).

**Run:**

```bash
python scripts/03_mediation_panel.py
```

**Seed:** `seed=12345`; 1,000 bootstrap resamples.

**Expected wall-clock runtime:** ~10 minutes.

## 4. `scripts/04_lightgbm_optimized.py`

**Produces:**

- Supp Fig 11 + Supp Table 21 (Optuna LightGBM hyperparameter optimization)

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `data/smk_status_sig_res.tsv` (used to pick the top-71 / top-150 / q<0.05 OTU subsets).

**Output files (under `outputs/`):**

- `lightgbm_optimized_hyperparameters.json` — winning Optuna hyperparameters + meta.
- `lightgbm_optimized_cv_results.tsv` — per-fold AUC + confusion-matrix diagnostics.
- `lightgbm_optimized_feature_importance.tsv` — top-20 gain importances.
- `lightgbm_optimization_history.tsv` — full Optuna trial history.
- `lightgbm_optimized_variant_screen.tsv` — baseline AUC across the 9 input variants.
- `lightgbm_final_model.pkl` — pickled refit model + feature columns + best params.
- `lightgbm_optimized.png` — Optuna history + per-fold CV AUC plots.
- `lightgbm_optimized_summary.md` — narrative summary.
- `logs/lightgbm_optimized_run.log` — full run log.

**Run:**

```bash
python scripts/04_lightgbm_optimized.py
```

**Seeds:** `TPE_SEED = 20260504` (Optuna), `CV_SEED = 123` (`StratifiedKFold(random_state=123)`), LGBM `random_state=20260504`.

**Expected wall-clock runtime:** ~30 minutes (200-trial / 30-minute Optuna cap).

## 5. `scripts/02b_maaslin2_male_only_otu.R` (Bug 1 fix)

**Produces:** Supp Table 11 (male-only OTU concordance, MaAsLin2 LOG framework) and the right-panel data of Supp Fig 7.

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/male_only_overlap_otu_concordance.tsv` (produced by script 01; carries the buggy arcsin-sqrt male coefs, but we only use its `feature` / `smoking_coef` / `age_coef` / `direction_concordant` columns — these are correct).

**Output files (under `outputs/`):**

- `male_only_overlap_otu_concordance_maaslin2.tsv` — corrected male-only MaAsLin2 LOG smoking coefficients for the 40 direction-concordant overlap OTUs.
- `male_only_concordance_maaslin2_summary.json` — summary stats (concordance counts, Spearman/Pearson old-vs-new male coef).
- `outputs/maaslin2_runs/male_only_otu_smoking/` — MaAsLin2 native run directory.

**Run:**

```bash
Rscript scripts/02b_maaslin2_male_only_otu.R
```

**Expected wall-clock runtime:** ~5 minutes.

## 6. `scripts/05_replot_male_only_figure.py` (Bug 1 figure re-render)

**Produces:** the post-Bug-1 version of `supp_male_only_sensitivity.{png,pdf}` (right panel only — left panel is overwritten again by script 06).

**Inputs (defaults):**

- `outputs/male_only_module_effects_by_age.tsv` (consumed for the left panel; later replaced by script 06).
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (from script 02b).

**Output files (under `outputs/`):**

- `supp_male_only_sensitivity.{png,pdf}` — re-rendered.
- `supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}` — backup of the buggy original.

**Run:**

```bash
python scripts/05_replot_male_only_figure.py
```

**Expected wall-clock runtime:** seconds.

## 7. `scripts/06_refit_male_quartile_module.py` (Bug 2 fix)

**Produces:** Supp Table 10 (male-only quartile sensitivity, both Approach A and Approach B) and re-renders Supp Fig 7 (left panel = Approach B forest, right panel = MaAsLin2 LOG concordance scatter).

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/otu_overlap_from_original_results.tsv` (from script 01).
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (from script 02b).

**Output files (under `outputs/`):**

- `male_only_module_effects_by_age_corrected_4level.tsv` — Approach B (used for Supp Table 10).
- `male_only_module_effects_by_age_corrected_everyday_vs_never.tsv` — Approach A (sensitivity).
- `supp_male_only_sensitivity.{png,pdf}` — re-rendered (Approach B left + MaAsLin2 right).

**Run:**

```bash
python scripts/06_refit_male_quartile_module.py
```

**Expected wall-clock runtime:** ~1 minute.

## 8. `scripts/07_refit_family_genus_smoking_ols.py` (Bug 3 fix)

**Produces:** Supp Tables 12, 14, 15, 17 (family/genus smoking-side and shared OLS — corrected n=5,926 instead of buggy 6,496) and re-renders Supp Fig 8-9.

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/family_age_results.tsv`, `outputs/genus_age_results.tsv` (from script 01 — clean; rejoined to the corrected smoking results to rebuild the shared TSVs).

**Output files (under `outputs/`):**

- `family_smoking_results_corrected.tsv`, `family_shared_results_corrected.tsv`
- `genus_smoking_results_corrected.tsv`, `genus_shared_results_corrected.tsv`
- `supp_family_overlap_scatter.{png,pdf}`, `supp_genus_overlap_scatter.{png,pdf}` — re-rendered.
- Buggy backups preserved as `*_buggy_smoke_binary.{tsv,png,pdf}`.

**Run:**

```bash
python scripts/07_refit_family_genus_smoking_ols.py
```

**Expected wall-clock runtime:** ~2 minutes.

## 9. `scripts/08_refit_mediation_panel.py` (Bug 4 fix)

**Produces:** Supp Fig 10 + Supp Table 20 (multi-mediator panel) with the corrected ASCVD smoker mapping per `cvrisk.R` (n=4,642 complete cases vs the buggy n=3,426 never-smoker-only cases).

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/otu_overlap_from_original_results.tsv` (from script 01).

**Output files (under `outputs/`):**

- `mediation_panel_results_corrected.tsv` — full ten-mediator results table (used for Supp Table 20).
- `mediation_panel_summary_corrected.md` — narrative summary with publication-ready Methods paragraph.
- `supp_mediation_panel.{png,pdf}` — forest plot of indirect effects (re-rendered).
- Buggy backups preserved as `supp_mediation_panel_smoker_buggy.{png,pdf}`.

**Run:**

```bash
python scripts/08_refit_mediation_panel.py
```

**Seed:** `seed=12345` (matches the published `gai_med.R` `set.seed(12345)`); 1,000 bootstrap resamples.

**Expected wall-clock runtime:** ~10 minutes.

## Suggested order of execution

```
1.  python  scripts/01_revision_analysis.py             # ~5 min   (rarefaction, age-side OLS, OTU overlap)
2.  Rscript scripts/02_maaslin2_reanalysis.R            # ~10 min  (MaAsLin2 family/genus cross-check)
3.  Rscript scripts/02b_maaslin2_male_only_otu.R        # ~5 min   (Bug 1 fix)
4.  python  scripts/05_replot_male_only_figure.py       # ~10 sec  (Bug 1 figure)
5.  python  scripts/06_refit_male_quartile_module.py    # ~1 min   (Bug 2 fix; final Supp Fig 7)
6.  python  scripts/07_refit_family_genus_smoking_ols.py # ~2 min  (Bug 3 fix)
7.  python  scripts/08_refit_mediation_panel.py         # ~10 min  (Bug 4 fix; manuscript-ready Supp Table 20)
8.  python  scripts/04_lightgbm_optimized.py            # ~30 min  (independent of bug fixes)
```

Script `03_mediation_panel.py` is intentionally omitted from the suggested order because it has been superseded by script 08. Total wall-clock: ~1.2 hours on a 2026-era laptop.

## Cross-checking results

A small set of reference output TSVs from our own run is provided under [`outputs-reference/`](../outputs-reference). Compare the corresponding files in your own `outputs/` against these to confirm a faithful re-run. Tiny floating-point differences are expected because of platform-dependent BLAS / LightGBM threading; the broad numbers (sample counts, fold AUCs to 2-3 decimals, top-feature ranks) should match.

## Random seed cheatsheet

| Component | Variable | Value |
| --- | --- | --- |
| Rarefaction sub-sampling (script 01) | `--seed` | `20260319` |
| Mediation bootstrap (script 03 / script 08) | `--seed` | `12345` |
| Optuna TPE sampler (script 04) | `TPE_SEED` | `20260504` |
| StratifiedKFold (script 04) | `CV_SEED` | `123` |
| LightGBM `random_state` (script 04) | `TPE_SEED` | `20260504` |
