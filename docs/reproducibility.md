# Reproducibility guide

This document maps each script in [`scripts/`](../scripts) to the figures and tables in Bao et al. (NTR-2026-092), and lists the inputs, outputs, expected runtime, and random seeds used. All commands assume the repository root as the current working directory.

The repository is split into:

- `scripts/main/` — main-paper analyses for Figures 1-6 and Supp Figs 1-5.
- `scripts/revision/` — peer-review-period sensitivity / robustness analyses for Supp Figs 6-11 and Supp Tables 9-21.

> **Bug fixes (May 2026).** A four-bug audit pass on the revision-analysis package produced four new fix scripts (`02b`, `06`, `07`, `08`) plus a figure-only re-render (`05`) that supersede sub-pipelines of `01` and `03` for the manuscript-ready outputs. The fix scripts emit `*_corrected*.tsv` files alongside the originals; the supplementary tables, figures, and `outputs-reference/` files in this repository all reflect the **post-fix** state. See [`audit/`](audit/) for full memos and the README "Audit history" section for a one-paragraph summary of each bug.

## Prerequisites

1. Place the data files under `./data/` as described in [`data/README.md`](../data/README.md).
2. Install Python and R dependencies as described in the project [`README.md`](../README.md).
3. Outputs default to `./outputs/` and run logs to `./logs/`. Both are created on first run.

# Reproducing the main-paper figures

These scripts produce the OTU-level signature lists, the LightGBM smoking classifier, the GSEA enrichment, the ASCVD risk score, and the mediation analysis underlying Figures 1-6 and Supp Figs 1-5 of the published manuscript.

> **Important: phyloseq-object prerequisite.** All `scripts/main/*.R` scripts expect an analysis-ready phyloseq object (`non_mach_abx_GPf.rel`) already in the R session. The construction of this object depends on auxiliary `.rda` files (~1.5-2 GB each) maintained in the corresponding author's working tree (`220914/221008.rda`, `221019.rda`, etc.) and is **NOT redistributed** in this repository. To reproduce, follow the steps documented at the top of `scripts/main/01_main_maaslin2.R`: obtain `GGMP7009_even10k.biom` from the GGMP processing repo (<https://github.com/SMUJYYXB/GGMP-Regional-variations>), obtain `GPf_metadata.tsv` from the corresponding author, build `GPf.rel`, and apply the antibiotic-use exclusion + derived-variable construction in lines 19-110 of `01_main_maaslin2.R`. See also `data/README.md`.

## M1. `scripts/main/01_main_maaslin2.R`

**Produces:** Figures 1, 2, 3 OTU lists (222 first-hand smoking + 117 second-hand + 330 age-within-never-smokers); Supp Figs 1, 2, 4, 5 cross-stratum panels.

**Inputs:** analysis-ready phyloseq object `non_mach_abx_GPf.rel` (see above).

**Outputs:** nine MaAsLin2 native run directories (one per `Maaslin2()` call) under the working directory; in-memory `ggplot2` objects for the violin / Venn panels.

**Run:**

```bash
Rscript scripts/main/01_main_maaslin2.R
```

**Seed:** none required (MaAsLin2 LM fits are deterministic).

**Expected wall-clock runtime:** ~30 minutes (9 MaAsLin2 fits with `cores=14`; lower `cores` on smaller machines).

## M2. `scripts/main/02_lightgbm_smk.ipynb`

**Produces:** Supp Fig 3 (LightGBM smoking classifier ROC, AUC ≈ 0.73; SHAP feature-importance panel).

**Inputs:** `data/predict_smk.txt` — tab-separated smoker-classifier feature matrix (7,009 x 72; 71 OTU relative-abundance columns + the `Districts` target column encoding binary smoker / never-smoker; see notebook header for derivation).

**Outputs:** ROC and SHAP panels rendered inline.

**Run:** open the notebook in JupyterLab / Jupyter Notebook with a PyCaret 2.x kernel and execute cells top-to-bottom.

**Seeds:** `dataset.sample(random_state=786)`, `setup(session_id=123)`, `tune_model(n_iter=20)`.

**Expected wall-clock runtime:** ~5-10 minutes interactive.

## M3. `scripts/main/03_gsea_disease.R`

**Produces:** Figure 5 (GSEA-style enrichment of smoking- and age-related OTU lists across atherosclerosis, fatty liver, T2DM, hepatic calculus, gout, MetS).

**Inputs:** `otu.mets.meta`, `otu.mets.count`, `dis_gmt2` already in the R session (constructed in the corresponding author's working tree from the analysis-ready phyloseq object plus the OTU lists from script M1; see script header for construction).

**Outputs:** six GSEA result objects (`mets.y3.*`) and six `gseaplot2` ggplot objects (`mets.gsea.p3.*`); `.pptx` exports under `./pic/`.

**Run:**

```bash
Rscript scripts/main/03_gsea_disease.R
```

**Seed:** `clusterProfiler::GSEA` uses an internal RNG; set a seed before each `GSEA()` call for fully reproducible p-values. `nPermSimple = 10000`, `eps = 1e-50`, `pvalueCutoff = 2`.

**Expected wall-clock runtime:** ~10 minutes.

## M4. `scripts/main/04_cvrisk_ascvd.R`

**Produces:** Figure 6a (|GAI| vs ASCVD 10-year risk score scatter); the per-subject `ascvd` column consumed by `05_mediation_main.R`.

**Inputs:** `non_mach_abx_GPf.rel.meta` (constructed by script M1) plus the `abs_age_gap_adjust` column from the LightGBM age-prediction model maintained in the corresponding author's tree.

**Outputs:** `cvrisk_df`, `med_ana_dat2` data.frames; `ascvd_gai_p1` ggplot object.

**Run:**

```bash
Rscript scripts/main/04_cvrisk_ascvd.R
```

**Seed:** none (`CVrisk::ascvd_10y_accaha` is deterministic).

**Expected wall-clock runtime:** ~1 minute.

## M5. `scripts/main/05_mediation_main.R`

**Produces:** Figure 6b (mediation diagram); `res` data.frame of ACME / ADE / proportion-mediated for the 6-mediator panel.

**Inputs:** `med_ana_dat2` (from script M4) plus `gai_plot_df`, `lifestyle_name`, `Biochemistry_name`, `gai_with_other.sig` already in the R session (see script header).

**Outputs:** `lifestyle_biochem_df`, `path` (lavaan SEM fit), `res` (panel mediation results table).

**Run:**

```bash
Rscript scripts/main/05_mediation_main.R
```

**Seeds:** `set.seed(12345)` at the top of the script and again immediately before the panel mediation step (matches the published `gai_med.R` setting). Bootstrap resamples: `sims = 1000`, `boot = TRUE`. lavaan SEM `bootstrap = 1000`.

**Expected wall-clock runtime:** ~10 minutes.

# Reproducing the revision-paper supplements

The revision-period scripts under `scripts/revision/` produce the supplementary figures and tables introduced during peer review. They are self-contained and only depend on the data files placed under `data/` (no R-session prerequisites).

## R1. `scripts/revision/01_revision_analysis.py`



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
python scripts/revision/01_revision_analysis.py
```

**Seed:** rarefaction sub-sampling uses `seed=20260319` (overridable via `--seed`).

**Expected wall-clock runtime:** ~5 minutes on a modern laptop.

## R2. `scripts/revision/02_maaslin2_reanalysis.R`

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
Rscript scripts/revision/02_maaslin2_reanalysis.R
```

**Seed:** none required; MaAsLin2 LM fits are deterministic given the inputs.

**Expected wall-clock runtime:** ~10 minutes.

## R3. `scripts/revision/03_mediation_panel.py` (superseded — kept for reproducibility of the original buggy run)

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
python scripts/revision/03_mediation_panel.py
```

**Seed:** `seed=12345`; 1,000 bootstrap resamples.

**Expected wall-clock runtime:** ~10 minutes.

## R4. `scripts/revision/04_lightgbm_optimized.py`

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
python scripts/revision/04_lightgbm_optimized.py
```

**Seeds:** `TPE_SEED = 20260504` (Optuna), `CV_SEED = 123` (`StratifiedKFold(random_state=123)`), LGBM `random_state=20260504`.

**Expected wall-clock runtime:** ~30 minutes (200-trial / 30-minute Optuna cap).

## R5. `scripts/revision/02b_maaslin2_male_only_otu.R` (Bug 1 fix)

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
Rscript scripts/revision/02b_maaslin2_male_only_otu.R
```

**Expected wall-clock runtime:** ~5 minutes.

## R6. `scripts/revision/05_replot_male_only_figure.py` (Bug 1 figure re-render)

**Produces:** the post-Bug-1 version of `supp_male_only_sensitivity.{png,pdf}` (right panel only — left panel is overwritten again by script 06).

**Inputs (defaults):**

- `outputs/male_only_module_effects_by_age.tsv` (consumed for the left panel; later replaced by script 06).
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (from script 02b).

**Output files (under `outputs/`):**

- `supp_male_only_sensitivity.{png,pdf}` — re-rendered.
- `supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}` — backup of the buggy original.

**Run:**

```bash
python scripts/revision/05_replot_male_only_figure.py
```

**Expected wall-clock runtime:** seconds.

## R7. `scripts/revision/06_refit_male_quartile_module.py` (Bug 2 fix)

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
python scripts/revision/06_refit_male_quartile_module.py
```

**Expected wall-clock runtime:** ~1 minute.

## R8. `scripts/revision/07_refit_family_genus_smoking_ols.py` (Bug 3 fix)

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
python scripts/revision/07_refit_family_genus_smoking_ols.py
```

**Expected wall-clock runtime:** ~2 minutes.

## R9. `scripts/revision/08_refit_mediation_panel.py` (Bug 4 fix)

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
python scripts/revision/08_refit_mediation_panel.py
```

**Seed:** `seed=12345` (matches the published `gai_med.R` `set.seed(12345)`); 1,000 bootstrap resamples.

**Expected wall-clock runtime:** ~10 minutes.

## Suggested order of execution

```
# main-paper (R session prerequisite: load construction of `non_mach_abx_GPf.rel`):
M1. Rscript scripts/main/01_main_maaslin2.R              # ~30 min  (MaAsLin2 OTU lists; Fig 1-3, Supp Figs 1, 2, 4, 5)
M2. jupyter notebook scripts/main/02_lightgbm_smk.ipynb  # interactive (Supp Fig 3)
M3. Rscript scripts/main/03_gsea_disease.R               # ~10 min  (Fig 5)
M4. Rscript scripts/main/04_cvrisk_ascvd.R               # ~1 min   (Fig 6a)
M5. Rscript scripts/main/05_mediation_main.R             # ~10 min  (Fig 6b)

# revision-period:
R1. python  scripts/revision/01_revision_analysis.py             # ~5 min   (rarefaction, age-side OLS, OTU overlap)
R2. Rscript scripts/revision/02_maaslin2_reanalysis.R            # ~10 min  (MaAsLin2 family/genus cross-check)
R5. Rscript scripts/revision/02b_maaslin2_male_only_otu.R        # ~5 min   (Bug 1 fix)
R6. python  scripts/revision/05_replot_male_only_figure.py       # ~10 sec  (Bug 1 figure)
R7. python  scripts/revision/06_refit_male_quartile_module.py    # ~1 min   (Bug 2 fix; final Supp Fig 7)
R8. python  scripts/revision/07_refit_family_genus_smoking_ols.py # ~2 min  (Bug 3 fix)
R9. python  scripts/revision/08_refit_mediation_panel.py         # ~10 min  (Bug 4 fix; manuscript-ready Supp Table 20)
R4. python  scripts/revision/04_lightgbm_optimized.py            # ~30 min  (independent of bug fixes)
```

Script `03_mediation_panel.py` is intentionally omitted from the suggested order because it has been superseded by script 08. Total wall-clock for the revision-period scripts: ~1.2 hours on a 2026-era laptop. Add ~1 hour more for the full main-paper run (assuming the analysis-ready phyloseq object is already constructed).

## Cross-checking results

A small set of reference output TSVs from our own run is provided under [`outputs-reference/`](../outputs-reference). Compare the corresponding files in your own `outputs/` against these to confirm a faithful re-run. Tiny floating-point differences are expected because of platform-dependent BLAS / LightGBM threading; the broad numbers (sample counts, fold AUCs to 2-3 decimals, top-feature ranks) should match.

## Random seed cheatsheet

| Component | Variable | Value |
| --- | --- | --- |
| Main-paper LightGBM (script M2) train/test split | `random_state` (`pd.DataFrame.sample`) | `786` |
| Main-paper LightGBM (script M2) PyCaret session | `session_id` | `123` |
| Main-paper mediation (script M5) | `set.seed(...)` | `12345` |
| Rarefaction sub-sampling (script R1) | `--seed` | `20260319` |
| Mediation bootstrap (scripts R3 / R9) | `--seed` | `12345` |
| Optuna TPE sampler (script R4) | `TPE_SEED` | `20260504` |
| StratifiedKFold (script R4) | `CV_SEED` | `123` |
| LightGBM `random_state` (script R4) | `TPE_SEED` | `20260504` |
