# Reproducibility guide

This document maps each script in [`scripts/`](../scripts) to the supplementary figures and tables in Bao et al. (NTR-2026-092), and lists the inputs, outputs, expected runtime, and random seeds used. All commands assume the repository root as the current working directory.

## Prerequisites

1. Place the data files under `./data/` as described in [`data/README.md`](../data/README.md).
2. Install Python and R dependencies as described in the project [`README.md`](../README.md).
3. Outputs default to `./outputs/` and run logs to `./logs/`. Both are created on first run.

## 1. `scripts/01_revision_analysis.py`

**Produces:**

- Supp Fig 6 + Supp Table 9 (rarefaction / sequencing depth)
- Supp Fig 7 + Supp Tables 10-11 (male-only sensitivity)
- Supp Figs 8-9 + Supp Tables 12-17 (family/genus OLS robustness)

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

## 3. `scripts/03_mediation_panel.py`

**Produces:**

- Supp Fig 10 + Supp Table 20 (multi-mediator panel)

**Inputs (defaults):**

- `data/GGMP7009_even10k.biom`
- `data/GPf_metadata.tsv`
- `outputs/otu_overlap_from_original_results.tsv` (produced by `01_revision_analysis.py`; **run script 01 first**).

**Output files (under `outputs/`):**

- `mediation_panel_results.tsv` — full ten-mediator results table.
- `mediation_panel_summary.md` — narrative summary with publication-ready Methods paragraph.
- `supp_mediation_panel.{png,pdf}` — forest plot of indirect effects.

**Run:**

```bash
python scripts/03_mediation_panel.py
```

**Seed:** `seed=12345` (matches the published `gai_med.R` `set.seed(12345)`); 1,000 bootstrap resamples.

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

## Suggested order of execution

```
1. python scripts/01_revision_analysis.py          # ~5 min
2. Rscript scripts/02_maaslin2_reanalysis.R        # ~10 min
3. python scripts/03_mediation_panel.py            # ~10 min  (depends on script 01 output)
4. python scripts/04_lightgbm_optimized.py         # ~30 min
```

Total wall-clock: ~1 hour on a 2026-era laptop.

## Cross-checking results

A small set of reference output TSVs from our own run is provided under [`outputs-reference/`](../outputs-reference). Compare the corresponding files in your own `outputs/` against these to confirm a faithful re-run. Tiny floating-point differences are expected because of platform-dependent BLAS / LightGBM threading; the broad numbers (sample counts, fold AUCs to 2-3 decimals, top-feature ranks) should match.

## Random seed cheatsheet

| Component | Variable | Value |
| --- | --- | --- |
| Rarefaction sub-sampling (script 01) | `--seed` | `20260319` |
| Mediation bootstrap (script 03) | `--seed` | `12345` |
| Optuna TPE sampler (script 04) | `TPE_SEED` | `20260504` |
| StratifiedKFold (script 04) | `CV_SEED` | `123` |
| LightGBM `random_state` (script 04) | `TPE_SEED` | `20260504` |
