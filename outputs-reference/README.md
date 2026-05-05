# Reference outputs

The TSV / JSON files in this directory are the **original-run results** from our own execution of the four analysis scripts on the GGMP dataset (May 2026). They are committed so that anyone re-running the pipeline can confirm a faithful reproduction by diffing their own `outputs/` files against these.

These are intentionally a small subset of the full output set — the large `*_maaslin2.tsv` files, the full pickle (`lightgbm_final_model.pkl`), the full Optuna trial history, and the analysis-metadata snapshot are **not** included to keep the repository lightweight.

## Files

| File | Source script | Description |
| --- | --- | --- |
| `sequencing_depth_summary.tsv` | `01_revision_analysis.py` | Per-subset sequencing-depth summary (analysis n, min/median/max reads). Underlies Supp Table 9. |
| `male_only_module_effects_by_age.tsv` | `01_revision_analysis.py` | Male-only adjusted effect of everyday smoking on the pro-aging module score, stratified by age quartile. Underlies Supp Table 11. |
| `family_smoking_results.tsv` | `01_revision_analysis.py` | Family-level OLS smoking-association results (one row per family). Underlies Supp Tables 12-13. |
| `genus_smoking_results.tsv` | `01_revision_analysis.py` | Genus-level OLS smoking-association results (one row per genus). Underlies Supp Tables 14-17. |
| `mediation_panel_results.tsv` | `03_mediation_panel.py` | Ten-mediator panel results (path-a/b/c coefficients, q-values, indirect effect with 95% bootstrap CI, proportion mediated). Underlies Supp Table 20. |
| `lightgbm_optimized_hyperparameters.json` | `04_lightgbm_optimized.py` | Optuna-winning LightGBM hyperparameters + meta (winning variant, CV AUC mean / SD, n_trials). |
| `lightgbm_optimized_cv_results.tsv` | `04_lightgbm_optimized.py` | Per-fold CV diagnostics (AUC + threshold-0.5 and Youden-thresholded confusion matrices, precision, recall, F1, MCC). Underlies Supp Table 21. |

## Diffing notes

Because LightGBM and several BLAS-backed numpy/scipy operations are platform-dependent, exact byte equality across machines is not expected. Acceptable matches:

- Sample counts (e.g. `n_complete`, `analysis_subset_n`) should match exactly.
- Coefficients / p-values / CV AUC should agree to ~3 decimal places.
- Top-ranked features (top-20 by gain in script 04, top-10 jointly significant mediators in script 03) should be identical or differ only in tail ranks.

If you see large divergences (e.g. CV AUC differing by >0.05, sample counts disagreeing), please open an issue with the OS / Python version / lightgbm version you are running.
