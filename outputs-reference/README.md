# Reference outputs

The TSV / JSON files in this directory are the **post-audit, manuscript-ready outputs** from our own execution of the analysis scripts on the GGMP dataset (May 2026). They are committed so that anyone re-running the pipeline can confirm a faithful reproduction by diffing their own `outputs/` files against these.

These are intentionally a small subset of the full output set — the large `*_maaslin2.tsv` files, the full pickle (`lightgbm_final_model.pkl`), the full Optuna trial history, and the analysis-metadata snapshot are **not** included to keep the repository lightweight.

> **Naming convention.** Where a fix script emits a `*_corrected*.tsv` or `*_corrected_4level.tsv` file alongside a buggy original, the version committed here is the **corrected** one, but it is committed under the **canonical, un-suffixed name** (e.g. `family_smoking_results.tsv` rather than `family_smoking_results_corrected.tsv`) — that is the name a re-runner who cherry-picks the manuscript-ready scripts (02, 02b, 04, 06, 07, 08) would naturally compare against. The "Source script" column in the table below names the script that actually produced the values.

## Files

| File | Source script | Description |
| --- | --- | --- |
| `sequencing_depth_summary.tsv` | `01_revision_analysis.py` | Per-subset sequencing-depth summary (analysis n, min/median/max reads). Underlies Supp Table 9. |
| `male_only_module_effects_by_age.tsv` | `06_refit_male_quartile_module.py` | Male-only adjusted effect of everyday smoking on the pro-aging module score, stratified by age quartile. Approach B (4-level smk_status, reference=never_smoker) — corrected for the smoke_binary string-NaN bug (n=2,801). Underlies Supp Table 10. Headline values: Q1 coef=0.098, Q2=0.114, Q3=0.031, Q4=0.065. |
| `male_only_overlap_otu_concordance_maaslin2.tsv` | `02b_maaslin2_male_only_otu.R` | Male-only MaAsLin2 LOG smoking coefficients for the 40 direction-concordant overlap OTUs (Bug 1 fix), aligned with the published full-sample LOG framework. Underlies Supp Table 11 and the right panel of Supp Fig 7. |
| `family_smoking_results.tsv` | `07_refit_family_genus_smoking_ols.py` | Family-level OLS smoking-association results (one row per family). Corrected for the smoke_binary string-NaN bug (n=5,926 vs the 6,496 in the buggy run). 7 family-level FDR-significant features. Underlies Supp Table 12. |
| `family_shared_results.tsv` | `07_refit_family_genus_smoking_ols.py` | Family-level shared smoking + age-within-never-smokers results (one row per family). Corrected smoking-side rejoined to the (clean) age-side TSV. Underlies Supp Table 14. |
| `genus_smoking_results.tsv` | `07_refit_family_genus_smoking_ols.py` | Genus-level OLS smoking-association results (one row per genus). Corrected (n=5,926). 15 genus-level FDR-significant features. Underlies Supp Table 15. |
| `genus_shared_results.tsv` | `07_refit_family_genus_smoking_ols.py` | Genus-level shared smoking + age-within-never-smokers results (one row per genus). Corrected. Underlies Supp Table 17. |
| `mediation_panel_results.tsv` | `08_refit_mediation_panel.py` | Ten-mediator panel results (path-a/b/c coefficients, q-values, indirect effect with 95% bootstrap CI, proportion mediated). Corrected for the ASCVD smoker-mapping bug — uses the published `cvrisk.R` mapping (n=4,642). SBP ranks #1 by absolute indirect effect with proportion mediated 39.7% (95% CI 32.4–48.7%); 8 of 10 mediators are jointly FDR-significant. Underlies Supp Table 20. |
| `lightgbm_optimized_hyperparameters.json` | `04_lightgbm_optimized.py` | Optuna-winning LightGBM hyperparameters + meta (winning variant, CV AUC mean / SD, n_trials). |
| `lightgbm_optimized_cv_results.tsv` | `04_lightgbm_optimized.py` | Per-fold CV diagnostics (AUC + threshold-0.5 and Youden-thresholded confusion matrices, precision, recall, F1, MCC). Underlies Supp Table 21. |

## Mapping: reference output → paper figure / table

| Reference output | Paper artefact |
| --- | --- |
| `sequencing_depth_summary.tsv` | Supp Table 9 |
| `male_only_module_effects_by_age.tsv` | Supp Table 10; left panel of Supp Fig 7 |
| `male_only_overlap_otu_concordance_maaslin2.tsv` | Supp Table 11; right panel of Supp Fig 7 |
| `family_smoking_results.tsv` | Supp Table 12; underlying scatter for Supp Fig 8 |
| `family_shared_results.tsv` | Supp Table 14 |
| `genus_smoking_results.tsv` | Supp Table 15; underlying scatter for Supp Fig 9 |
| `genus_shared_results.tsv` | Supp Table 17 |
| `mediation_panel_results.tsv` | Supp Table 20; Supp Fig 10 |
| `lightgbm_optimized_cv_results.tsv` | Supp Table 21; Supp Fig 11 |
| `lightgbm_optimized_hyperparameters.json` | Methods section + Supp Fig 11 caption |

## Diffing notes

Because LightGBM and several BLAS-backed numpy/scipy operations are platform-dependent, exact byte equality across machines is not expected. Acceptable matches:

- Sample counts (e.g. `n_complete`, `analysis_subset_n`) should match exactly.
- Coefficients / p-values / CV AUC should agree to ~3 decimal places.
- Top-ranked features (top-20 by gain in script 04, top-10 jointly significant mediators in script 08) should be identical or differ only in tail ranks.

If you see large divergences (e.g. CV AUC differing by >0.05, sample counts disagreeing), please open an issue with the OS / Python version / lightgbm version you are running.
