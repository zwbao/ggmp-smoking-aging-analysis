# Changelog

All notable changes to this repository are documented in this file. The repository tracks the analysis code for Bao et al. (NTR-2026-092); the corresponding manuscript reference and DOI will be added once the article is assigned a DOI.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) — though for a code-release-with-paper repo "version" really tracks the snapshot used for a particular submission round.

## [1.0.0] — 2026-05-05

First public release accompanying the peer-review revision (NTR-2026-092). Built up over three commits:

### Initial release (`78d0973`, 2026-05-05)
- Added the four original revision-analysis scripts: `01_revision_analysis.py` (rarefaction supplement, male-only sensitivity, family/genus OLS robustness), `02_maaslin2_reanalysis.R` (MaAsLin2 family/genus cross-check), `03_mediation_panel.py` (ten-mediator cardiometabolic panel), and `04_lightgbm_optimized.py` (Optuna-tuned LightGBM smoking classifier).
- Added repo scaffolding: `README.md`, `LICENSE` (MIT), `.gitignore`, `requirements.txt`, `requirements-r.txt`, `data/README.md`, `outputs-reference/` (small reference TSV/JSON outputs), `docs/reproducibility.md`.

### URL placeholder fix (`cd99763`, 2026-05-05)
- Replaced the `<this-repo>` placeholder in the README with the actual repository URL (`https://github.com/zwbao/ggmp-smoking-aging-analysis`).

### Four-bug audit fix (`7f71292`, 2026-05-05)
A cross-pipeline audit (`docs/audit/cross_pipeline_audit.md`) identified four sample-filter / encoding bugs in the original revision-analysis package. All four were fixed before the first public release; the affected supplementary tables, figures, and `outputs-reference/` files in this repository all reflect the post-fix state.

- **Bug 1 — male-only OTU transform mismatch.** Added `02b_maaslin2_male_only_otu.R` and `05_replot_male_only_figure.py`. Re-fits the male-only OTU smoking model in the published MaAsLin2 LOG framework (instead of Python OLS on `arcsin(sqrt(.))`) so the y=x diagonal of Supp Fig 7 right panel is meaningful. Memo: `docs/audit/male_only_otu_concordance_fix_memo.md`.
- **Bug 2 — male-only quartile `smoke_binary` string-NaN bug.** Added `06_refit_male_quartile_module.py`. Re-fits the male-only per-quartile pro-aging-module regression two ways (binary `everyday vs never`, 4-level `smk_status` factor); the manuscript uses Approach B (n=2,801, Q1=0.098, Q2=0.114, Q3=0.031, Q4=0.065). Re-renders Supp Fig 7. Memo: `docs/audit/cross_pipeline_audit.md` (Finding 1.1).
- **Bug 3 — family/genus OLS `smoke_binary` string-NaN bug (same root cause).** Added `07_refit_family_genus_smoking_ols.py`. Re-fits the family- and genus-level OLS smoking models with the corrected sample filter (n=5,926 instead of the buggy 6,496) and rejoins the (clean) age-side TSVs. Re-renders Supp Figs 8–9. Memo: `docs/audit/final_bug_fixes_memo.md` (Bug 3).
- **Bug 4 — ASCVD smoker mapping.** Added `08_refit_mediation_panel.py`. Adopts the published `cvrisk.R` smoker indicator (everyday + not_everyday → 1, never_smoker → 0, former_smoker excluded); the corrected mediation panel runs on n=4,642 complete cases (vs the buggy 3,426 never-smoker-only frame). SBP retains its #1 rank by absolute indirect effect; proportion mediated is 39.7% (95% CI 32.4–48.7%); 8 of 10 mediators are jointly FDR-significant. Memo: `docs/audit/final_bug_fixes_memo.md` (Bug 4).

### Pre-release polish (this commit, 2026-05-05)
- Added `CITATION.cff` and this `CHANGELOG.md`.
- Synced `outputs-reference/` with the post-audit corrected outputs and added three more reference TSVs: `male_only_overlap_otu_concordance_maaslin2.tsv` (Supp Fig 7B / Supp Table 11), `family_shared_results.tsv` (Supp Table 14), `genus_shared_results.tsv` (Supp Table 17).
- Added an explicit "reference output → paper figure / table" mapping table to `outputs-reference/README.md`.
- Cleaned up audit memos for a public-repo audience (removed leftover `revision-analysis/` working-tree paths and internal-process language; kept the technical content unchanged).
- Removed the unused R packages `dplyr`, `readr`, and `phyloseq` from `requirements-r.txt`; added `jsonlite` (used by `02b`).
- Clarified in `data/README.md` that `current_smoker` is a derived label, not a raw GGMP `smk_status` value.
- Various small README corrections (canonical-home URL near the top, post-publication DOI note, accurate package list).
