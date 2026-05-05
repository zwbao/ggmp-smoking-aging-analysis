# Smoking, aging, and the gut microbiome — analysis code for Bao et al., NTR-2026-092

**Canonical home:** <https://github.com/zwbao/ggmp-smoking-aging-analysis>

This repository contains the analysis code for the supplementary results of:

> Bao Z, Yang Z, Sun R, Meng R, Wu W, Li MD. **Associations of smoking, aging, and their interplay with the gut microbiome and chronic disease risk profiles.** *Nicotine & Tobacco Research* (2026). DOI: TBD (will be updated once the article DOI is assigned).

It is provided so that anyone with access to the public Guangdong Gut Microbiome Project (GGMP) data can re-run the analyses underlying the published manuscript. The original GGMP raw-processing workflow is *not* included here — see [Data dependencies](#data-dependencies).

## Repository scope

This repository covers:

- **Main-paper analyses** for Figures 1-6 and Supplementary Figures 1-5 (`scripts/main/`): MaAsLin2 OTU-level smoking and age models, the LightGBM smoking classifier (Supp Fig 3), GSEA-style chronic-disease enrichment (Fig 5), the ACC/AHA 2013 PCE ASCVD computation (Fig 6a), and the cardiometabolic mediation panel (Fig 6b).
- **Revision-period sensitivity / robustness analyses** for Supplementary Figures 6-11 and Supplementary Tables 9-21 (`scripts/revision/`): rarefaction supplement, male-only sensitivity, family- and genus-level OLS robustness, MaAsLin2 family/genus cross-check, ten-mediator cardiometabolic panel, and Optuna LightGBM hyperparameter optimization, plus the four post-audit fix scripts.

Raw FASTQ → BIOM-table processing is in the original GGMP repository (<https://github.com/SMUJYYXB/GGMP-Regional-variations>) and is not duplicated here.

## Audit history

In May 2026 a four-bug cross-pipeline audit was performed on the revision-analysis package. All four bugs have been fixed; the supplementary tables, figures, and the `outputs-reference/` files in this repository reflect the **post-fix** state. Buggy outputs are preserved in the local working tree for diagnostic record but are not redistributed here. Full memos: [`docs/audit/`](docs/audit/).

1. **Male-only OTU concordance — transform-space mismatch.** The right panel of Supp Fig 7 plotted male-only OTU smoking coefficients computed by Python OLS on `arcsin(sqrt(rel_abundance))` against full-sample MaAsLin2 `LOG`-transform coefficients from `smk_status_sig_res.tsv`. The two coefficient sets lived in different statistical spaces, so the y=x diagonal was meaningless and male-only points appeared spuriously attenuated by ~17x. Fix: re-run the male-only OTU smoking model in the published MaAsLin2 LOG framework (`scripts/revision/02b_maaslin2_male_only_otu.R`) and re-render the figure (`scripts/revision/05_replot_male_only_figure.py`). See [`docs/audit/male_only_otu_concordance_fix_memo.md`](docs/audit/male_only_otu_concordance_fix_memo.md).

2. **Male-only quartile encoding — `smoke_binary` string-NaN bug.** `clean_metadata` builds `smoke_binary` via a nested `np.where(...)` whose inner third branch coerces to the literal Python string `'nan'` rather than `NaN`. Downstream `.notna()` filters therefore let 523 former-smoker / not-everyday male rows leak into the male-only quartile model with `smoke_everyday=0`, contaminating the never-smoker baseline. Fix: `scripts/revision/06_refit_male_quartile_module.py` re-fits per-quartile two ways (binary `everyday vs never`, 4-level `smk_status` factor); the manuscript uses Approach B (4-level factor, n=2,801). See [`docs/audit/cross_pipeline_audit.md`](docs/audit/cross_pipeline_audit.md) (Finding 1.1).

3. **Family/genus OLS encoding — same `smoke_binary` string-NaN bug.** The same string-NaN bug propagated into `run_taxonomic_robustness`, contaminating the family- and genus-level smoking-side OLS models with the same 570 rows of former/non-daily smokers (n inflated from the intended 5,926 → buggy 6,496). Fix: `scripts/revision/07_refit_family_genus_smoking_ols.py` re-fits both levels with the corrected `smk_status.isin(["everyday", "never_smoker"])` filter and rejoins the (untouched) age-side TSVs. The age-side TSVs (`family_age_results.tsv`, `genus_age_results.tsv`) used a different code path and were unaffected. See [`docs/audit/final_bug_fixes_memo.md`](docs/audit/final_bug_fixes_memo.md) (Bug 3).

4. **ASCVD smoker mapping — non-existent category.** `run_mediation_panel.py` mapped the GGMP `smk_status` column with `{"never_smoker": 0, "current_smoker": 1}`, but `current_smoker` does not exist in GGMP — it is a *derived* category in the original `code.R`/`cvrisk.R` (`everyday` ∪ `not_everyday` → `current_smoker`). Under the buggy mapping, all `everyday`, `not_everyday`, and `former_smoker` rows were silently sent to NaN, so ASCVD existed only for never-smokers (n=3,426). Fix: `scripts/revision/08_refit_mediation_panel.py` adopts the published `cvrisk.R` definition exactly (`never_smoker → 0`, `everyday → 1`, `not_everyday → 1`, `former_smoker → excluded`); the corrected mediation panel runs on n=4,642 complete cases. See [`docs/audit/final_bug_fixes_memo.md`](docs/audit/final_bug_fixes_memo.md) (Bug 4) and [`docs/audit/mediation_panel_summary_corrected.md`](docs/audit/mediation_panel_summary_corrected.md).

## Abstract

**Background:** Smoking and aging are both linked to chronic disease and to variation in the gut microbiome, but their joint relationship with microbiome composition remains incompletely characterized. This study examined the overlap between smoking- and age-associated gut microbial signals and their associations with cardiovascular risk profiles. **Methods:** This study was a secondary analysis of the published Guangdong Gut Microbiome Project (GGMP) resource. We analyzed 6,676 participants from a filtered project-specific metadata subset together with the public GGMP OTU table. Associations of smoking phenotypes and age with gut microbiota were examined using multivariable linear models. Machine-learning and mediation analyses were treated as exploratory. Additional sensitivity analyses were performed in males only, and robustness analyses were repeated at the family and genus levels. **Results:** We identified 222 OTUs associated with first-hand smoking and 117 OTUs associated with second-hand smoke exposure. Among never smokers, age was associated with 330 OTUs. Eighty-five OTUs overlapped between smoking- and age-associated signals. Differences in the abundance of these OTU groups between never smokers and daily smokers were more pronounced in younger than older age strata. The smoking- and age-related OTU groups were associated with several cardiometabolic markers, and exploratory mediation analysis suggested that systolic blood pressure may partly account for the association between the pro-aging microbial pattern and ASCVD risk score. In a male-only sensitivity analysis, the central smoking-aging microbial signal remained directionally consistent. Family- and genus-level robustness analyses also supported persistence of selected associations. **Conclusions:** In this cross-sectional secondary analysis, smoking and aging were associated with overlapping gut microbiome patterns, with stronger smoking-related deviations in younger adults. These findings support future longitudinal studies of microbiome-linked cardiovascular risk in smoking-exposed populations.

## What this repo contains

Fourteen self-contained analysis scripts split across two subdirectories:

- `scripts/main/` — main-paper analyses producing Figures 1-6 and Supp Figs 1-5.
- `scripts/revision/` — peer-review-period sensitivity / robustness analyses producing Supp Figs 6-11 and Supp Tables 9-21. Within `scripts/revision/`, scripts `01`-`04` are the original revision-analysis pipeline; scripts `02b` and `05`-`08` are the post-audit fix scripts (see [Audit history](#audit-history)) and supersede the corresponding sub-pipelines in `01` and `03` for the manuscript-ready results.

| Script | Language | Paper figure / table | Purpose |
| --- | --- | --- | --- |
| [`scripts/main/01_main_maaslin2.R`](scripts/main/01_main_maaslin2.R) | R | Fig 1, 2, 3 OTU lists; Supp Figs 1, 2, 4, 5 | MaAsLin2 OTU-level fits for first-hand smoking (222 OTUs), second-hand smoking (117 OTUs), and age within never-smokers (330 OTUs); plus 9 stratified `Maaslin2()` calls (overall + Q1/Q2/Q3 + pack-year variants) and the violin-plot / Venn-diagram inputs of Fig 2 and the cross-stratum Supp Figs. |
| [`scripts/main/02_lightgbm_smk.ipynb`](scripts/main/02_lightgbm_smk.ipynb) | Python (Jupyter, PyCaret 2.x) | Supp Fig 3 | Original PyCaret LightGBM smoking classifier (AUC ≈ 0.73) with SHAP feature-importance panel. |
| [`scripts/main/03_gsea_disease.R`](scripts/main/03_gsea_disease.R) | R | Fig 5 | GSEA-style enrichment of smoking- and age-related OTU lists across six chronic-disease phenotypes (atherosclerosis, fatty liver, T2DM, hepatic calculus, gout, MetS) with `DESeq2`-derived OTU rankings. |
| [`scripts/main/04_cvrisk_ascvd.R`](scripts/main/04_cvrisk_ascvd.R) | R | Fig 6a | ACC/AHA 2013 PCE 10-year ASCVD risk score per subject; Spearman correlation of |GAI| vs ASCVD. Also defines the canonical main-paper `current_smoker` derived smoker indicator (re-imported by the audit-fix `08_refit_mediation_panel.py`). |
| [`scripts/main/05_mediation_main.R`](scripts/main/05_mediation_main.R) | R | Fig 6b | `mediation::mediate(..., sims = 1000, boot = TRUE)` panel of 6 candidate cardiometabolic mediators (FBG, HDL, LDL, HbA1c, SBP, UA), `set.seed(12345)`. |
| [`scripts/revision/01_revision_analysis.py`](scripts/revision/01_revision_analysis.py) | Python | Supp Fig 6 + Supp Table 9; Supp Tables 13, 16 | Rarefaction / sequencing-depth supplement; family/genus age-side OLS (clean — never had a bug); produces `otu_overlap_from_original_results.tsv` consumed downstream. The male-only quartile and family/genus smoking OLS sub-pipelines have been superseded by scripts `06` and `07`. |
| [`scripts/revision/02_maaslin2_reanalysis.R`](scripts/revision/02_maaslin2_reanalysis.R) | R | Supp Tables 18-19 | MaAsLin2 cross-check of the family- and genus-level smoking and age associations. |
| [`scripts/revision/02b_maaslin2_male_only_otu.R`](scripts/revision/02b_maaslin2_male_only_otu.R) | R | Supp Fig 7 right panel + Supp Table 11 | MaAsLin2 male-only OTU-level re-fit for the 40 direction-concordant overlap OTUs in the same LOG framework as the published full-sample model (Bug 1 fix). |
| [`scripts/revision/03_mediation_panel.py`](scripts/revision/03_mediation_panel.py) | Python | (superseded by `08`) | Ten-mediator panel under the original `gai_med.R` LM specification with bootstrap percentile 95% CIs. Superseded by script `08` (Bug 4 fix). |
| [`scripts/revision/04_lightgbm_optimized.py`](scripts/revision/04_lightgbm_optimized.py) | Python | Supp Fig 11 + Supp Table 21 | Optuna TPE LightGBM smoking classifier (200-trial budget) with full per-fold diagnostics. |
| [`scripts/revision/05_replot_male_only_figure.py`](scripts/revision/05_replot_male_only_figure.py) | Python | Supp Fig 7 right panel re-render | Re-render the right panel of Supp Fig 7, consuming the `02b` output. |
| [`scripts/revision/06_refit_male_quartile_module.py`](scripts/revision/06_refit_male_quartile_module.py) | Python | Supp Fig 7 left panel + Supp Table 10 | Re-fit the male-only per-quartile pro-aging-module regression with the corrected sample filter (Bug 2 fix). |
| [`scripts/revision/07_refit_family_genus_smoking_ols.py`](scripts/revision/07_refit_family_genus_smoking_ols.py) | Python | Supp Figs 8-9 + Supp Tables 12, 14, 15, 17 | Re-fit family- and genus-level OLS smoking models with the corrected sample filter (Bug 3 fix); rejoins the (untouched) age-side TSVs. |
| [`scripts/revision/08_refit_mediation_panel.py`](scripts/revision/08_refit_mediation_panel.py) | Python | Supp Fig 10 + Supp Table 20 | Re-fit the multi-mediator panel with the correct ASCVD smoker mapping per `cvrisk.R` (Bug 4 fix). |

A reference set of small output TSVs from our own run is provided under [`outputs-reference/`](outputs-reference/) so that re-runs can be cross-checked.

> **A note on Figure 4 (correlations heatmap).** The published Figure 4 (correlation heatmap of smoking- and age-related OTUs) is built from a separate visualization step in the corresponding author's working tree that consumes the OTU lists produced by `scripts/main/01_main_maaslin2.R` (`non_mach_abx_GPf.rel.smk_status3.mas2`, etc.). That visualization step is not redistributed in this minimal A1 migration; it is straightforward to reconstruct from the OTU lists using any standard correlation-heatmap implementation (`corrplot`, `pheatmap`, etc.).

## Data dependencies

The scripts require three inputs. **None of them are redistributed in this repository.**

1. **Public GGMP OTU table — `GGMP7009_even10k.biom`.**
   Obtain from the published GGMP processing repository:
   <https://github.com/SMUJYYXB/GGMP-Regional-variations>
   Place the file at `data/GGMP7009_even10k.biom`.

2. **Raw 16S FASTQ — ENA accession `PRJEB18535`.**
   <https://www.ebi.ac.uk/ena/browser/view/PRJEB18535>
   *Only required if you wish to reproduce the OTU table from raw reads.* The scripts in this repository all start from the public BIOM table (item 1) and do **not** reprocess raw reads.

3. **Project-specific metadata — `GPf_metadata.tsv`.**
   This is the filtered 6,676-sample metadata subset used in Bao et al. It is derived from the GGMP raw metadata (Supplementary Table S13 of He et al., 2018, *Nature Medicine*) plus project-specific filtering (antibiotic-use exclusion, complete-case retention for the variables used in the present manuscript).
   This file is **not redistributed in this repository** because the derivations are part of the published manuscript's authorship contribution. Users should obtain it from the corresponding author or follow the GGMP data-access procedure. See [`data/README.md`](data/README.md) for the full list of columns the scripts depend on.
   Place the file at `data/GPf_metadata.tsv`.

`scripts/revision/01_revision_analysis.py` and `scripts/revision/04_lightgbm_optimized.py` additionally read the OTU-level smoking and age association tables (`smk_status_sig_res.tsv` / `age_sig_res.tsv`, plus the subordinate `smk_amount_categ_sig_res.tsv` / `smk_y_categ_sig_res.tsv`). These are part of the manuscript's Supplementary Tables and can be regenerated from the main MaAsLin2 calls in the original GGMP repository. Place them at `data/smk_status_sig_res.tsv`, `data/age_sig_res.tsv`, etc.

## Software dependencies

- **Python ≥ 3.9** with the packages listed in [`requirements.txt`](requirements.txt) (notably `biom-format`, `pandas`, `scikit-bio`, `statsmodels`, `scikit-learn`, `lightgbm`, `optuna`).
- **R ≥ 4.2** with the packages listed in [`requirements-r.txt`](requirements-r.txt) (Bioconductor packages `Maaslin2` and `biomformat`; CRAN packages `data.table` and `jsonlite`).

Install Python deps with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Install R deps with:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Maaslin2", "biomformat"))
install.packages(c("data.table", "jsonlite"))
```

## Quick start

```bash
git clone https://github.com/zwbao/ggmp-smoking-aging-analysis.git
cd ggmp-smoking-aging-analysis

# 1. Set up the Python environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Place data files (see Data dependencies above) under ./data/:
#       data/GGMP7009_even10k.biom
#       data/GPf_metadata.tsv
#       data/smk_status_sig_res.tsv
#       data/age_sig_res.tsv
#       data/smk_amount_categ_sig_res.tsv
#       data/smk_y_categ_sig_res.tsv

# 3. Run scripts in order. All revision-period scripts default to relative paths
#    and write to ./outputs/. Main-paper R scripts assume an analysis-ready
#    phyloseq object `non_mach_abx_GPf.rel` is already in the R session
#    (see scripts/main/01_main_maaslin2.R header for construction).

# --- main-paper analyses (Figures 1-6 + Supp Fig 3) -----------------------
Rscript scripts/main/01_main_maaslin2.R            # ~30 min  (main MaAsLin2)
jupyter notebook scripts/main/02_lightgbm_smk.ipynb # interactive (Supp Fig 3)
Rscript scripts/main/03_gsea_disease.R             # ~10 min  (Fig 5)
Rscript scripts/main/04_cvrisk_ascvd.R             # ~1 min   (Fig 6a)
Rscript scripts/main/05_mediation_main.R           # ~10 min  (Fig 6b)

# --- revision-period sensitivity / robustness analyses --------------------
# Scripts 01-04 are the original pipeline. Scripts 02b and 05-08 are the post-audit
# fix scripts and produce the manuscript-ready outputs (see Audit history).
python  scripts/revision/01_revision_analysis.py            # ~5 minutes
Rscript scripts/revision/02_maaslin2_reanalysis.R           # ~10 minutes
Rscript scripts/revision/02b_maaslin2_male_only_otu.R       # ~5 minutes  (Bug 1 fix)
python  scripts/revision/03_mediation_panel.py              # ~10 minutes (superseded by 08)
python  scripts/revision/04_lightgbm_optimized.py           # ~30 minutes
python  scripts/revision/05_replot_male_only_figure.py      # ~10 seconds (Bug 1 fix figure)
python  scripts/revision/06_refit_male_quartile_module.py   # ~1 minute   (Bug 2 fix)
python  scripts/revision/07_refit_family_genus_smoking_ols.py  # ~2 minutes (Bug 3 fix)
python  scripts/revision/08_refit_mediation_panel.py        # ~10 minutes (Bug 4 fix)
```

Total wall-clock runtime on a modern laptop (M-series Mac, 2026): ~1.2 hours. See [`docs/reproducibility.md`](docs/reproducibility.md) for a full step-by-step reproducibility guide.

## Mapping: script → paper figure / table

### Main-paper analyses

| Output | Source script |
| --- | --- |
| Fig 1, 2, 3 OTU lists; Supp Figs 1, 2, 4, 5 (cross-stratum overlap and abundance trends) | `scripts/main/01_main_maaslin2.R` |
| Supp Fig 3 (LightGBM smoking classifier ROC + SHAP) | `scripts/main/02_lightgbm_smk.ipynb` |
| Fig 5 (GSEA enrichment across chronic-disease phenotypes) | `scripts/main/03_gsea_disease.R` |
| Fig 6a (\|GAI\| vs ASCVD risk) | `scripts/main/04_cvrisk_ascvd.R` |
| Fig 6b (mediation diagram) | `scripts/main/05_mediation_main.R` |
| Fig 4 (correlations heatmap) | bespoke visualization step in the corresponding author's tree; consumes OTU lists from `scripts/main/01_main_maaslin2.R` (see "What this repo contains" note above) |

### Revision-period analyses

| Output | Source script |
| --- | --- |
| Supp Fig 6 + Supp Table 9 (rarefaction / sequencing depth) | `scripts/revision/01_revision_analysis.py` |
| Supp Fig 7 left panel + Supp Table 10 (male-only quartile sensitivity) | `scripts/revision/06_refit_male_quartile_module.py` (corrected; supersedes the 01 sub-pipeline) |
| Supp Fig 7 right panel + Supp Table 11 (male-only OTU concordance) | `scripts/revision/02b_maaslin2_male_only_otu.R` + `scripts/revision/05_replot_male_only_figure.py` (corrected) |
| Supp Fig 8 + Supp Tables 12, 14 (family-level smoking + shared) | `scripts/revision/07_refit_family_genus_smoking_ols.py` (corrected; supersedes the 01 sub-pipeline) |
| Supp Table 13 (family-level age within never-smokers) | `scripts/revision/01_revision_analysis.py` (clean — never had a bug) |
| Supp Fig 9 + Supp Tables 15, 17 (genus-level smoking + shared) | `scripts/revision/07_refit_family_genus_smoking_ols.py` (corrected) |
| Supp Table 16 (genus-level age within never-smokers) | `scripts/revision/01_revision_analysis.py` (clean) |
| Supp Tables 18-19 (MaAsLin2 family/genus cross-check) | `scripts/revision/02_maaslin2_reanalysis.R` |
| Supp Fig 10 + Supp Table 20 (multi-mediator panel) | `scripts/revision/08_refit_mediation_panel.py` (corrected; supersedes `03_mediation_panel.py`) |
| Supp Fig 11 + Supp Table 21 (Optuna LightGBM) | `scripts/revision/04_lightgbm_optimized.py` |

A more detailed mapping (including expected file names under `outputs/`) is in [`docs/reproducibility.md`](docs/reproducibility.md).

## Citation

If you use this code, please cite:

> Bao Z, Yang Z, Sun R, Meng R, Wu W, Li MD. **Associations of smoking, aging, and their interplay with the gut microbiome and chronic disease risk profiles.** *Nicotine & Tobacco Research* (2026). DOI: TBD.

## License

This code is released under the [MIT License](LICENSE). The GGMP data are subject to their own licenses and access procedures (see [Data dependencies](#data-dependencies)).

## Contact

**Ming D. Li, PhD** (corresponding author)
State Key Laboratory for Diagnosis and Treatment of Infectious Diseases,
The First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou, China.
Email: ml2km@zju.edu.cn

For questions about the code specifically, please open a GitHub issue.
