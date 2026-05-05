# Smoking, aging, and the gut microbiome — analysis code for Bao et al., NTR-2026-092

This repository contains the analysis code for the supplementary results of:

> Bao Z, Yang Z, Sun R, Meng R, Wu W, Li MD. **Associations of smoking, aging, and their interplay with the gut microbiome and chronic disease risk profiles.** *Nicotine & Tobacco Research* (2026). DOI: TBD.

It is provided so that anyone with access to the public Guangdong Gut Microbiome Project (GGMP) data can re-run the four secondary / sensitivity analyses introduced during peer review. The original GGMP raw-processing workflow is *not* included here — see [Data dependencies](#data-dependencies).

## Abstract

**Background:** Smoking and aging are both linked to chronic disease and to variation in the gut microbiome, but their joint relationship with microbiome composition remains incompletely characterized. This study examined the overlap between smoking- and age-associated gut microbial signals and their associations with cardiovascular risk profiles. **Methods:** This study was a secondary analysis of the published Guangdong Gut Microbiome Project (GGMP) resource. We analyzed 6,676 participants from a filtered project-specific metadata subset together with the public GGMP OTU table. Associations of smoking phenotypes and age with gut microbiota were examined using multivariable linear models. Machine-learning and mediation analyses were treated as exploratory. Additional sensitivity analyses were performed in males only, and robustness analyses were repeated at the family and genus levels. **Results:** We identified 222 OTUs associated with first-hand smoking and 117 OTUs associated with second-hand smoke exposure. Among never smokers, age was associated with 330 OTUs. Eighty-five OTUs overlapped between smoking- and age-associated signals. Differences in the abundance of these OTU groups between never smokers and daily smokers were more pronounced in younger than older age strata. The smoking- and age-related OTU groups were associated with several cardiometabolic markers, and exploratory mediation analysis suggested that systolic blood pressure may partly account for the association between the pro-aging microbial pattern and ASCVD risk score. In a male-only sensitivity analysis, the central smoking-aging microbial signal remained directionally consistent. Family- and genus-level robustness analyses also supported persistence of selected associations. **Conclusions:** In this cross-sectional secondary analysis, smoking and aging were associated with overlapping gut microbiome patterns, with stronger smoking-related deviations in younger adults. These findings support future longitudinal studies of microbiome-linked cardiovascular risk in smoking-exposed populations.

## What this repo contains

Four self-contained analysis scripts that produce the supplementary figures and tables introduced during peer review:

| Script | Language | Purpose |
| --- | --- | --- |
| [`scripts/01_revision_analysis.py`](scripts/01_revision_analysis.py) | Python | Rarefaction / sequencing-depth supplement, male-only sensitivity analysis on the pro-aging module, and family/genus OLS robustness re-tests. |
| [`scripts/02_maaslin2_reanalysis.R`](scripts/02_maaslin2_reanalysis.R) | R | MaAsLin2 cross-check of the family- and genus-level smoking and age associations. |
| [`scripts/03_mediation_panel.py`](scripts/03_mediation_panel.py) | Python | Ten-mediator multi-mediator panel under the original `gai_med.R` LM specification, with bootstrap percentile 95% CIs. |
| [`scripts/04_lightgbm_optimized.py`](scripts/04_lightgbm_optimized.py) | Python | Optuna TPE LightGBM smoking classifier (200-trial budget) with full per-fold diagnostics. |

A reference set of small output TSVs from our own run is provided under [`outputs-reference/`](outputs-reference/) so that re-runs can be cross-checked.

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

`scripts/01_revision_analysis.py` and `scripts/04_lightgbm_optimized.py` additionally read the OTU-level smoking and age association tables (`smk_status_sig_res.tsv` / `age_sig_res.tsv`, plus the subordinate `smk_amount_categ_sig_res.tsv` / `smk_y_categ_sig_res.tsv`). These are part of the manuscript's Supplementary Tables and can be regenerated from the main MaAsLin2 calls in the original GGMP repository. Place them at `data/smk_status_sig_res.tsv`, `data/age_sig_res.tsv`, etc.

## Software dependencies

- **Python ≥ 3.9** with the packages listed in [`requirements.txt`](requirements.txt) (notably `biom-format`, `pandas`, `scikit-bio`, `statsmodels`, `lightgbm`, `optuna`).
- **R ≥ 4.2** with the packages listed in [`requirements-r.txt`](requirements-r.txt) (notably the Bioconductor packages `Maaslin2`, `biomformat`, and `phyloseq`).

Install Python deps with:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Install R deps with:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Maaslin2", "biomformat", "phyloseq"))
install.packages(c("dplyr", "readr", "data.table"))
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

# 3. Run scripts in order. All scripts default to relative paths and write to ./outputs/.
python scripts/01_revision_analysis.py            # ~5 minutes
Rscript scripts/02_maaslin2_reanalysis.R          # ~10 minutes
python scripts/03_mediation_panel.py              # ~10 minutes
python scripts/04_lightgbm_optimized.py           # ~30 minutes
```

Total wall-clock runtime on a modern laptop (M-series Mac, 2026): ~1 hour. See [`docs/reproducibility.md`](docs/reproducibility.md) for a full step-by-step reproducibility guide.

## Mapping: script → paper figure / table

| Output | Source script |
| --- | --- |
| Supp Fig 6 + Supp Table 9 (rarefaction / sequencing depth) | `01_revision_analysis.py` |
| Supp Fig 7 + Supp Tables 10-11 (male-only sensitivity) | `01_revision_analysis.py` |
| Supp Figs 8-9 + Supp Tables 12-17 (family/genus OLS robustness) | `01_revision_analysis.py` |
| Supp Tables 18-19 (MaAsLin2 family/genus cross-check) | `02_maaslin2_reanalysis.R` |
| Supp Fig 10 + Supp Table 20 (multi-mediator panel) | `03_mediation_panel.py` |
| Supp Fig 11 + Supp Table 21 (Optuna LightGBM) | `04_lightgbm_optimized.py` |

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
