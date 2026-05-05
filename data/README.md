# Data dependencies

This directory holds the input data files. **None of them are tracked in git** (see the project `.gitignore`); cloners must obtain them as described below before running any of the scripts.

## 1. Public GGMP OTU table — `GGMP7009_even10k.biom`

- **Source:** the published GGMP processing repository
  <https://github.com/SMUJYYXB/GGMP-Regional-variations>
- **Description:** rarefied (even depth = 10,000 reads) OTU table covering the GGMP cohort. The analyzed subset in this manuscript is restricted to 6,676 samples that intersect with the project-specific metadata (item 3).
- **Place at:** `data/GGMP7009_even10k.biom`
- **Used by:** revision-period scripts `01_revision_analysis.py`, `02_maaslin2_reanalysis.R`, `03_mediation_panel.py`, `04_lightgbm_optimized.py`. Main-paper scripts (`scripts/main/01_main_maaslin2.R` and downstream) consume this file indirectly via the analysis-ready phyloseq object `non_mach_abx_GPf.rel`; see the next section.

## 2. Raw 16S FASTQ — ENA `PRJEB18535`

- **Source:** European Nucleotide Archive
  <https://www.ebi.ac.uk/ena/browser/view/PRJEB18535>
- **Used by:** *none of the scripts in this repository.* This entry is documented for completeness only — the analyses here all start from the public BIOM table (item 1) and do **not** reprocess raw reads. If you wish to reproduce the OTU table from raw FASTQs, follow the original GGMP processing workflow at <https://github.com/SMUJYYXB/GGMP-Regional-variations>.

## 3. Project-specific metadata — `GPf_metadata.tsv`

- **Source:** derived from the GGMP raw metadata (Supplementary Table S13 of He et al., 2018, *Nature Medicine*) plus project-specific filtering used by Bao et al. (n = 6,676 after antibiotic-use exclusion and complete-case filtering for the variables used here).
- **Important: this file is NOT redistributed in this repository.** The derivations are part of the published manuscript's authorship contribution. Users should obtain it from the corresponding author (see the project `README.md` Contact section) or follow the GGMP data-access procedure documented in the original GGMP publication.
- **Place at:** `data/GPf_metadata.tsv`
- **Used by:** revision-period scripts `01_revision_analysis.py`, `02_maaslin2_reanalysis.R`, `03_mediation_panel.py`, `04_lightgbm_optimized.py`. Main-paper scripts consume this file indirectly via the analysis-ready phyloseq object `non_mach_abx_GPf.rel`; see the next section.

### Analysis-ready phyloseq object — `non_mach_abx_GPf.rel`

Some main-paper scripts (`scripts/main/01_main_maaslin2.R`, `03_gsea_disease.R`, `04_cvrisk_ascvd.R`, `05_mediation_main.R`) expect a pre-constructed analysis-ready phyloseq object stored as a `.rda` file (`non_mach_abx_GPf.rel`). This object's construction is documented inside `scripts/main/01_main_maaslin2.R` (lines 19-110) and depends on:

1. The GGMP raw-processing pipeline at <https://github.com/SMUJYYXB/GGMP-Regional-variations>, used to build the relative-abundance phyloseq object `GPf.rel` from `GGMP7009_even10k.biom` plus `GPf_metadata.tsv`.
2. The project-specific antibiotic-use exclusion described in the manuscript Methods (`subset_samples(GPf.rel, antibiotics == "n")`).
3. The derived-variable construction (collapsed `smk_status3`, `age_categ`, `packyear`, `liquor`, `cvd`, `tumor`, `all_dis`, `health`) coded in lines 19-110 of `01_main_maaslin2.R`.

The intermediate `.rda` files used by the corresponding author (`220914/221008.rda`, `221019.rda`, `221108.rda`, etc.; ~1.5-2 GB each) are NOT redistributed here because they bundle too much intermediate session state. Reproduction users should follow the construction steps above; alternatively, request the relevant `.rda` from the corresponding author.

### Key columns the scripts depend on

The scripts read a tab-separated file with the following columns. Numeric columns are coerced via `pd.to_numeric(..., errors="coerce")` / `as.numeric(...)` so missing values can be encoded as `NA`, blanks, or any non-numeric string.

**Identifiers and demographics**

- `ID` — sample identifier matching the BIOM table column IDs.
- `age` — age in years (numeric).
- `gender` — `m` / `f`.
- `smk_status` — categorical smoking status. The raw GGMP values are `everyday`, `not_everyday`, `former_smoker`, and `never_smoker`. The label `current_smoker` is **not** a raw value — it is derived in the original GGMP `code.R` / `cvrisk.R` by collapsing `everyday` ∪ `not_everyday` → `current_smoker`. See [`docs/audit/final_bug_fixes_memo.md`](../docs/audit/final_bug_fixes_memo.md) (Bug 4) for why this distinction matters for the ASCVD mediation panel.
- `smk_y` — smoking pack-year category (used by the OTU-level main analysis).
- `smk_amount_*` — smoking-amount categorical encodings (used by the OTU-level main analysis).
- `smk_sec_*` — second-hand-smoke exposure encodings (used by the OTU-level main analysis).

**Anthropometrics**

- `anthrop_BMI` — body-mass index.
- `anthrop_SBP`, `anthrop_DBP` — systolic / diastolic blood pressure.
- `anthrop_waist` — waist circumference.

**Biochemistry**

- `biochem_FBG` — fasting blood glucose.
- `biochem_TCHO` — total cholesterol.
- `biochem_TG` — triglycerides.
- `biochem_UA` — uric acid.
- `biochem_HDL`, `biochem_LDL` — high- / low-density lipoprotein.

**Stool / region / age stratum**

- `Bristol_stool_type` — Bristol stool scale (1-7).
- `county_level_code` — geographic stratum code.
- `age_categ` — age quartile labelled as `Quantile 1` / `Quantile 2` / `Quantile 3` / `Quantile 4`.

**ASCVD-related fields (used by `03_mediation_panel.py`)**

- `bp_control_medication1`, `bp_control_medication2`, `bp_control_medication3` — blood-pressure-control medication flags (`y` / `n`); any `y` triggers `bp_med = TRUE` in the ACC/AHA 2013 PCE risk equation.
- `bg_diagnosis` — diabetes diagnosis flag (`y` / `n`); used as the diabetes indicator in the PCE risk equation.

If a precomputed `ascvd_10y` (or any column containing `ascvd` / `cvrisk` / `risk_10y` in its name) is present in `GPf_metadata.tsv` and has > 100 numeric values, `03_mediation_panel.py` will reuse it directly; otherwise the script computes 10-year ASCVD risk in-script via a port of `CVrisk::ascvd_10y_accaha`.

## 4. Existing OTU-level association tables (regenerable)

`scripts/revision/01_revision_analysis.py` additionally reads the OTU-level smoking and age association tables produced by the main MaAsLin2 analyses in the original GGMP repository (these tables are themselves the published result of `scripts/main/01_main_maaslin2.R`):

- `data/smk_status_sig_res.tsv` — OTU-level results for the everyday-vs-never smoking model (used as the smoking arm of the OTU overlap).
- `data/age_sig_res.tsv` — OTU-level results for the age model within never smokers.
- `data/smk_amount_categ_sig_res.tsv` — OTU-level results for the categorical pack-year exposure model.
- `data/smk_y_categ_sig_res.tsv` — OTU-level results for the categorical years-smoking exposure model.

These are part of the manuscript's Supplementary Tables. They can be regenerated from the main MaAsLin2 calls in the original GGMP processing repository (<https://github.com/SMUJYYXB/GGMP-Regional-variations>) using the same MaAsLin2 settings documented in the manuscript's Methods.

`scripts/revision/04_lightgbm_optimized.py` reads only `data/smk_status_sig_res.tsv` (to pick the top-71 / top-150 / q<0.05 OTU subsets for the variant screen).

## 5. Main-paper LightGBM feature matrix (`predict_smk.txt`)

`scripts/main/02_lightgbm_smk.ipynb` reads `data/predict_smk.txt` — a tab-separated feature matrix (7,009 rows x 72 columns; 71 OTU relative-abundance feature columns plus the binary smoker `Districts` target column). This file is regenerable from `data/GGMP7009_even10k.biom` + `data/GPf_metadata.tsv` + `data/smk_status_sig_res.tsv` per the OTU-selection logic documented in the manuscript Methods, and is NOT redistributed here.

## Summary checklist

Before running any script, verify that the following files exist:

```
data/GGMP7009_even10k.biom
data/GPf_metadata.tsv
data/smk_status_sig_res.tsv
data/age_sig_res.tsv
data/smk_amount_categ_sig_res.tsv
data/smk_y_categ_sig_res.tsv
```
