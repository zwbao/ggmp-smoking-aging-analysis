# ============================================================================
# 01_main_maaslin2.R
#
# Main-paper MaAsLin2 OTU-level analyses for Bao et al. (NTR-2026-092).
#
# Produces (in conjunction with downstream visualization steps maintained
# in the corresponding author's working tree):
#   - 222 first-hand smoking OTU list (Figure 1, 2)
#   - 117 second-hand smoke OTU list (Figure 1)
#   - 330 age OTU list within never-smokers (Figure 1, 3)
#   - per-age-stratum smoking and pack-year MaAsLin2 fits (Figure 2)
#   - violin plots of summed / signed / index relative abundance of
#     smoking-related OTUs across age strata (Figure 2)
#   - Venn-diagram inputs of OTU overlap across age strata (Supp Figs 1, 2,
#     4, 5 source data)
#
# This script is migrated as-is from `220914/code.R` (lines 1-501) and is the
# canonical script that produced the OTU-level signature lists used by
# Figures 1-3 of the published manuscript.
#
# ----------------------------------------------------------------------------
# Inputs (expected — NOT redistributed in this repository)
# ----------------------------------------------------------------------------
#
# This script expects an analysis-ready phyloseq object named
#   non_mach_abx_GPf.rel
# already present in the R session. The construction depends on auxiliary
# `.rda` files maintained in the corresponding author's working tree (see
# `220914/221008.rda` and friends, ~1.5-2 GB each) and is NOT redistributed
# here.
#
# To reproduce the construction step:
#   1. Obtain the public GGMP OTU table from the GGMP processing repo
#      <https://github.com/SMUJYYXB/GGMP-Regional-variations>
#      (file: `GGMP7009_even10k.biom`; place at `data/GGMP7009_even10k.biom`).
#   2. Obtain the project-specific filtered metadata `GPf_metadata.tsv` from
#      the corresponding author (see `data/README.md`).
#   3. Build a phyloseq object `GPf.rel` (relative-abundance-transformed)
#      from those two files following the GGMP pipeline.
#   4. Apply the antibiotic-use exclusion and derived-variable construction
#      steps coded in lines 19-110 of this script: subset to
#      `antibiotics == "n"`, derive `smk_status3` (everyday + not_everyday →
#      current_smoker), `age_categ` (Young / Middle-Age / Elderly), `packyear`,
#      `packyear_categ`, `liquor`, `cvd`, `tumor`, `all_dis`, `health`. The
#      resulting phyloseq object is what we call `non_mach_abx_GPf.rel`.
#
# A working alternative is to load the pre-built `.rda` from the corresponding
# author's tree (e.g. `load("220914/221008.rda")`).
#
# ----------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------
#
# This script writes nine MaAsLin2 native run directories under the working
# directory:
#   non_mach_abx_GPf.rel.smk_status3/       (overall smoking)
#   non_mach_abx_GPf.rel.q1.smk_status3/    (Young stratum smoking)
#   non_mach_abx_GPf.rel.q2.smk_status3/    (Middle-Age stratum smoking)
#   non_mach_abx_GPf.rel.q3.smk_status3/    (Elderly stratum smoking)
#   non_mach_abx_c_GPf.rel.packyear/        (overall pack-year, current smokers)
#   non_mach_abx_c_GPf.rel.q1.packyear/     (Young pack-year)
#   non_mach_abx_c_GPf.rel.q2.packyear/     (Middle-Age pack-year)
#   non_mach_abx_c_GPf.rel.q3.packyear/     (Elderly pack-year)
#   non_mach_abx_c_GPf.rel.packyear_categ/  (categorical pack-year)
#
# Each directory contains MaAsLin2's standard files (`all_results.tsv`,
# `significant_results.tsv`, heatmap, etc.).
#
# The script also constructs in-memory `ggplot2` objects (`fig2_change_p1`,
# `fig2_change_p1.up`, `fig2_change_p1.down`, `fig2_change_p1.smk_status3_index`,
# `smk_status3_q_overlap_p1`, `packyear_q_overlap_p1`,
# `packyear_categ_q_overlap_p1`) used to assemble Figure 2 and the related
# Supp Fig 1, 2, 4, 5 panels.
#
# Note: the published Figure 4 (OTU-OTU correlation heatmap) is built by a
# separate bespoke step in the corresponding author's working tree that is
# NOT part of this minimal A1 migration. The OTU lists consumed by that step
# are produced here (see `non_mach_abx_GPf.rel.smk_status3.mas2`,
# `non_mach_abx_GPf.rel.q1.smk_status3.mas2`, etc.). For Supp Fig 3
# (LightGBM smoking classifier) see `02_lightgbm_smk.ipynb`.
#
# ----------------------------------------------------------------------------
# How to run
# ----------------------------------------------------------------------------
#
#   # interactive (recommended): construct `non_mach_abx_GPf.rel` first
#   # (see step 4 above), then:
#   Rscript scripts/main/01_main_maaslin2.R
#
#   # If an .rda holding the analysis-ready environment is available locally:
#   #   load("/path/to/221108.rda")
#   #   source("scripts/main/01_main_maaslin2.R")
#
# This script does not save figures to disk — it only saves the MaAsLin2
# native run directories. Figure rendering is left to the downstream
# visualization step (PowerPoint / vector export from RStudio interactive
# session). If you want to dump the in-memory ggplot objects to disk, do so
# manually after sourcing this script.
#
# ----------------------------------------------------------------------------
# Reproducibility notes
# ----------------------------------------------------------------------------
#
# - MaAsLin2 LM fits are deterministic given the inputs; no random seed is
#   needed.
# - cores=14 reflects the corresponding author's workstation; lower this on
#   smaller machines.
# ============================================================================

library(tidyverse)
library(ggstatsplot)
library(ggpubr)
library(scales)
library(ggsci)
library(patchwork)
library(phyloseq)
library(ggrepel)
library(Maaslin2)
library(microbiome)
library(microbiomeutilities)
library(export)
library(reshape2)
library(ggstatsplot)
library(vegan)
library(ggside)
library(ggVennDiagram)

# ---------------------------------------------------------------------------
# Step 1. Construct the analysis-ready phyloseq object.
# ---------------------------------------------------------------------------
# `GPf.rel` is expected to be a relative-abundance phyloseq object covering
# the full GGMP cohort with project-specific metadata already merged. See
# the header comment for how to construct it.

non_mach_abx_GPf.rel <- subset_samples(GPf.rel, antibiotics == "n")
non_mach_abx_GPf.rel.otu <- abundances(non_mach_abx_GPf.rel)
non_mach_abx_GPf.rel.meta <- meta(non_mach_abx_GPf.rel)
non_mach_abx_GPf.rel.tax <- data.frame(tax_table(non_mach_abx_GPf.rel))

# Derive `smk_status3`: collapse `everyday` ∪ `not_everyday` → `current_smoker`.
# This is the canonical main-paper smoker definition (see also `cvrisk.R`).
non_mach_abx_GPf.rel.meta$smk_status3 <- str_replace_all(
  non_mach_abx_GPf.rel.meta$smk_status,
  c("not_everyday" = "current_smoker", "everyday" = "current_smoker")
)

non_mach_abx_GPf.rel.meta$age_categ <- cut(
  non_mach_abx_GPf.rel.meta$age,
  breaks = c(-Inf, 39, 59, Inf),
  labels = c("Young", "Middle-Age", "Elderly")
)

non_mach_abx_GPf.rel.meta$tmp <- 0
non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$smk_amount == 1 & !is.na(non_mach_abx_GPf.rel.meta$smk_amount_day), ]$tmp <-
  non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$smk_amount == 1 & !is.na(non_mach_abx_GPf.rel.meta$smk_amount_day), ]$smk_amount_day

non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$smk_amount == 2 & !is.na(non_mach_abx_GPf.rel.meta$smk_amount_week), ]$tmp <-
  as.numeric(non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$smk_amount == 2 & !is.na(non_mach_abx_GPf.rel.meta$smk_amount_week), ]$smk_amount_week) / 7

non_mach_abx_GPf.rel.meta$packyear <- as.numeric(non_mach_abx_GPf.rel.meta$tmp) * as.numeric(non_mach_abx_GPf.rel.meta$smk_y) / 20

non_mach_abx_GPf.rel.meta$packyear_categ <- cut(
  non_mach_abx_GPf.rel.meta$packyear,
  breaks = c(-Inf, 5, 10, 20, Inf),
  labels = c("<5", "5-10", "10-20", ">20"),
  right = FALSE
)

non_mach_abx_GPf.rel.meta$smk_quit_year2 <- non_mach_abx_GPf.rel.meta$smk_quit_days / 365

non_mach_abx_GPf.rel.meta$liquor <- ifelse(
  non_mach_abx_GPf.rel.meta$high_alcohol_liquor > 0 |
    non_mach_abx_GPf.rel.meta$low_alcohol_liquor > 0 |
    non_mach_abx_GPf.rel.meta$beer > 0 |
    non_mach_abx_GPf.rel.meta$yellow_rice_wine > 0 |
    non_mach_abx_GPf.rel.meta$rice_wine > 0 |
    non_mach_abx_GPf.rel.meta$wine > 0,
  "y", "n"
)

non_mach_abx_GPf.rel.meta$cvd <- ifelse(
  non_mach_abx_GPf.rel.meta$dis_atherosclerosis == "y" |
    non_mach_abx_GPf.rel.meta$heart_myocardial_infarction == "y" |
    non_mach_abx_GPf.rel.meta$heart_atrial_fibrillation == "y" |
    non_mach_abx_GPf.rel.meta$heart_angina_pectoris == "y" |
    non_mach_abx_GPf.rel.meta$stroke_ischemic == "y" |
    non_mach_abx_GPf.rel.meta$stroke_hemorrhagic == "y",
  "y", "n"
)

non_mach_abx_GPf.rel.meta$tumor <- ifelse(
  is.na(non_mach_abx_GPf.rel.meta$malignant_tumor_disease) |
    non_mach_abx_GPf.rel.meta$malignant_tumor_disease == "a",
  "n", "y"
)

# `all_dis` flags any of ~40 chronic-disease columns; `health` then filters
# down to the analysis-ready healthy reference subset.
non_mach_abx_GPf.rel.meta$all_dis <- ifelse(
  non_mach_abx_GPf.rel.meta$dis_T1DM == "y" |
    non_mach_abx_GPf.rel.meta$dis_T2DM == "y" |
    non_mach_abx_GPf.rel.meta$dis_fatty_liver == "y" |
    non_mach_abx_GPf.rel.meta$dis_psoriasis == "y" |
    non_mach_abx_GPf.rel.meta$dis_AD == "y" |
    non_mach_abx_GPf.rel.meta$dis_PD == "y" |
    non_mach_abx_GPf.rel.meta$dis_ASD == "y" |
    non_mach_abx_GPf.rel.meta$dis_MS == "y" |
    non_mach_abx_GPf.rel.meta$dis_atherosclerosis == "y" |
    non_mach_abx_GPf.rel.meta$dis_LE == "y" |
    non_mach_abx_GPf.rel.meta$dis_ARDS == "y" |
    non_mach_abx_GPf.rel.meta$dis_gastritis == "y" |
    non_mach_abx_GPf.rel.meta$dis_hepatic_calculus == "y" |
    non_mach_abx_GPf.rel.meta$dis_cholecystitis == "y" |
    non_mach_abx_GPf.rel.meta$dis_colitis == "y" |
    non_mach_abx_GPf.rel.meta$dis_IBS == "y" |
    non_mach_abx_GPf.rel.meta$dis_kidneyStone == "y" |
    non_mach_abx_GPf.rel.meta$dis_gout == "y" |
    non_mach_abx_GPf.rel.meta$dis_AS == "y" |
    non_mach_abx_GPf.rel.meta$dis_RA == "y" |
    non_mach_abx_GPf.rel.meta$dis_neurosis == "y" |
    non_mach_abx_GPf.rel.meta$dis_CFS == "y" |
    non_mach_abx_GPf.rel.meta$dis_constipation_symptom == "y" |
    non_mach_abx_GPf.rel.meta$dis_constipation_days == "y" |
    non_mach_abx_GPf.rel.meta$dis_diarrhea_symptom == "y" |
    non_mach_abx_GPf.rel.meta$dis_diarrhea_days == "y" |
    non_mach_abx_GPf.rel.meta$MetS == "y" |
    non_mach_abx_GPf.rel.meta$heart_myocardial_infarction == "y" |
    non_mach_abx_GPf.rel.meta$heart_atrial_fibrillation == "y" |
    non_mach_abx_GPf.rel.meta$heart_bypass_surgery == "y" |
    non_mach_abx_GPf.rel.meta$heart_stent_surgery == "y" |
    non_mach_abx_GPf.rel.meta$heart_angina_pectoris == "y" |
    non_mach_abx_GPf.rel.meta$stroke_ischemic == "y" |
    non_mach_abx_GPf.rel.meta$stroke_hemorrhagic == "y" |
    non_mach_abx_GPf.rel.meta$copd == "y" |
    non_mach_abx_GPf.rel.meta$asthma == "y" |
    non_mach_abx_GPf.rel.meta$osteoarticular_disease == "y" |
    non_mach_abx_GPf.rel.meta$waist_neck_disease == "y" |
    non_mach_abx_GPf.rel.meta$digestive_system_disease == "y" |
    non_mach_abx_GPf.rel.meta$urinary_system_disease == "y" |
    non_mach_abx_GPf.rel.meta$tumor == "y",
  "y", "n"
)

non_mach_abx_GPf.rel.meta$health <- ifelse(
  non_mach_abx_GPf.rel.meta$all_dis == "n" &
    non_mach_abx_GPf.rel.meta$biochem_FBG < 6.1 &
    non_mach_abx_GPf.rel.meta$anthrop_BMI < 24,
  "y", "n"
)

non_mach_abx_GPf.rel.meta <- non_mach_abx_GPf.rel.meta[!is.na(non_mach_abx_GPf.rel.meta$health), ]
non_mach_abx_GPf.rel.otu <- non_mach_abx_GPf.rel.otu[, row.names(non_mach_abx_GPf.rel.meta)]
sample_data(non_mach_abx_GPf.rel) <- non_mach_abx_GPf.rel.meta

# ---------------------------------------------------------------------------
# Step 2. MaAsLin2: smk_status3 (overall) — produces 222 first-hand smoking OTUs
# ---------------------------------------------------------------------------

mas <- Maaslin2(
  input_data = non_mach_abx_GPf.rel.otu,
  input_metadata = non_mach_abx_GPf.rel.meta,
  output = "non_mach_abx_GPf.rel.smk_status3",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("smk_status3", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH",
  standardize = FALSE,
  reference = c("smk_status3,never_smoker", "county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14,
  plot_heatmap = TRUE,
  plot_scatter = FALSE
)

non_mach_abx_GPf.rel.smk_status3.mas <- mas$results
non_mach_abx_GPf.rel.smk_status3.mas2 <- non_mach_abx_GPf.rel.smk_status3.mas[
  non_mach_abx_GPf.rel.smk_status3.mas$metadata == "smk_status3" &
    non_mach_abx_GPf.rel.smk_status3.mas$qval < 0.05, ]

# ---------------------------------------------------------------------------
# Step 3. Stratified MaAsLin2 (Q1/Q2/Q3 = Young / Middle-Age / Elderly)
# ---------------------------------------------------------------------------

## q1
non_mach_abx_GPf.rel.q1.meta <- non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$age_categ == "Young", ]
non_mach_abx_GPf.rel.q1.otu <- non_mach_abx_GPf.rel.otu[, row.names(non_mach_abx_GPf.rel.q1.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_GPf.rel.q1.otu,
  input_metadata = non_mach_abx_GPf.rel.q1.meta,
  output = "non_mach_abx_GPf.rel.q1.smk_status3",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("smk_status3", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("smk_status3,never_smoker", "county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_GPf.rel.q1.smk_status3.mas <- mas$results
non_mach_abx_GPf.rel.q1.smk_status3.mas2 <- non_mach_abx_GPf.rel.q1.smk_status3.mas[
  non_mach_abx_GPf.rel.q1.smk_status3.mas$metadata == "smk_status3" &
    non_mach_abx_GPf.rel.q1.smk_status3.mas$qval < 0.05, ]

## q2
non_mach_abx_GPf.rel.q2.meta <- non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$age_categ == "Middle-Age", ]
non_mach_abx_GPf.rel.q2.otu <- non_mach_abx_GPf.rel.otu[, row.names(non_mach_abx_GPf.rel.q2.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_GPf.rel.q2.otu,
  input_metadata = non_mach_abx_GPf.rel.q2.meta,
  output = "non_mach_abx_GPf.rel.q2.smk_status3",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("smk_status3", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("smk_status3,never_smoker", "county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_GPf.rel.q2.smk_status3.mas <- mas$results
non_mach_abx_GPf.rel.q2.smk_status3.mas2 <- non_mach_abx_GPf.rel.q2.smk_status3.mas[
  non_mach_abx_GPf.rel.q2.smk_status3.mas$metadata == "smk_status3" &
    non_mach_abx_GPf.rel.q2.smk_status3.mas$qval < 0.05, ]

## q3
non_mach_abx_GPf.rel.q3.meta <- non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$age_categ == "Elderly", ]
non_mach_abx_GPf.rel.q3.otu <- non_mach_abx_GPf.rel.otu[, row.names(non_mach_abx_GPf.rel.q3.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_GPf.rel.q3.otu,
  input_metadata = non_mach_abx_GPf.rel.q3.meta,
  output = "non_mach_abx_GPf.rel.q3.smk_status3",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("smk_status3", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("smk_status3,never_smoker", "county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_GPf.rel.q3.smk_status3.mas <- mas$results
non_mach_abx_GPf.rel.q3.smk_status3.mas2 <- non_mach_abx_GPf.rel.q3.smk_status3.mas[
  non_mach_abx_GPf.rel.q3.smk_status3.mas$metadata == "smk_status3" &
    non_mach_abx_GPf.rel.q3.smk_status3.mas$qval < 0.05, ]

# ---------------------------------------------------------------------------
# Step 4. MaAsLin2: pack-year (continuous) on current smokers only
# ---------------------------------------------------------------------------

non_mach_abx_c_GPf.rel <- subset_samples(non_mach_abx_GPf.rel, smk_status3 == "current_smoker")
non_mach_abx_c_GPf.rel.otu <- abundances(non_mach_abx_c_GPf.rel)
non_mach_abx_c_GPf.rel.meta <- meta(non_mach_abx_c_GPf.rel)

mas <- Maaslin2(
  input_data = non_mach_abx_c_GPf.rel.otu,
  input_metadata = non_mach_abx_c_GPf.rel.meta,
  output = "non_mach_abx_c_GPf.rel.packyear",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("packyear", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_c_GPf.rel.packyear.mas <- mas$results
non_mach_abx_c_GPf.rel.packyear.mas2 <- non_mach_abx_c_GPf.rel.packyear.mas[
  non_mach_abx_c_GPf.rel.packyear.mas$metadata == "packyear" &
    non_mach_abx_c_GPf.rel.packyear.mas$qval < 0.05, ]

## q1
non_mach_abx_c_GPf.rel.q1.meta <- non_mach_abx_c_GPf.rel.meta[non_mach_abx_c_GPf.rel.meta$age_categ == "Young", ]
non_mach_abx_c_GPf.rel.q1.otu <- non_mach_abx_c_GPf.rel.otu[, row.names(non_mach_abx_c_GPf.rel.q1.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_c_GPf.rel.q1.otu,
  input_metadata = non_mach_abx_c_GPf.rel.q1.meta,
  output = "non_mach_abx_c_GPf.rel.q1.packyear",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("packyear", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_c_GPf.rel.q1.packyear.mas <- mas$results
non_mach_abx_c_GPf.rel.q1.packyear.mas2 <- non_mach_abx_c_GPf.rel.q1.packyear.mas[
  non_mach_abx_c_GPf.rel.q1.packyear.mas$metadata == "packyear" &
    non_mach_abx_c_GPf.rel.q1.packyear.mas$qval < 0.05, ]

## q2
non_mach_abx_c_GPf.rel.q2.meta <- non_mach_abx_c_GPf.rel.meta[non_mach_abx_c_GPf.rel.meta$age_categ == "Middle-Age", ]
non_mach_abx_c_GPf.rel.q2.otu <- non_mach_abx_c_GPf.rel.otu[, row.names(non_mach_abx_c_GPf.rel.q2.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_c_GPf.rel.q2.otu,
  input_metadata = non_mach_abx_c_GPf.rel.q2.meta,
  output = "non_mach_abx_c_GPf.rel.q2.packyear",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("packyear", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_c_GPf.rel.q2.packyear.mas <- mas$results
non_mach_abx_c_GPf.rel.q2.packyear.mas2 <- non_mach_abx_c_GPf.rel.q2.packyear.mas[
  non_mach_abx_c_GPf.rel.q2.packyear.mas$metadata == "packyear" &
    non_mach_abx_c_GPf.rel.q2.packyear.mas$qval < 0.05, ]

## q3
non_mach_abx_c_GPf.rel.q3.meta <- non_mach_abx_c_GPf.rel.meta[non_mach_abx_c_GPf.rel.meta$age_categ == "Elderly", ]
non_mach_abx_c_GPf.rel.q3.otu <- non_mach_abx_c_GPf.rel.otu[, row.names(non_mach_abx_c_GPf.rel.q3.meta)]

mas <- Maaslin2(
  input_data = non_mach_abx_c_GPf.rel.q3.otu,
  input_metadata = non_mach_abx_c_GPf.rel.q3.meta,
  output = "non_mach_abx_c_GPf.rel.q3.packyear",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("packyear", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_c_GPf.rel.q3.packyear.mas <- mas$results
non_mach_abx_c_GPf.rel.q3.packyear.mas2 <- non_mach_abx_c_GPf.rel.q3.packyear.mas[
  non_mach_abx_c_GPf.rel.q3.packyear.mas$metadata == "packyear" &
    non_mach_abx_c_GPf.rel.q3.packyear.mas$qval < 0.05, ]

# ---------------------------------------------------------------------------
# Step 5. MaAsLin2: categorical pack-year exposure
# ---------------------------------------------------------------------------

mas <- Maaslin2(
  input_data = non_mach_abx_c_GPf.rel.otu,
  input_metadata = non_mach_abx_c_GPf.rel.meta,
  output = "non_mach_abx_c_GPf.rel.packyear_categ",
  min_abundance = 0.0, min_prevalence = 0.0,
  normalization = "TSS", transform = "LOG", analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("packyear_categ", "gender", "age", "county_level_code", "Bristol_stool_type"),
  correction = "BH", standardize = FALSE,
  reference = c("packyear_categ,<5", "county_level_code,G440282", "Bristol_stool_type,4"),
  cores = 14, plot_heatmap = TRUE, plot_scatter = FALSE
)

non_mach_abx_c_GPf.rel.packyear_categ.mas <- mas$results
non_mach_abx_c_GPf.rel.packyear_categ.mas2 <- non_mach_abx_c_GPf.rel.packyear_categ.mas[
  non_mach_abx_c_GPf.rel.packyear_categ.mas$metadata == "packyear_categ" &
    non_mach_abx_c_GPf.rel.packyear_categ.mas$qval < 0.05, ]

# ---------------------------------------------------------------------------
# Step 6. Re-load saved MaAsLin2 result tables (from the .tsv files just written)
# ---------------------------------------------------------------------------

non_mach_abx_GPf.rel.smk_status3.mas <- read.table("./non_mach_abx_GPf.rel.smk_status3/all_results.tsv", header = TRUE)
non_mach_abx_GPf.rel.smk_status3.mas2 <- non_mach_abx_GPf.rel.smk_status3.mas[
  non_mach_abx_GPf.rel.smk_status3.mas$metadata == "smk_status3" &
    non_mach_abx_GPf.rel.smk_status3.mas$qval < 0.05, ]

non_mach_abx_c_GPf.rel.packyear.mas <- read.table("./non_mach_abx_c_GPf.rel.packyear/all_results.tsv", header = TRUE)
non_mach_abx_c_GPf.rel.packyear.mas2 <- non_mach_abx_c_GPf.rel.packyear.mas[
  non_mach_abx_c_GPf.rel.packyear.mas$metadata == "packyear" &
    non_mach_abx_c_GPf.rel.packyear.mas$qval < 0.05, ]

non_mach_abx_c_GPf.rel.packyear_categ.mas <- read.table("./non_mach_abx_c_GPf.rel.packyear_categ/all_results.tsv", header = TRUE)
non_mach_abx_c_GPf.rel.packyear_categ.mas2 <- non_mach_abx_c_GPf.rel.packyear_categ.mas[
  non_mach_abx_c_GPf.rel.packyear_categ.mas$metadata == "packyear_categ" &
    non_mach_abx_c_GPf.rel.packyear_categ.mas$qval < 0.05, ]

# ---------------------------------------------------------------------------
# Step 7. Aggregate smoking-related OTU abundance and visualize across age strata
# ---------------------------------------------------------------------------

### all smk_status3-related OTUs

non_mach_abx_GPf.rel.meta$smk_status3_related.sum <- apply(
  non_mach_abx_GPf.rel.otu[str_replace(non_mach_abx_GPf.rel.smk_status3.mas2$feature, "X", ""), ],
  2, sum
)

tmp <- melt(non_mach_abx_GPf.rel.meta,
            id.vars = c("age_categ", "smk_status3"),
            measure.vars = c("smk_status3_related.sum"))
tmp <- tmp[tmp$smk_status %in% c("never_smoker", "current_smoker"), ]

(fig2_change_p1 <- ggviolin(tmp, x = "age_categ", y = "value", fill = "smk_status3",
                            palette = c("#ababab", "#747474"),
                            add = "boxplot", add.params = list(fill = "smk_status3")) +
    stat_compare_means(method = "wilcox.test", aes(group = smk_status3),
                       label = "p.signif", hide.ns = TRUE,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns"))) +
    scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance") + xlab("") +
    theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 2) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 18),
          legend.position = "top",
          legend.title = element_blank()) +
    theme(legend.position = "right"))

### up
non_mach_abx_GPf.rel.meta$smk_status3_related.up.sum <- apply(
  non_mach_abx_GPf.rel.otu[str_replace(
    non_mach_abx_GPf.rel.smk_status3.mas2[non_mach_abx_GPf.rel.smk_status3.mas2$coef > 0, ]$feature,
    "X", ""), ],
  2, sum
)

tmp <- melt(non_mach_abx_GPf.rel.meta,
            id.vars = c("age_categ", "smk_status3"),
            measure.vars = c("smk_status3_related.up.sum"))
tmp <- tmp[tmp$smk_status %in% c("never_smoker", "current_smoker"), ]

(fig2_change_p1.up <- ggviolin(tmp, x = "age_categ", y = "value", fill = "smk_status3",
                               palette = c("#ababab", "#747474"),
                               add = "boxplot", add.params = list(fill = "smk_status3")) +
    stat_compare_means(method = "wilcox.test", aes(group = smk_status3),
                       label = "p.signif", hide.ns = TRUE,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns"))) +
    scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance") + xlab("") +
    theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 2) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 18),
          legend.position = "top",
          legend.title = element_blank()) +
    theme(legend.position = "right"))

### down
non_mach_abx_GPf.rel.meta$smk_status3_related.down.sum <- apply(
  non_mach_abx_GPf.rel.otu[str_replace(
    non_mach_abx_GPf.rel.smk_status3.mas2[non_mach_abx_GPf.rel.smk_status3.mas2$coef < 0, ]$feature,
    "X", ""), ],
  2, sum
)

tmp <- melt(non_mach_abx_GPf.rel.meta,
            id.vars = c("age_categ", "smk_status3"),
            measure.vars = c("smk_status3_related.down.sum"))
tmp <- tmp[tmp$smk_status %in% c("never_smoker", "current_smoker"), ]

(fig2_change_p1.down <- ggviolin(tmp, x = "age_categ", y = "value", fill = "smk_status3",
                                 palette = c("#ababab", "#747474"),
                                 add = "boxplot", add.params = list(fill = "smk_status3")) +
    stat_compare_means(method = "wilcox.test", aes(group = smk_status3),
                       label = "p.signif", hide.ns = TRUE,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns"))) +
    scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance") + xlab("") +
    theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 2) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 18),
          legend.position = "top",
          legend.title = element_blank()) +
    theme(legend.position = "right"))

# ---------------------------------------------------------------------------
# Step 8. Smoking-index (signed, magnitude-weighted aggregate)
# ---------------------------------------------------------------------------

non_mach_abx_GPf.rel.smk_status3.mas3 <- non_mach_abx_GPf.rel.smk_status3.mas2[
  non_mach_abx_GPf.rel.smk_status3.mas2$value == "current_smoker", ]
smk_status3_related.rel.otu <- data.frame(t(non_mach_abx_GPf.rel.otu))[, non_mach_abx_GPf.rel.smk_status3.mas3$feature]

non_mach_abx_GPf.rel.smk_status3.mas3 <- non_mach_abx_GPf.rel.smk_status3.mas2[
  non_mach_abx_GPf.rel.smk_status3.mas2$value == "current_smoker", ]

non_mach_abx_GPf.rel.meta$smk_status3_index <- 0

for (i in 1:dim(non_mach_abx_GPf.rel.meta)[1]) {
  non_mach_abx_GPf.rel.meta[i, "smk_status3_index"] <- sum(
    log10(smk_status3_related.rel.otu[i, ] * 10000 + 1) *
      non_mach_abx_GPf.rel.smk_status3.mas3$coef *
      log10(non_mach_abx_GPf.rel.smk_status3.mas3$qval)
  )
}

tmp <- melt(non_mach_abx_GPf.rel.meta,
            id.vars = c("age_categ", "smk_status3"),
            measure.vars = c("smk_status3_index"))
tmp <- tmp[tmp$smk_status %in% c("never_smoker", "current_smoker"), ]

(fig2_change_p1.smk_status3_index <- ggviolin(tmp, x = "age_categ", y = "value", fill = "smk_status3",
                                              palette = c("#ababab", "#747474"),
                                              add = "boxplot", add.params = list(fill = "smk_status3")) +
    stat_compare_means(method = "wilcox.test", aes(group = smk_status3),
                       label = "p.signif", hide.ns = TRUE,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "ns"))) +
    scale_y_continuous(labels = scales::percent) + ylab("Relative Abundance") + xlab("") +
    theme_bw(base_size = 14, base_line_size = 1, base_rect_size = 2) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(size = 18),
          legend.position = "top",
          legend.title = element_blank()) +
    theme(legend.position = "right"))

ggbetweenstats(
  data = non_mach_abx_GPf.rel.meta[non_mach_abx_GPf.rel.meta$smk_status3 == "never_smoker", ],
  x = age_categ,
  y = smk_status3_index
) + xlab("")

# ---------------------------------------------------------------------------
# Step 9. Cross-stratum overlap (Venn diagrams)
# ---------------------------------------------------------------------------

tmp <- list(q1 = non_mach_abx_GPf.rel.q1.smk_status3.mas2$feature,
            q2 = non_mach_abx_GPf.rel.q2.smk_status3.mas2$feature,
            q3 = non_mach_abx_GPf.rel.q3.smk_status3.mas2$feature)
(smk_status3_q_overlap_p1 <- ggVennDiagram(tmp))

tmp <- list(q1 = non_mach_abx_c_GPf.rel.q1.packyear.mas2$feature,
            q2 = non_mach_abx_c_GPf.rel.q2.packyear.mas2$feature,
            q3 = non_mach_abx_c_GPf.rel.q3.packyear.mas2$feature)
(packyear_q_overlap_p1 <- ggVennDiagram(tmp))

tmp <- list(q1 = non_mach_abx_c_GPf.rel.q1.packyear_categ.mas2$feature,
            q2 = non_mach_abx_c_GPf.rel.q2.packyear_categ.mas2$feature,
            q3 = non_mach_abx_c_GPf.rel.q3.packyear_categ.mas2$feature)
(packyear_categ_q_overlap_p1 <- ggVennDiagram(tmp))

# Optional: persist the workspace for downstream visualization scripts.
# save.image("221108.rda")
