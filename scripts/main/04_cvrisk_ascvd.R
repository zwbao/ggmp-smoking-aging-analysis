# ============================================================================
# 04_cvrisk_ascvd.R
#
# ASCVD 10-year cardiovascular risk score per ACC/AHA 2013 Pooled Cohort
# Equations (PCE), and association with the absolute gut-microbiome age gap
# (|GAI|). Produces the data underlying **Figure 6a** (violin / scatter
# panels of |GAI| vs ASCVD risk score) of Bao et al. (NTR-2026-092).
#
# This script is migrated as-is from `220914/cvrisk.R`.
#
# This script also defines the canonical main-paper `current_smoker` derived
# smoker indicator, which is the same definition that the audit-fix
# `scripts/revision/08_refit_mediation_panel.py` re-imports for the
# corrected mediation panel (Bug 4 fix).
#
# ----------------------------------------------------------------------------
# Inputs (expected — NOT redistributed)
# ----------------------------------------------------------------------------
#
# This script expects the following objects already present in the R session:
#
#   - `non_mach_abx_GPf.rel.meta` : analysis-ready metadata data.frame from
#     the phyloseq object built by `01_main_maaslin2.R`. Required columns:
#     `gender` (m/f), `age`, `biochem_TCHO`, `biochem_HDL`, `anthrop_SBP`,
#     `bp_control_medication1`/`2`/`3` (y/n), `smk_status3`
#     (already collapsed to `current_smoker`/`never_smoker`/`former_smoker`),
#     `bg_diagnosis` (y/n), `anthrop_BMI`, `abs_age_gap_adjust` (the absolute
#     gut-microbiome age gap from the LightGBM age-prediction model used
#     elsewhere in the manuscript), `age_categ`.
#
# Construction is documented in the header of `01_main_maaslin2.R`.
#
# ----------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------
#
#   - `cvrisk_df`  : per-subject ASCVD inputs and computed `ascvd` 10-year
#                    risk score. Filtered to `smk_status %in% c("never_smoker",
#                    "current_smoker")` (excludes `former_smoker`).
#   - `med_ana_dat2` : metadata data.frame for the ASCVD-eligible subset,
#                      augmented with the computed `ascvd` column. This frame
#                      is the input to the main-paper mediation analysis (see
#                      `05_mediation_main.R`) and is the published n=4,642
#                      complete-cases frame after further filtering on
#                      cardiometabolic mediator availability.
#   - `ascvd_gai_p1` : ggscatter ggplot object of |GAI| vs ASCVD score (Figure 6a).
#   - PowerPoint file `./pic/h_predicted_age_p1.pptx` (interactive export).
#
# ----------------------------------------------------------------------------
# How to run
# ----------------------------------------------------------------------------
#
#   # interactive (recommended)
#   load("/path/to/221108.rda")  # or run scripts/main/01_main_maaslin2.R first
#   source("scripts/main/04_cvrisk_ascvd.R")
#
# ----------------------------------------------------------------------------
# Reproducibility notes
# ----------------------------------------------------------------------------
#
# - `CVrisk::ascvd_10y_accaha` uses a deterministic closed-form equation;
#   no random seed is needed.
# - Total cholesterol and HDL are converted from mmol/L to mg/dL via
#   the `* 38.67` factor expected by the PCE.
# - `bp_med` falls back to `FALSE` for samples missing all three
#   BP-medication flags (lines 13-14).
# ============================================================================

# install.packages("CVrisk")
library(CVrisk)
library(ggpubr)
library(ggstatsplot)
library(export)

# ---------------------------------------------------------------------------
# Step 1. Build the ASCVD input frame.
# ---------------------------------------------------------------------------

cvrisk_df <- non_mach_abx_GPf.rel.meta[
  c("gender", "age", "biochem_TCHO",
    "biochem_HDL", "anthrop_SBP",
    "bp_control_medication1",
    "bp_control_medication2",
    "bp_control_medication3",
    "smk_status3", "bg_diagnosis", "anthrop_BMI")
]

# bp_med = TRUE if any of the three control-medication flags is "y".
cvrisk_df$bp_control_medication <- ifelse(
  cvrisk_df$bp_control_medication1 == "y" |
    cvrisk_df$bp_control_medication2 == "y" |
    cvrisk_df$bp_control_medication3 == "y",
  "y", "n"
)
cvrisk_df$bp_med <- ifelse(cvrisk_df$bp_control_medication == "y", TRUE, FALSE)
cvrisk_df[is.na(cvrisk_df$bp_med), ]$bp_med <- FALSE

# Map gender labels to the values expected by `ascvd_10y_accaha`.
cvrisk_df[cvrisk_df$gender == "f", ]$gender <- "female"
cvrisk_df[cvrisk_df$gender == "m", ]$gender <- "male"

# Restrict to the canonical main-paper smoker definition: current vs never.
# `current_smoker` here is the *derived* category from `smk_status3` produced
# by `01_main_maaslin2.R` (everyday + not_everyday → current_smoker).
# `former_smoker` is excluded; this is the same filter the audit-fix
# `scripts/revision/08_refit_mediation_panel.py` adopts for Bug 4.
cvrisk_df <- cvrisk_df[cvrisk_df$smk_status %in% c("never_smoker", "current_smoker"), ]
cvrisk_df$smk_status <- factor(cvrisk_df$smk_status,
                               levels = c("never_smoker", "current_smoker"))

cvrisk_df$smk <- ifelse(cvrisk_df$smk_status == "never_smoker", 0, 1)
cvrisk_df$bg <- ifelse(cvrisk_df$bg_diagnosis == "n", 0, 1)

# Convert mmol/L to mg/dL for the ACC/AHA 2013 PCE.
cvrisk_df$biochem_TCHO2 <- cvrisk_df$biochem_TCHO * 38.67
cvrisk_df$biochem_HDL2  <- cvrisk_df$biochem_HDL  * 38.67

# ---------------------------------------------------------------------------
# Step 2. Smoke test on row 1 (matches the original script's sanity check).
# ---------------------------------------------------------------------------

ascvd_10y_accaha(
  race = "other",
  gender = cvrisk_df$gender[1],
  age = cvrisk_df$age[1],
  totchol = cvrisk_df$biochem_TCHO2[1],
  hdl = cvrisk_df$biochem_HDL2[1],
  sbp = cvrisk_df$anthrop_SBP[1],
  bp_med = cvrisk_df$bp_med[1],
  smoker = cvrisk_df$smk[1],
  diabetes = cvrisk_df$bg[1]
)

# ---------------------------------------------------------------------------
# Step 3. Vectorized ASCVD computation across the full eligible sample.
# ---------------------------------------------------------------------------

cvrisk_df[is.na(cvrisk_df$bg), ]$bg <- 0

cvrisk_df$ascvd <- ascvd_10y_accaha(
  race = "other",
  gender = cvrisk_df$gender,
  age = cvrisk_df$age,
  totchol = cvrisk_df$biochem_TCHO2,
  hdl = cvrisk_df$biochem_HDL2,
  sbp = cvrisk_df$anthrop_SBP,
  bp_med = cvrisk_df$bp_med,
  smoker = cvrisk_df$smk,
  diabetes = cvrisk_df$bg
)

# ---------------------------------------------------------------------------
# Step 4. Merge ASCVD back into the full metadata frame for downstream use.
# ---------------------------------------------------------------------------

med_ana_dat2 <- non_mach_abx_GPf.rel.meta[row.names(cvrisk_df), ]
med_ana_dat2$ascvd <- cvrisk_df$ascvd

# Spearman correlation of |GAI| with ASCVD score.
cor.test(abs(med_ana_dat2$abs_age_gap_adjust),
         med_ana_dat2$ascvd, method = "spearman")

# ---------------------------------------------------------------------------
# Step 5. Figure 6a — scatter of |GAI| vs ASCVD score.
# ---------------------------------------------------------------------------

(ascvd_gai_p1 <- ggscatter(
  med_ana_dat2,
  x = "ascvd", y = "abs_age_gap_adjust",
  ggtheme = theme_bw(),
  fill = "#237af2", alpha = 0.1, shape = 21, size = 3,
  add = "reg.line",
  add.params = list(color = "black", fill = "lightgray"),
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.coeff.args = list(method = "spearman", label.x = 20, label.sep = "\n")
) +
    labs_pubr() +
    ylab("Abs(GAI)") + xlab("ASCVD risk score"))

# Note: `h_predicted_age_p1` is constructed in the corresponding author's
# working tree (predicted-age scatter for healthy reference), not in this
# script. The `graph2ppt` line is kept for fidelity but will fail unless
# `h_predicted_age_p1` is bound in the session.
graph2ppt(x = h_predicted_age_p1, file = "./pic/h_predicted_age_p1.pptx",
          width = 4, height = 4)

# ASCVD across age strata (sanity / supplementary panel).
ggbetweenstats(
  data = med_ana_dat2,
  x = age_categ,
  y = ascvd
) + xlab("")
