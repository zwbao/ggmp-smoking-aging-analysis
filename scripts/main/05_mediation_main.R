# ============================================================================
# 05_mediation_main.R
#
# Main-paper mediation analysis: |GAI| → cardiometabolic mediators → ASCVD.
# Produces the data underlying **Figure 6b** (mediation diagram) of Bao et al.
# (NTR-2026-092).
#
# This script is migrated as-is from `220914/gai_med.R` and uses
# `mediation::mediate(..., sims = 1000, boot = TRUE)` per the manuscript
# Methods. Six candidate cardiometabolic mediators are screened
# (`biochem_FBG`, `biochem_HDL`, `biochem_LDL`, `biochem_HbA1c`, `anthrop_SBP`,
# `biochem_UA`); the published main-paper mediator is systolic blood pressure
# (`anthrop_SBP`).
#
# **Note on the parallel revision-period script.** The revision-analysis
# `scripts/revision/08_refit_mediation_panel.py` re-runs an expanded
# 10-mediator panel under the same ASCVD smoker mapping, in Python with
# percentile bootstrap CIs, and is the source for Supp Table 20 / Supp Fig
# 10. The present script is the original main-paper specification used to
# generate Figure 6b.
#
# ----------------------------------------------------------------------------
# Inputs (expected — NOT redistributed)
# ----------------------------------------------------------------------------
#
# This script expects the following objects already present in the R session:
#
#   - `med_ana_dat2`     : the ASCVD-augmented metadata frame produced by
#                          `04_cvrisk_ascvd.R`. Required columns:
#                          `abs_age_gap_adjust`, `ascvd`, `biochem_FBG`,
#                          `biochem_HDL`, `biochem_LDL`, `biochem_HbA1c`,
#                          `anthrop_SBP`, `biochem_UA`, plus the lifestyle
#                          variables used in the exploratory SEM
#                          (`fam_sauce_avg`, `act_static_time`, `vegetables`,
#                          `grains`, `fruit_drinks`, `carbonated_beverage`).
#   - `gai_with_other.sig`: data.frame of "GAI vs other variable" Spearman
#                           correlations with FDR-significant rows already
#                           filtered. Used in the lifestyle-biochem screening
#                           step at the top of the script.
#   - `lifestyle_name`, `Biochemistry_name` : character vectors of column
#                           names used by the lifestyle-biochem screening loop.
#   - `gai_plot_df`      : per-subject data.frame indexed by sample with the
#                          lifestyle and biochemistry columns named in
#                          `lifestyle_name` / `Biochemistry_name`.
#
# All of these are constructed in the corresponding author's working tree
# from `non_mach_abx_GPf.rel.meta` plus the LightGBM age-prediction outputs.
# To reproduce, run `scripts/main/01_main_maaslin2.R` and
# `scripts/main/04_cvrisk_ascvd.R` first; then build `gai_plot_df`,
# `lifestyle_name`, `Biochemistry_name`, and `gai_with_other.sig` per the
# manuscript Methods (the full construction is documented in the
# corresponding author's tree).
#
# ----------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------
#
#   - `lifestyle_biochem_df` : screening grid of Spearman correlations
#                              between lifestyle and biochemistry variables
#                              that are individually associated with |GAI|.
#   - `results`              : single mediator-screen `mediation::mediate`
#                              result object (LDL example, kept for fidelity).
#   - `path`                 : `lavaan::sem` SEM fit object (exploratory).
#   - `res`                  : data.frame of ACME / ADE / total-effect /
#                              proportion-mediated for the 6-mediator panel.
#                              Each row is one candidate mediator. This is
#                              the source for the Figure 6b mediation diagram.
#
# ----------------------------------------------------------------------------
# How to run
# ----------------------------------------------------------------------------
#
#   # interactive (recommended)
#   load("/path/to/221108.rda")  # or run scripts 01 and 04 first
#   source("scripts/main/05_mediation_main.R")
#
# ----------------------------------------------------------------------------
# Reproducibility notes
# ----------------------------------------------------------------------------
#
# - Random seed is set explicitly to `12345` before the panel mediation step,
#   per the manuscript Methods. This matches the seed reused by
#   `scripts/revision/08_refit_mediation_panel.py`.
# - Bootstrap resamples: `sims = 1000`, `boot = TRUE`.
# - lavaan::sem is fit with `se = "bootstrap", bootstrap = 1000`.
# ============================================================================

# Set reproducibility seed up-front (the original `gai_med.R` set it later;
# moved up here per the migration spec).
set.seed(12345)

library(mediation)
library(lavaan)
library(semPlot)
library(plyr)

# ---------------------------------------------------------------------------
# Step 0. Lifestyle x biochemistry screening grid.
# ---------------------------------------------------------------------------
# Builds a table of Spearman correlations for every (lifestyle, biochemistry)
# pair where both variables are individually FDR-significant against |GAI|
# (`gai_with_other.sig`).

lifestyle_biochem_df <- data.frame(
  id = 1:(sum(lifestyle_name %in% gai_with_other.sig$id) *
            sum(Biochemistry_name %in% gai_with_other.sig$id)),
  cor_r = 0,
  cor_p = 0
)
n <- 1
for (i in lifestyle_name[lifestyle_name %in% gai_with_other.sig$id]) {
  for (j in Biochemistry_name[Biochemistry_name %in% gai_with_other.sig$id]) {
    lifestyle_biochem_df[n, "id"] <- paste(i, j, sep = "_")
    tmp1 <- cor.test(gai_plot_df[, i], gai_plot_df[, j], method = "spearman")
    lifestyle_biochem_df[n, "cor_r"] <- tmp1$estimate
    lifestyle_biochem_df[n, "cor_p"] <- tmp1$p.value
    n <- n + 1
  }
}

# Variables of interest reminded in the original script (kept for context):
fam_sauce_avg
act_static_time
biochem_LDL
med_ana_dat2
fruit_drinks

# ---------------------------------------------------------------------------
# Step 1. Single-mediator example (LDL) — kept for fidelity with the original.
# ---------------------------------------------------------------------------

tmp <- med_ana_dat2[!is.na(med_ana_dat2$fruit_drinks) &
                      !is.na(med_ana_dat2$biochem_LDL) &
                      !is.na(med_ana_dat2$abs_age_gap_adjust), ]

fit.totaleffect <- lm(biochem_LDL ~ fam_sauce_avg + act_static_time + vegetables + grains, tmp)
fit.mediator   <- lm(abs_age_gap_adjust ~ fam_sauce_avg + act_static_time + vegetables + grains, tmp)
fit.dv         <- lm(biochem_LDL ~ fam_sauce_avg + act_static_time + vegetables + grains + abs_age_gap_adjust, tmp)
results <- mediate(fit.mediator, fit.dv,
                   treat = "vegetables",
                   mediator = "abs_age_gap_adjust", boot = TRUE)
summary(results)

# ---------------------------------------------------------------------------
# Step 2. Exploratory SEM (lavaan).
# ---------------------------------------------------------------------------

tmp <- med_ana_dat2[!is.na(med_ana_dat2$fruit_drinks) &
                      !is.na(med_ana_dat2$biochem_LDL) &
                      !is.na(med_ana_dat2$abs_age_gap_adjust) &
                      !is.na(med_ana_dat2$fam_sauce_avg) &
                      !is.na(med_ana_dat2$act_static_time) &
                      !is.na(med_ana_dat2$vegetables) &
                      !is.na(med_ana_dat2$grains) &
                      !is.na(med_ana_dat2$biochem_HDL), ]

model <- '
abs_age_gap_adjust ~ fam_sauce_avg + act_static_time + vegetables + grains
biochem_LDL ~ abs_age_gap_adjust + fam_sauce_avg + act_static_time + vegetables + grains
biochem_HDL ~ abs_age_gap_adjust + fam_sauce_avg + act_static_time + vegetables + grains
ascvd ~ biochem_LDL + abs_age_gap_adjust + biochem_HDL
'

# Alternative SEM specification (sauce + carbonated_beverage):
# model <- '
# abs_age_gap_adjust ~ fam_sauce_avg + carbonated_beverage
# biochem_LDL ~ abs_age_gap_adjust + fam_sauce_avg + carbonated_beverage
# ascvd ~ biochem_LDL + abs_age_gap_adjust
# '

path <- sem(model = model, data = tmp,
            se = "bootstrap", bootstrap = 1000)
summary(path, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

fitMeasures(path, c("chisq", "df", "pvalue", "gfi", "cfi", "rmr", "srmr", "rmsea"))

semPaths(path, what = "std", layout = "tree",
         residuals = FALSE, edge.label.cex = 1)

# ---------------------------------------------------------------------------
# Step 3. Six-mediator panel (Figure 6b) via `mediation::mediate` bootstrap.
# ---------------------------------------------------------------------------

setRefClass("TS", fields = list(m = "character", x = "character",
                                xm = "formula", xy = "formula"))
test <- new("TS")

batch_mediate <- function(ms, y, xs, df) {
  test$m <- ms
  test$x <- xs
  test$xm <- as.formula(paste0(ms, "~", xs))
  model_xm <- lm(test$xm, data = df)
  test$xy <- as.formula(paste0(y, "~", ms, "+", xs))
  model_xy <- lm(test$xy, data = df)
  model_med <- mediate(model_xm, model_xy,
                       treat = test$x, mediator = test$m,
                       sims = 1000, boot = TRUE)
  med_sum <- summary(model_med)
  ACME <- med_sum$d0
  ACME_P <- med_sum$d0.p
  ADE <- med_sum$z0
  ADE_P <- med_sum$z0.p
  TE <- (ACME + ADE)
  Prop_Mediated <- ACME / TE
  result <- data.frame(
    Characteristics = ms,
    ACME = ACME, ACME_P = ACME_P,
    ADE = ADE, ADE_P = ADE_P,
    Total_Effect = TE,
    Prop_Mediated = Prop_Mediated
  )
  return(result)
}

# Candidate cardiometabolic mediators screened in the main-paper analysis.
cols <- c("biochem_FBG", "biochem_HDL", "biochem_LDL",
          "biochem_HbA1c", "anthrop_SBP", "biochem_UA")

# Reset seed immediately before the panel for full reproducibility.
set.seed(12345)

tmp <- med_ana_dat2[!is.na(med_ana_dat2$biochem_FBG) &
                      !is.na(med_ana_dat2$biochem_LDL) &
                      !is.na(med_ana_dat2$abs_age_gap_adjust) &
                      !is.na(med_ana_dat2$biochem_HbA1c) &
                      !is.na(med_ana_dat2$anthrop_SBP) &
                      !is.na(med_ana_dat2$biochem_UA) &
                      !is.na(med_ana_dat2$ascvd) &
                      !is.na(med_ana_dat2$biochem_HDL), ]

res <- lapply(cols, batch_mediate,
              y = "ascvd", x = "abs_age_gap_adjust", df = tmp)
res <- ldply(res, data.frame)
# `res` is the source data frame for Figure 6b (mediation diagram).
