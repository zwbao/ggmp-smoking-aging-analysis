#!/usr/bin/env Rscript
# MaAsLin2 male-only OTU-level smoking re-analysis for the 40 direction-concordant
# overlap OTUs flagged in male_only_overlap_otu_concordance.tsv.
#
# Background: the existing male_only_overlap_otu_concordance.tsv right-panel
# coefficients were computed via Python OLS on arcsin-sqrt-transformed relative
# abundance — this is in a different transform space than the published
# full-sample MaAsLin2 LOG coefficients in smk_status_sig_res.tsv. Numerical
# magnitudes therefore differ by ~17x, making the male-only coefficients look
# spuriously attenuated when plotted against the full-sample coefficients.
# This script re-runs the male-only OTU-level smoking model in the *same*
# MaAsLin2 LOG framework (TSS + LOG + LM, BH) so the two sets of coefficients
# live in the same statistical space and the y=x diagonal becomes meaningful.
#
# Settings replicate 220914/code.R: min_prevalence=0, min_abundance=0,
# normalization="TSS", transform="LOG", analysis_method="LM", correction="BH",
# standardize=FALSE.
#
# Inputs:
#   * GGMP7009_even10k.biom      - public OTU table (rarefied to 10k reads)
#   * GPf_metadata.tsv           - project metadata
#   * outputs/male_only_overlap_otu_concordance.tsv  - 40 concordant overlap OTUs
#                                                       (we copy smoking_coef,
#                                                       age_coef and direction
#                                                       columns through)
# Outputs:
#   * outputs/male_only_overlap_otu_concordance_maaslin2.tsv
#   * outputs/male_only_concordance_maaslin2_summary.json
#   * outputs/maaslin2_runs/male_only_otu_smoking/    (Maaslin2 working dir)

suppressPackageStartupMessages({
  library(biomformat)
  library(Maaslin2)
  library(data.table)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Paths (relative to the repository root / current working directory)
# ---------------------------------------------------------------------------
data_dir     <- "data"
biom_path    <- file.path(data_dir, "GGMP7009_even10k.biom")
meta_path    <- file.path(data_dir, "GPf_metadata.tsv")
out_dir      <- "outputs"
maaslin_root <- file.path(out_dir, "maaslin2_runs")
log_dir      <- "logs"
overlap_path <- file.path(out_dir, "male_only_overlap_otu_concordance.tsv")

dir.create(out_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(maaslin_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,      recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---------------------------------------------------------------------------
# Load metadata and cast columns. Match the casts in run_maaslin2_reanalysis.R
# ---------------------------------------------------------------------------
log_msg("Loading metadata: ", meta_path)
meta <- fread(meta_path, sep = "\t", data.table = FALSE,
              colClasses = "character", check.names = FALSE)
log_msg("Metadata: ", nrow(meta), " rows x ", ncol(meta), " cols")

meta$ID           <- as.character(meta$ID)
meta$age          <- suppressWarnings(as.numeric(meta$age))
meta$bmi          <- suppressWarnings(as.numeric(meta$anthrop_BMI))
meta$bristol_num  <- suppressWarnings(as.numeric(meta$Bristol_stool_type))
meta$bristol_cat  <- as.character(round(meta$bristol_num))
meta$gender       <- as.character(meta$gender)
meta$county       <- as.character(meta$county_level_code)
meta$smk_status   <- as.character(meta$smk_status)
meta$age_categ    <- as.character(meta$age_categ)

# ---------------------------------------------------------------------------
# Load BIOM table
# ---------------------------------------------------------------------------
log_msg("Loading BIOM: ", biom_path)
biom <- read_biom(biom_path)
counts_mat <- as.matrix(biom_data(biom))
storage.mode(counts_mat) <- "double"
otu_ids    <- rownames(counts_mat)
sample_ids <- colnames(counts_mat)
log_msg("BIOM: ", nrow(counts_mat), " features x ", ncol(counts_mat), " samples")

# Intersect samples to metadata subset
sample_keep <- intersect(meta$ID, sample_ids)
log_msg("Sample intersection (metadata vs BIOM): ", length(sample_keep))
meta <- meta[match(sample_keep, meta$ID), , drop = FALSE]
counts_mat <- counts_mat[, sample_keep, drop = FALSE]
rownames(meta) <- meta$ID

# ---------------------------------------------------------------------------
# Subset to males with complete cases on the smoking model covariates.
# Matches the analytic frame of the buggy run_revision_analysis.py
# male_only_sensitivity (which is also the analytic frame in the buggy TSV's
# `n` column = 2801): gender == "m" & smk_status non-missing & non-missing
# age, BMI, Bristol, age_categ. The published full-sample MaAsLin2 model
# treats smk_status as a 4-level factor (everyday / not_everyday /
# former_smoker / never_smoker) with reference = never_smoker; we follow
# that exactly so the male-only `smk_status everyday` coefficient is in the
# same statistical space as the full-sample `smk_status everyday`
# coefficient in smk_status_sig_res.tsv. Should reproduce n = 2,801.
# ---------------------------------------------------------------------------
male_keep <- with(meta,
  gender == "m" &
    !is.na(smk_status) & nzchar(smk_status) & smk_status != "NA" &
    !is.na(age) &
    !is.na(bmi) &
    !is.na(bristol_cat) & bristol_cat != "NA" &
    !is.na(age_categ) & nzchar(age_categ) & age_categ != "NA")
male_meta <- meta[male_keep, , drop = FALSE]
log_msg("Male complete-case subset: n = ", nrow(male_meta))
log_msg("smk_status breakdown: ",
        paste(names(table(male_meta$smk_status)),
              table(male_meta$smk_status), sep = "=", collapse = ", "))

# 4-level factor matching the published spec (reference = never_smoker)
male_meta$smk_status   <- factor(male_meta$smk_status,
                                 levels = c("never_smoker", "everyday",
                                            "not_everyday", "former_smoker"))
male_meta$bristol_cat  <- factor(male_meta$bristol_cat)
male_meta$county       <- factor(male_meta$county)
rownames(male_meta) <- male_meta$ID

# ---------------------------------------------------------------------------
# Load the 40 direction-concordant overlap OTUs and keep only their counts.
# We do NOT modify the buggy concordance TSV — we read its `feature` column
# and copy the `smoking_coef`, `age_coef`, `direction_concordant` columns
# (which are correct and live in the published MaAsLin2 LOG space) through
# to the new output.
# ---------------------------------------------------------------------------
log_msg("Loading 40 concordant overlap OTUs from: ", overlap_path)
overlap <- fread(overlap_path, sep = "\t", data.table = FALSE,
                 colClasses = "character", check.names = FALSE)
overlap$feature       <- as.character(overlap$feature)
overlap$smoking_coef  <- suppressWarnings(as.numeric(overlap$smoking_coef))
overlap$age_coef      <- suppressWarnings(as.numeric(overlap$age_coef))
overlap$direction_concordant <- as.logical(overlap$direction_concordant)
log_msg("Concordant overlap OTU rows: ", nrow(overlap))

feature_ids <- overlap$feature
feature_ids_in_biom <- intersect(feature_ids, otu_ids)
log_msg("Of ", length(feature_ids), " overlap OTUs, ",
        length(feature_ids_in_biom), " present in BIOM.")
stopifnot(length(feature_ids_in_biom) == length(feature_ids))

# Subset OTU table to males only and to the 40 concordant OTUs
male_counts <- counts_mat[feature_ids_in_biom, male_meta$ID, drop = FALSE]
log_msg("Male x 40-OTU count matrix: ", nrow(male_counts), " features x ",
        ncol(male_counts), " samples; row-sum range = ",
        min(rowSums(male_counts)), "..", max(rowSums(male_counts)))

# MaAsLin2 input shape: rows = samples, cols = features
# Note: MaAsLin2 / R's data.frame() applies make.names() to feature ids that
# start with a digit, prefixing them with "X" (e.g. 1044419 -> X1044419). We
# capture the original-ID -> safe-ID mapping so we can map MaAsLin2 results
# back to the original BIOM/overlap feature ids.
safe_feature_ids <- make.names(feature_ids_in_biom)
id_map <- data.frame(orig_feature = feature_ids_in_biom,
                     safe_feature = safe_feature_ids,
                     stringsAsFactors = FALSE)

input_data <- as.data.frame(t(male_counts))
colnames(input_data) <- safe_feature_ids
input_meta <- male_meta[, c("smk_status", "age", "bmi",
                            "bristol_cat", "county"), drop = FALSE]
# Rename to match the spec'd fixed_effects column names exactly:
#   smk_status, age, BMI, Bristol_stool_type, county_level_code
colnames(input_meta) <- c("smk_status", "age", "BMI",
                          "Bristol_stool_type", "county_level_code")

# ---------------------------------------------------------------------------
# Run MaAsLin2 with original-spec settings (TSS + LOG + LM, BH, no filters)
# ---------------------------------------------------------------------------
maaslin_out <- file.path(maaslin_root, "male_only_otu_smoking")
dir.create(maaslin_out, recursive = TRUE, showWarnings = FALSE)

log_msg("Running MaAsLin2 [male-only, 40 OTUs] on ",
        nrow(input_data), " samples x ", ncol(input_data), " features")
fit <- Maaslin2(
  input_data       = input_data,
  input_metadata   = input_meta,
  output           = maaslin_out,
  min_abundance    = 0.0,
  min_prevalence   = 0.0,
  normalization    = "TSS",
  transform        = "LOG",
  analysis_method  = "LM",
  max_significance = 1.0,           # keep all features in results table
  fixed_effects    = c("smk_status", "age", "BMI",
                       "Bristol_stool_type", "county_level_code"),
  reference        = c("smk_status,never_smoker",
                       "Bristol_stool_type,4",
                       "county_level_code,G440282"),
  correction       = "BH",
  standardize      = FALSE,
  cores            = 4,
  plot_heatmap     = FALSE,
  plot_scatter     = FALSE
)
res <- fit$results

# Pull the smk_status everyday coefficient + p-value for each OTU.
smk_rows <- subset(res, metadata == "smk_status" & value == "everyday")
log_msg("MaAsLin2 returned smk_status=everyday rows for ",
        nrow(smk_rows), " of ", length(feature_ids), " OTUs")

# Map safe_feature back to the original BIOM feature id
smk_rows$safe_feature <- as.character(smk_rows$feature)
smk_rows <- merge(smk_rows, id_map, by = "safe_feature", all.x = TRUE)

# Compute 95% CIs from coef/stderr (MaAsLin2 doesn't return them by default).
smk_rows$ci_low  <- smk_rows$coef - 1.96 * smk_rows$stderr
smk_rows$ci_high <- smk_rows$coef + 1.96 * smk_rows$stderr

# Keep one row per feature, using the original (BIOM) feature id
smk_rows <- smk_rows[, c("orig_feature", "coef", "stderr", "pval",
                         "ci_low", "ci_high")]
colnames(smk_rows) <- c("feature",
                        "male_smoking_coef_maaslin2",
                        "male_stderr_maaslin2",
                        "male_pval_maaslin2",
                        "male_ci_low_maaslin2",
                        "male_ci_high_maaslin2")

# ---------------------------------------------------------------------------
# Apply BH-FDR across the 40 OTUs only (apples-to-apples with the full-sample
# 40-OTU comparison; we are not re-using the qval that MaAsLin2 computes
# internally because that is computed across all features it tested, here it
# is also 40 so it should match — but we recompute explicitly to be safe).
# ---------------------------------------------------------------------------
smk_rows$male_qval_maaslin2 <- p.adjust(smk_rows$male_pval_maaslin2,
                                        method = "BH")

# ---------------------------------------------------------------------------
# Merge with the input concordance frame so we keep `smoking_coef`,
# `age_coef`, `direction_concordant`, and add `same_direction_as_full_maaslin2`.
# ---------------------------------------------------------------------------
keep_cols <- c("feature", "smoking_coef", "age_coef", "direction_concordant")
overlap_keep <- overlap[, keep_cols, drop = FALSE]

out_df <- merge(overlap_keep, smk_rows, by = "feature", all.x = TRUE,
                sort = FALSE)
out_df$same_direction_as_full_maaslin2 <- sign(out_df$male_smoking_coef_maaslin2) ==
                                          sign(out_df$smoking_coef)

# Reorder columns per spec
out_df <- out_df[, c("feature", "smoking_coef", "age_coef",
                     "direction_concordant",
                     "male_smoking_coef_maaslin2", "male_stderr_maaslin2",
                     "male_pval_maaslin2", "male_qval_maaslin2",
                     "male_ci_low_maaslin2", "male_ci_high_maaslin2",
                     "same_direction_as_full_maaslin2")]

out_path <- file.path(out_dir, "male_only_overlap_otu_concordance_maaslin2.tsv")
fwrite(out_df, file = out_path, sep = "\t", na = "NA")
log_msg("Wrote: ", out_path)

# ---------------------------------------------------------------------------
# Summary stats
# ---------------------------------------------------------------------------
n_total       <- nrow(out_df)
n_same_dir    <- sum(out_df$same_direction_as_full_maaslin2, na.rm = TRUE)
n_nominal_sig <- sum(out_df$male_pval_maaslin2 < 0.05, na.rm = TRUE)
n_fdr_sig     <- sum(out_df$male_qval_maaslin2 < 0.05, na.rm = TRUE)

# Spearman correlation between the buggy male coef (in old TSV) and the new
# MaAsLin2 male coef — should be high if the bug was purely a transform issue.
buggy <- fread(overlap_path, sep = "\t", data.table = FALSE,
               colClasses = "character", check.names = FALSE)
buggy$feature           <- as.character(buggy$feature)
buggy$male_smoking_coef <- suppressWarnings(as.numeric(buggy$male_smoking_coef))
cor_df <- merge(buggy[, c("feature", "male_smoking_coef")],
                out_df[, c("feature", "male_smoking_coef_maaslin2")],
                by = "feature")
spearman_old_vs_new <- suppressWarnings(
  cor(cor_df$male_smoking_coef, cor_df$male_smoking_coef_maaslin2,
      method = "spearman", use = "pairwise.complete.obs")
)
pearson_old_vs_new <- suppressWarnings(
  cor(cor_df$male_smoking_coef, cor_df$male_smoking_coef_maaslin2,
      method = "pearson", use = "pairwise.complete.obs")
)

summary_obj <- list(
  n_otus_total                  = as.integer(n_total),
  n_same_direction_maaslin2     = as.integer(n_same_dir),
  n_nominal_p_lt_0_05_maaslin2  = as.integer(n_nominal_sig),
  n_fdr_q_lt_0_05_maaslin2      = as.integer(n_fdr_sig),
  spearman_old_vs_new_male_coef = spearman_old_vs_new,
  pearson_old_vs_new_male_coef  = pearson_old_vs_new,
  male_n_samples                = as.integer(nrow(input_data)),
  male_everyday_n               = as.integer(sum(input_meta$smk_status == "everyday")),
  male_never_n                  = as.integer(sum(input_meta$smk_status == "never_smoker")),
  male_former_n                 = as.integer(sum(input_meta$smk_status == "former_smoker")),
  male_not_everyday_n           = as.integer(sum(input_meta$smk_status == "not_everyday")),
  buggy_n_same_direction        = 36L,
  buggy_n_nominal_p_lt_0_05     = 17L,
  buggy_n_fdr_q_lt_0_05         = 7L
)

summary_path <- file.path(out_dir, "male_only_concordance_maaslin2_summary.json")
write_json(summary_obj, path = summary_path, pretty = TRUE,
           auto_unbox = TRUE, digits = 8)
log_msg("Wrote: ", summary_path)

log_msg("=== Summary ===")
log_msg("n total OTUs                  : ", n_total)
log_msg("n same direction (MaAsLin2)   : ", n_same_dir, " (was 36 in buggy run)")
log_msg("n nominal P<0.05 (MaAsLin2)   : ", n_nominal_sig, " (was 17 in buggy run)")
log_msg("n FDR q<0.05 (MaAsLin2)       : ", n_fdr_sig, " (was 7 in buggy run)")
log_msg("Spearman(old buggy, new)      : ", round(spearman_old_vs_new, 4))
log_msg("Pearson(old buggy, new)       : ", round(pearson_old_vs_new, 4))
log_msg("Male-only OTU MaAsLin2 re-run done.")
