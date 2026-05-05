#!/usr/bin/env Rscript
# =============================================================================
# 02_maaslin2_reanalysis.R
#
# MaAsLin2 family/genus cross-check for the GGMP smoking x aging revision
# (Bao et al., NTR-2026-092). Confirms the family/genus higher-taxonomy
# robustness conclusions using the actual MaAsLin2 R package (the framework
# used in the main paper) and writes outputs alongside the Python OLS TSVs
# without overwriting them. Produces Supp Tables 18-19.
#
# Settings follow the original main-paper code (`code.R` in the original GGMP
# processing repository, https://github.com/SMUJYYXB/GGMP-Regional-variations) exactly:
#   normalization="TSS", transform="LOG", analysis_method="LM",
#   min_abundance=0.0, min_prevalence=0.0, correction="BH",
#   standardize=FALSE
# Covariates match the Python OLS sensitivity script (01_revision_analysis.py):
#   smoking model: smoke_binary + age + bmi + gender + bristol_cat + county;
#   age model (within never_smoker): age + bmi + gender + bristol_cat + county.
#
# Run from the repository root: `Rscript scripts/02_maaslin2_reanalysis.R`.
# Inputs are read from ./data/ and outputs/logs are written to
# ./outputs/ and ./logs/ relative to the current working directory.
# =============================================================================

suppressPackageStartupMessages({
  library(biomformat)
  library(Maaslin2)
  library(data.table)
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

dir.create(out_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(maaslin_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,      recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# Backfill Greengenes-style taxonomy: drop the f__/g__ prefix, fall back to
# Unclassified_<higher_rank> when the requested level is empty (matches the
# Python `taxonomy_label` rule).
taxonomy_label <- function(tax_vec, level) {
  level_idx <- c(domain = 1, phylum = 2, class = 3, order = 4,
                 family = 5, genus = 6, species = 7)[[level]]
  if (length(tax_vec) < 7) tax_vec <- c(tax_vec, rep("", 7 - length(tax_vec)))
  vals <- vapply(tax_vec, function(x) {
    if (is.na(x) || is.null(x)) return("")
    parts <- strsplit(x, "__", fixed = TRUE)[[1]]
    if (length(parts) >= 2) trimws(parts[2]) else trimws(x)
  }, character(1))
  cur <- vals[level_idx]
  if (!is.na(cur) && nzchar(cur)) return(cur)
  for (back in seq.int(level_idx - 1, 1)) {
    prev <- vals[back]
    if (!is.na(prev) && nzchar(prev)) return(paste0("Unclassified_", prev))
  }
  "Unclassified"
}

# ---------------------------------------------------------------------------
# Load metadata and cast columns
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
meta$smoke_binary <- ifelse(meta$smk_status == "everyday", "everyday",
                     ifelse(meta$smk_status == "never_smoker", "never_smoker", NA_character_))

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

# Pull taxonomy slot (list of vectors)
obs_md <- observation_metadata(biom)
tax_list <- if (is.data.frame(obs_md)) {
  setNames(lapply(seq_len(nrow(obs_md)), function(i) as.character(unlist(obs_md[i, ]))),
           rownames(obs_md))
} else {
  obs_md
}
# Sanity: align taxonomy list with feature ids
tax_list <- tax_list[otu_ids]

# Intersect samples to metadata subset
sample_keep <- intersect(meta$ID, sample_ids)
log_msg("Sample intersection (metadata vs BIOM): ", length(sample_keep))
meta <- meta[match(sample_keep, meta$ID), , drop = FALSE]
counts_mat <- counts_mat[, sample_keep, drop = FALSE]
rownames(meta) <- meta$ID

# ---------------------------------------------------------------------------
# Collapse counts to a taxonomic level
# ---------------------------------------------------------------------------
collapse_to_level <- function(counts, tax_list, level) {
  labels <- vapply(tax_list, taxonomy_label, level = level, FUN.VALUE = character(1))
  stopifnot(length(labels) == nrow(counts))
  agg <- rowsum(counts, group = labels, reorder = TRUE)
  agg
}

# ---------------------------------------------------------------------------
# Run a MaAsLin2 model
# ---------------------------------------------------------------------------
run_maaslin <- function(input_data, input_metadata, fixed_effects, reference,
                        out_subdir, label) {
  out <- file.path(maaslin_root, out_subdir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log_msg("Running MaAsLin2 [", label, "] on ", nrow(input_data),
          " features x ", ncol(input_data), " samples")
  res <- Maaslin2(
    input_data       = input_data,
    input_metadata   = input_metadata,
    output           = out,
    min_abundance    = 0.0,
    min_prevalence   = 0.0,
    normalization    = "TSS",
    transform        = "LOG",
    analysis_method  = "LM",
    max_significance = 1.0,           # keep all features in results table
    fixed_effects    = fixed_effects,
    reference        = reference,
    correction       = "BH",
    standardize      = FALSE,
    cores            = 4,
    plot_heatmap     = FALSE,
    plot_scatter     = FALSE
  )
  res$results
}

# ---------------------------------------------------------------------------
# Per-level pipeline
# ---------------------------------------------------------------------------
run_level <- function(level) {
  log_msg("=== Level: ", level, " ===")
  agg <- collapse_to_level(counts_mat, tax_list, level)
  log_msg("Collapsed to ", nrow(agg), " ", level, " bins")

  # ---- Smoking model (everyday vs never) ----
  smk_keep <- with(meta,
    !is.na(smoke_binary) & !is.na(age) & !is.na(bmi) &
      !is.na(bristol_cat) & bristol_cat != "NA" &
      !is.na(gender) & nzchar(gender) &
      !is.na(county) & nzchar(county))
  smk_meta <- meta[smk_keep, , drop = FALSE]
  smk_meta$smoke_binary <- factor(smk_meta$smoke_binary,
                                  levels = c("never_smoker", "everyday"))
  smk_meta$gender       <- factor(smk_meta$gender)
  smk_meta$bristol_cat  <- factor(smk_meta$bristol_cat)
  smk_meta$county       <- factor(smk_meta$county)
  rownames(smk_meta) <- smk_meta$ID
  smk_counts <- agg[, smk_meta$ID, drop = FALSE]
  smk_counts <- smk_counts[rowSums(smk_counts) > 0, , drop = FALSE]

  smk_input_data <- t(smk_counts)  # rows=samples for MaAsLin2
  smk_input_meta <- smk_meta[, c("smoke_binary", "age", "bmi", "gender",
                                 "bristol_cat", "county"), drop = FALSE]

  log_msg("Smoking subset: ", nrow(smk_input_data), " samples (",
          sum(smk_meta$smoke_binary == "everyday"), " everyday, ",
          sum(smk_meta$smoke_binary == "never_smoker"), " never_smoker)")

  smk_res <- run_maaslin(
    input_data     = smk_input_data,
    input_metadata = smk_input_meta,
    fixed_effects  = c("smoke_binary", "age", "bmi", "gender",
                       "bristol_cat", "county"),
    reference      = c("smoke_binary,never_smoker", "gender,f",
                       "bristol_cat,4", "county,G440282"),
    out_subdir     = paste0(level, "_smoking"),
    label          = paste0(level, " smoking")
  )
  smk_out <- subset(smk_res, metadata == "smoke_binary" & value == "everyday")
  smk_out <- smk_out[order(smk_out$qval, smk_out$pval), ]
  smk_out$n_samples <- nrow(smk_input_data)
  fwrite(smk_out,
         file = file.path(out_dir, paste0(level, "_smoking_results_maaslin2.tsv")),
         sep = "\t", na = "NA")

  # ---- Age model within never_smoker ----
  age_keep <- with(meta,
    smk_status == "never_smoker" &
      !is.na(age) & !is.na(bmi) &
      !is.na(bristol_cat) & bristol_cat != "NA" &
      !is.na(gender) & nzchar(gender) &
      !is.na(county) & nzchar(county))
  age_meta <- meta[age_keep, , drop = FALSE]
  age_meta$gender      <- factor(age_meta$gender)
  age_meta$bristol_cat <- factor(age_meta$bristol_cat)
  age_meta$county      <- factor(age_meta$county)
  rownames(age_meta) <- age_meta$ID
  age_counts <- agg[, age_meta$ID, drop = FALSE]
  age_counts <- age_counts[rowSums(age_counts) > 0, , drop = FALSE]

  age_input_data <- t(age_counts)
  age_input_meta <- age_meta[, c("age", "bmi", "gender",
                                 "bristol_cat", "county"), drop = FALSE]

  log_msg("Age subset (never smokers): ", nrow(age_input_data), " samples")

  age_res <- run_maaslin(
    input_data     = age_input_data,
    input_metadata = age_input_meta,
    fixed_effects  = c("age", "bmi", "gender", "bristol_cat", "county"),
    reference      = c("gender,f", "bristol_cat,4", "county,G440282"),
    out_subdir     = paste0(level, "_age"),
    label          = paste0(level, " age")
  )
  age_out <- subset(age_res, metadata == "age")
  age_out <- age_out[order(age_out$qval, age_out$pval), ]
  age_out$n_samples <- nrow(age_input_data)
  fwrite(age_out,
         file = file.path(out_dir, paste0(level, "_age_results_maaslin2.tsv")),
         sep = "\t", na = "NA")

  # ---- Shared (smoking ∩ age) ----
  shared <- merge(smk_out[, c("feature", "coef", "pval", "qval")],
                  age_out[, c("feature", "coef", "pval", "qval")],
                  by = "feature",
                  suffixes = c("_smoking", "_age"))
  shared$direction_concordant <- sign(shared$coef_smoking) == sign(shared$coef_age)
  shared$both_fdr_sig <- shared$qval_smoking < 0.05 & shared$qval_age < 0.05
  shared <- shared[order(shared$qval_smoking + shared$qval_age), ]
  fwrite(shared,
         file = file.path(out_dir, paste0(level, "_shared_results_maaslin2.tsv")),
         sep = "\t", na = "NA")

  list(
    level                  = level,
    smoking_tested         = nrow(smk_out),
    age_tested             = nrow(age_out),
    smoking_fdr_sig        = sum(smk_out$qval < 0.05, na.rm = TRUE),
    age_fdr_sig            = sum(age_out$qval < 0.05, na.rm = TRUE),
    shared_features        = nrow(shared),
    shared_both_fdr_sig    = sum(shared$both_fdr_sig, na.rm = TRUE),
    shared_concordant_sig  = sum(shared$both_fdr_sig & shared$direction_concordant,
                                 na.rm = TRUE),
    smoking_n_samples      = nrow(smk_input_data),
    age_n_samples          = nrow(age_input_data)
  )
}

summaries <- lapply(c("family", "genus"), run_level)
names(summaries) <- c("family", "genus")

summary_df <- data.frame(
  level                 = c("family", "genus"),
  smoking_tested        = sapply(summaries, `[[`, "smoking_tested"),
  age_tested            = sapply(summaries, `[[`, "age_tested"),
  smoking_fdr_sig       = sapply(summaries, `[[`, "smoking_fdr_sig"),
  age_fdr_sig           = sapply(summaries, `[[`, "age_fdr_sig"),
  shared_features       = sapply(summaries, `[[`, "shared_features"),
  shared_both_fdr_sig   = sapply(summaries, `[[`, "shared_both_fdr_sig"),
  shared_concordant_sig = sapply(summaries, `[[`, "shared_concordant_sig"),
  smoking_n_samples     = sapply(summaries, `[[`, "smoking_n_samples"),
  age_n_samples         = sapply(summaries, `[[`, "age_n_samples")
)
fwrite(summary_df,
       file = file.path(out_dir, "maaslin2_summary.tsv"),
       sep = "\t")
log_msg("MaAsLin2 reanalysis complete.")
print(summary_df)
