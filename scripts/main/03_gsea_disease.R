# ============================================================================
# 03_gsea_disease.R
#
# GSEA-style enrichment of smoking- and age-related OTU lists across chronic
# disease phenotypes. Produces the data underlying **Figure 5** of Bao et al.
# (NTR-2026-092).
#
# This script is migrated as-is from `220914/gsea_dis.R` and uses
# `clusterProfiler::GSEA` with `DESeq2`-derived OTU rankings under a healthy
# (`health == "y"`) vs disease-positive case-control design.
#
# ----------------------------------------------------------------------------
# Inputs (expected — NOT redistributed)
# ----------------------------------------------------------------------------
#
# This script expects the following objects already present in the R session:
#
#   - `otu.mets.meta`  : project-specific metadata data.frame (rows = samples,
#                        with columns `Bristol_stool_type`, `health`, `age`,
#                        `gender`, `county_level_code`, and disease flags
#                        `dis_atherosclerosis`, `dis_fatty_liver`, `dis_T2DM`,
#                        `dis_hepatic_calculus`, `dis_gout`, `MetS`).
#   - `otu.mets.count` : raw OTU count matrix (rows = samples, columns = OTU
#                        IDs). +1 pseudocount is added internally before DESeq2.
#   - `dis_gmt2`       : a TERM2GENE 2-column data.frame
#                        (term name, OTU ID) used as the GSEA gene-set
#                        definition. Built from the smoking- and age-related
#                        OTU lists produced by `01_main_maaslin2.R`.
#
# These objects are constructed in the corresponding author's working tree
# (see `220914/221008.rda`, `220914/221019.rda`, etc., ~1.5-2 GB each) and
# are NOT redistributed in this repository. To reproduce:
#   1. Run `scripts/main/01_main_maaslin2.R` to obtain the smoking-, second-
#      hand-, and age-related OTU lists.
#   2. Build `dis_gmt2` as a long data.frame with one row per (gene-set,
#      OTU) pair, e.g.:
#        dis_gmt2 <- bind_rows(
#          data.frame(term = "smk_status3_up",  gene = up_otus),
#          data.frame(term = "smk_status3_down", gene = down_otus),
#          data.frame(term = "age_up",          gene = age_up_otus),
#          ...
#        )
#   3. Build `otu.mets.meta` and `otu.mets.count` from the analysis-ready
#      phyloseq object (`non_mach_abx_GPf.rel`) plus the count BIOM table
#      (`GGMP7009_even10k.biom`).
#
# ----------------------------------------------------------------------------
# Outputs
# ----------------------------------------------------------------------------
#
# Six in-memory GSEA result objects:
#   mets.y3.dis_atherosclerosis, mets.y3.dis_fatty_liver, mets.y3.dis_T2DM,
#   mets.y3.dis_hepatic_calculus, mets.y3.dis_gout, mets.y3.MetS
# and six matching `gseaplot2` ggplot objects assembled into Figure 5.
#
# The script also writes `.pptx` files via `export::graph2ppt` to a `./pic/`
# directory (kept for fidelity with the original; not strictly necessary).
#
# ----------------------------------------------------------------------------
# How to run
# ----------------------------------------------------------------------------
#
#   # interactive (recommended)
#   load("/path/to/221108.rda")  # or run scripts/main/01_main_maaslin2.R first
#   source("scripts/main/03_gsea_disease.R")
#
# ----------------------------------------------------------------------------
# Reproducibility notes
# ----------------------------------------------------------------------------
#
# - GSEA permutation count is `nPermSimple = 10000`, `eps = 1e-50`,
#   `pvalueCutoff = 2` (i.e. report all gene sets, do not filter by p).
# - DESeq2 fits are deterministic given the inputs.
# - `clusterProfiler::GSEA` uses an internal RNG; for fully reproducible
#   p-values across runs set `set.seed(...)` before each `GSEA()` call.
# ============================================================================

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(export)

# Filter to samples with non-missing Bristol stool type and define the healthy
# reference subset.
otu.mets.meta2 <- otu.mets.meta[!is.na(otu.mets.meta$Bristol_stool_type), ]
otu.mets.meta2.h <- otu.mets.meta2[otu.mets.meta2$health == "y", ]

### dis_atherosclerosis -----------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$dis_atherosclerosis) &
                                     otu.mets.meta2$dis_atherosclerosis == "y", ]
otu.mets.meta.dis_atherosclerosis <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.dis_atherosclerosis <- otu.mets.count[row.names(otu.mets.meta.dis_atherosclerosis), ]
otu.mets.count.dis_atherosclerosis <- otu.mets.count.dis_atherosclerosis + 1

mets.format.dis_atherosclerosis <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.dis_atherosclerosis)),
  otu.mets.meta.dis_atherosclerosis,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.dis_atherosclerosis <- DESeq(mets.format.dis_atherosclerosis)
mets.result.dis_atherosclerosis <- results(mets.analyze.dis_atherosclerosis,
                                           contrast = c("health", "n", "y"))
mets.result.dis_atherosclerosis <- mets.result.dis_atherosclerosis[order(mets.result.dis_atherosclerosis$log2FoldChange), ]
mets.geneList.dis_atherosclerosis <- mets.result.dis_atherosclerosis[, 2]
names(mets.geneList.dis_atherosclerosis) <- as.character(row.names(mets.result.dis_atherosclerosis))
mets.geneList.dis_atherosclerosis <- sort(mets.geneList.dis_atherosclerosis, decreasing = TRUE)

mets.y3.dis_atherosclerosis <- GSEA(mets.geneList.dis_atherosclerosis,
                                    TERM2GENE = dis_gmt2,
                                    pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.dis_atherosclerosis <- gseaplot2(mets.y3.dis_atherosclerosis, 1:6,
                                               pvalue_table = TRUE, title = "Atherosclerosis"))

### dis_fatty_liver ---------------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$dis_fatty_liver) &
                                     otu.mets.meta2$dis_fatty_liver == "y", ]
otu.mets.meta.dis_fatty_liver <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.dis_fatty_liver <- otu.mets.count[row.names(otu.mets.meta.dis_fatty_liver), ]
otu.mets.count.dis_fatty_liver <- otu.mets.count.dis_fatty_liver + 1

mets.format.dis_fatty_liver <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.dis_fatty_liver)),
  otu.mets.meta.dis_fatty_liver,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.dis_fatty_liver <- DESeq(mets.format.dis_fatty_liver)
mets.result.dis_fatty_liver <- results(mets.analyze.dis_fatty_liver,
                                       contrast = c("health", "n", "y"))
mets.result.dis_fatty_liver <- mets.result.dis_fatty_liver[order(mets.result.dis_fatty_liver$log2FoldChange), ]
mets.geneList.dis_fatty_liver <- mets.result.dis_fatty_liver[, 2]
names(mets.geneList.dis_fatty_liver) <- as.character(row.names(mets.result.dis_fatty_liver))
mets.geneList.dis_fatty_liver <- sort(mets.geneList.dis_fatty_liver, decreasing = TRUE)

mets.y3.dis_fatty_liver <- GSEA(mets.geneList.dis_fatty_liver,
                                TERM2GENE = dis_gmt2,
                                pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.dis_fatty_liver <- gseaplot2(mets.y3.dis_fatty_liver, 1:6,
                                           pvalue_table = TRUE, title = "Fatty liver"))

### dis_T2DM ----------------------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$dis_T2DM) &
                                     otu.mets.meta2$dis_T2DM == "y", ]
otu.mets.meta.dis_T2DM <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.dis_T2DM <- otu.mets.count[row.names(otu.mets.meta.dis_T2DM), ]
otu.mets.count.dis_T2DM <- otu.mets.count.dis_T2DM + 1

mets.format.dis_T2DM <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.dis_T2DM)),
  otu.mets.meta.dis_T2DM,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.dis_T2DM <- DESeq(mets.format.dis_T2DM)
mets.result.dis_T2DM <- results(mets.analyze.dis_T2DM, contrast = c("health", "n", "y"))
mets.result.dis_T2DM <- mets.result.dis_T2DM[order(mets.result.dis_T2DM$log2FoldChange), ]
mets.geneList.dis_T2DM <- mets.result.dis_T2DM[, 2]
names(mets.geneList.dis_T2DM) <- as.character(row.names(mets.result.dis_T2DM))
mets.geneList.dis_T2DM <- sort(mets.geneList.dis_T2DM, decreasing = TRUE)

mets.y3.dis_T2DM <- GSEA(mets.geneList.dis_T2DM,
                         TERM2GENE = dis_gmt2,
                         pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.dis_T2DM <- gseaplot2(mets.y3.dis_T2DM, 1:6,
                                    pvalue_table = TRUE, title = "T2DM"))

### dis_hepatic_calculus ----------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$dis_hepatic_calculus) &
                                     otu.mets.meta2$dis_hepatic_calculus == "y", ]
otu.mets.meta.dis_hepatic_calculus <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.dis_hepatic_calculus <- otu.mets.count[row.names(otu.mets.meta.dis_hepatic_calculus), ]
otu.mets.count.dis_hepatic_calculus <- otu.mets.count.dis_hepatic_calculus + 1

mets.format.dis_hepatic_calculus <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.dis_hepatic_calculus)),
  otu.mets.meta.dis_hepatic_calculus,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.dis_hepatic_calculus <- DESeq(mets.format.dis_hepatic_calculus)
mets.result.dis_hepatic_calculus <- results(mets.analyze.dis_hepatic_calculus,
                                            contrast = c("health", "n", "y"))
mets.result.dis_hepatic_calculus <- mets.result.dis_hepatic_calculus[order(mets.result.dis_hepatic_calculus$log2FoldChange), ]
mets.geneList.dis_hepatic_calculus <- mets.result.dis_hepatic_calculus[, 2]
names(mets.geneList.dis_hepatic_calculus) <- as.character(row.names(mets.result.dis_hepatic_calculus))
mets.geneList.dis_hepatic_calculus <- sort(mets.geneList.dis_hepatic_calculus, decreasing = TRUE)

mets.y3.dis_hepatic_calculus <- GSEA(mets.geneList.dis_hepatic_calculus,
                                     TERM2GENE = dis_gmt2,
                                     pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.dis_hepatic_calculus <- gseaplot2(mets.y3.dis_hepatic_calculus, 1:6,
                                                pvalue_table = TRUE, title = "Hepatic calculus"))

### dis_gout ----------------------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$dis_gout) &
                                     otu.mets.meta2$dis_gout == "y", ]
otu.mets.meta.dis_gout <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.dis_gout <- otu.mets.count[row.names(otu.mets.meta.dis_gout), ]
otu.mets.count.dis_gout <- otu.mets.count.dis_gout + 1

mets.format.dis_gout <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.dis_gout)),
  otu.mets.meta.dis_gout,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.dis_gout <- DESeq(mets.format.dis_gout)
mets.result.dis_gout <- results(mets.analyze.dis_gout, contrast = c("health", "n", "y"))
mets.result.dis_gout <- mets.result.dis_gout[order(mets.result.dis_gout$log2FoldChange), ]
mets.geneList.dis_gout <- mets.result.dis_gout[, 2]
names(mets.geneList.dis_gout) <- as.character(row.names(mets.result.dis_gout))
mets.geneList.dis_gout <- sort(mets.geneList.dis_gout, decreasing = TRUE)

mets.y3.dis_gout <- GSEA(mets.geneList.dis_gout,
                         TERM2GENE = dis_gmt2,
                         pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.dis_gout <- gseaplot2(mets.y3.dis_gout, 1:6,
                                    pvalue_table = TRUE, title = "Gout"))

### MetS --------------------------------------------------------------------

otu.mets.meta2.n <- otu.mets.meta2[!is.na(otu.mets.meta2$MetS) &
                                     otu.mets.meta2$MetS == "y", ]
otu.mets.meta.MetS <- rbind(otu.mets.meta2.h, otu.mets.meta2.n)
otu.mets.count.MetS <- otu.mets.count[row.names(otu.mets.meta.MetS), ]
otu.mets.count.MetS <- otu.mets.count.MetS + 1

mets.format.MetS <- DESeqDataSetFromMatrix(
  data.frame(t(otu.mets.count.MetS)),
  otu.mets.meta.MetS,
  design = ~ health + age + gender + county_level_code + Bristol_stool_type
)
mets.analyze.MetS <- DESeq(mets.format.MetS)
mets.result.MetS <- results(mets.analyze.MetS, contrast = c("health", "n", "y"))
mets.result.MetS <- mets.result.MetS[order(mets.result.MetS$log2FoldChange), ]
mets.geneList.MetS <- mets.result.MetS[, 2]
names(mets.geneList.MetS) <- as.character(row.names(mets.result.MetS))
mets.geneList.MetS <- sort(mets.geneList.MetS, decreasing = TRUE)

mets.y3.MetS <- GSEA(mets.geneList.MetS,
                     TERM2GENE = dis_gmt2,
                     pvalueCutoff = 2, nPermSimple = 10000, eps = 1e-50)

(mets.gsea.p3.MetS <- gseaplot2(mets.y3.MetS, 1:6,
                                pvalue_table = TRUE, title = "MetS"))

# ---------------------------------------------------------------------------
# Re-render the panels selecting two top-ranked gene sets per disease for
# the published Figure 5.
# ---------------------------------------------------------------------------

(mets.gsea.p3.dis_atherosclerosis <- gseaplot2(mets.y3.dis_atherosclerosis,
                                               c(2, 3), subplots = 1:2,
                                               pvalue_table = FALSE,
                                               title = "Atherosclerosis"))
graph2ppt(mets.gsea.p3.dis_atherosclerosis,
          width = 5, height = 3, "./pic/mets.gsea.p3.dis_atherosclerosis.pptx")

(mets.gsea.p3.dis_fatty_liver <- gseaplot2(mets.y3.dis_fatty_liver,
                                           c(3, 2), subplots = 1:2,
                                           pvalue_table = FALSE,
                                           title = "Fatty liver"))
graph2ppt(mets.gsea.p3.dis_fatty_liver,
          width = 5, height = 3, "./pic/mets.gsea.p3.dis_fatty_liver.pptx")

(mets.gsea.p3.dis_T2DM <- gseaplot2(mets.y3.dis_T2DM,
                                    c(1, 3), subplots = 1:2,
                                    pvalue_table = FALSE, title = "T2DM"))
graph2ppt(mets.gsea.p3.dis_T2DM,
          width = 5, height = 3, "./pic/mets.gsea.p3.dis_T2DM.pptx")

(mets.gsea.p3.dis_hepatic_calculus <- gseaplot2(mets.y3.dis_hepatic_calculus,
                                                c(4, 3), subplots = 1:2,
                                                pvalue_table = FALSE,
                                                title = "Hepatic calculus"))
graph2ppt(mets.gsea.p3.dis_hepatic_calculus,
          width = 5, height = 3, "./pic/mets.gsea.p3.dis_hepatic_calculus.pptx")

(mets.gsea.p3.dis_gout <- gseaplot2(mets.y3.dis_gout,
                                    c(2, 3), subplots = 1:2,
                                    pvalue_table = FALSE, title = "Gout"))
graph2ppt(mets.gsea.p3.dis_gout,
          width = 5, height = 3, "./pic/mets.gsea.p3.dis_gout.pptx")

(mets.gsea.p3.MetS <- gseaplot2(mets.y3.MetS,
                                c(2, 4), subplots = 1:2,
                                pvalue_table = FALSE, title = "MetS"))
graph2ppt(mets.gsea.p3.MetS,
          width = 5, height = 3, "./pic/mets.gsea.p3.MetS.pptx")
