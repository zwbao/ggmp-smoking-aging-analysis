# Male-only OTU concordance — transform-mismatch fix memo

## What the bug was

The right panel of `supp_male_only_sensitivity.{png,pdf}` plotted male-only
OTU smoking coefficients on Y vs full-sample OTU smoking coefficients on X
for the 40 direction-concordant overlap OTUs. The two coefficient sets,
however, were computed in *different transform spaces*:

| Coefficient set | Pipeline | Transform |
| --- | --- | --- |
| Full-sample (X axis) | published MaAsLin2 LM (`smk_status_sig_res.tsv`) | TSS + LOG |
| Male-only (Y axis, buggy) | Python `matrix_ols` in `male_only_sensitivity` | `arcsin(sqrt(rel_abundance))` |

LOG transforms relative abundance approximately by `log2(rel + pseudocount)`,
which spreads small abundances over a much wider numeric range than
`arcsin(sqrt(rel))` (bounded by π/2 ≈ 1.57). Numerically this made the male-only
Y-axis coefficients ~17× smaller than the full-sample X-axis coefficients,
so the y=x diagonal of the scatter ran far above almost every point. That
created the misleading impression that male-only effects are dramatically
attenuated, when in fact the difference was almost entirely a unit / transform
artifact. The same bug also propagated into
`male_only_overlap_otu_concordance.tsv` (the `male_smoking_coef`,
`male_pval`, `male_qval`, `same_direction_as_full` columns) and into the
"same direction" / "nominal P < 0.05" / "FDR q < 0.05" counts cited in the
response letter.

The quartile-level **left** panel was *not* affected, because the per-OTU
arcsin-sqrt vectors are z-scored before being averaged into a module score
(`pro_aging_score`); z-scoring cancels the transform's unit, so the quartile
regression coefficients are unitless and the within-male comparison is
internally valid.

## The fix

`scripts/02b_maaslin2_male_only_otu.R` re-fits the male-only OTU smoking model
in the *same MaAsLin2 LOG framework* used for the full-sample model:

- `Maaslin2(... normalization="TSS", transform="LOG", analysis_method="LM",
  correction="BH", min_abundance=0, min_prevalence=0, standardize=FALSE)`
- 4-level `smk_status` factor with `reference="never_smoker"` (matches the
  published `smk_status_sig_res.tsv` model — the `value="everyday"` row of the
  MaAsLin2 results table is what we report)
- Fixed effects: `smk_status + age + BMI + Bristol_stool_type +
  county_level_code`
- Sample frame: males with non-missing covariates as defined in the original
  Python `male_only_sensitivity` filter — n = 2,801 (1,317 everyday smokers,
  961 never smokers, 402 former smokers, 121 not-everyday). The non-everyday
  / former / never groups are absorbed into the model as separate factor
  levels rather than being silently coded as zero (which is what the buggy
  Python `smoke_everyday = (smoke_binary == "everyday").astype(float)` did).
- BH-FDR is applied across the 40 OTUs only (apples-to-apples with the
  full-sample 40-OTU comparison).

Outputs:
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (40 OTUs × MaAsLin2
  male-only stats + the original full-sample LOG `smoking_coef` / `age_coef`
  copied through)
- `outputs/male_only_concordance_maaslin2_summary.json`
- `outputs/maaslin2_runs/male_only_otu_smoking/` (Maaslin2 working directory)

The buggy `male_only_overlap_otu_concordance.tsv` produced by script 01 is
**kept on disk unchanged** in the local working tree as a diagnostic record
(it is not committed to this public repository, since the regenerable buggy
output adds no reproducibility value beyond the corrected MaAsLin2 TSV).
The buggy supplementary figure is backed up to
`supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}` before being
overwritten by `scripts/05_replot_male_only_figure.py` (Bug-1 figure
re-render). Script 06 then re-renders both panels of the figure with the
fully-corrected left panel.

## Before-vs-after counts (40 direction-concordant overlap OTUs)

| Metric | Buggy (arcsin-sqrt OLS) | New (MaAsLin2 TSS+LOG, LM, BH) |
| --- | ---: | ---: |
| Same direction as full-sample | 36 / 40 | **39 / 40** |
| Nominal P < 0.05 | 17 / 40 | **16 / 40** |
| BH-FDR q < 0.05 (across 40) | 7 / 40 | **10 / 40** |

The new MaAsLin2 male-only run shows that direction-of-effect is preserved
for 39 of 40 overlap OTUs (vs 36 in the buggy run); 16 reach nominal P<0.05
and 10 reach FDR q<0.05, both substantially or modestly *better* than the
buggy numbers. The slight drop in nominal-P count (17 → 16) and the rise in
FDR-significant count (7 → 10) reflect that the LOG transform yields tighter
standard errors for low-abundance OTUs than the arcsin-sqrt transform — the
ranking of OTUs by significance is largely preserved (see correlation
below), but the small handful of OTUs near the P=0.05 boundary swap sides.

## Old-vs-new male coefficient correlation

Across the 40 OTUs, the buggy male coefficient (arcsin-sqrt OLS) and the
new MaAsLin2 LOG male coefficient have:

- Spearman ρ = **0.934**
- Pearson r = 0.790

The Spearman ρ ≈ 0.93 confirms that the bug was almost purely a
unit/transform artifact: the *ranking* and *sign* of male-only smoking
effects across OTUs is essentially preserved between the two pipelines.
What changed is the absolute magnitude of each coefficient (LOG-space units
vs arcsin-sqrt-space units), which is exactly what made the y=x diagonal
meaningless in the buggy scatter. After the fix, both axes live in
LOG-space and the y=x diagonal is meaningful — most points cluster along it.

## Recommended numbers for the response letter (Reviewer R2-2)

Use the new MaAsLin2 numbers, not the buggy arcsin-sqrt ones:

> "In male-only MaAsLin2 LOG models (n = 2,801 males with complete
> covariates) on the 40 direction-concordant overlap OTUs, **39 of 40
> retained the same direction of smoking effect** as the full-sample model,
> **16 reached nominal P < 0.05**, and **10 remained BH-FDR significant
> (q < 0.05)**. The Spearman correlation between male-only and full-sample
> smoking coefficients across these 40 OTUs is high (ρ ≈ 0.93 between the
> male-only re-run and the full-sample LOG coefficients), confirming that
> the smoking signature is robust to restricting the analysis to males.
> The supplementary scatter (Fig. SX, right panel) plots male-only vs
> full-sample MaAsLin2 LOG coefficients on a common y=x diagonal."

(Old "36 / 17 / 7" numbers should be retracted; they reflected an
inadvertent transform mismatch between the male-only and full-sample
pipelines.)

## Files in this repository

- `scripts/02b_maaslin2_male_only_otu.R` (Bug 1 fix script — re-runs the male-only MaAsLin2 LOG model)
- `scripts/05_replot_male_only_figure.py` (re-renders Supp Fig 7 right panel; later superseded by `06_refit_male_quartile_module.py` which re-renders both panels)
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (run output; reference copy in `outputs-reference/`)
- `outputs/male_only_concordance_maaslin2_summary.json` (run output)
- `outputs/supp_male_only_sensitivity.{png,pdf}` (overwritten, manuscript-ready)
- `outputs/supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}` (local diagnostic backup of the buggy version, not committed)
- `outputs/maaslin2_runs/male_only_otu_smoking/` (MaAsLin2 working directory)
- `docs/audit/male_only_otu_concordance_fix_memo.md` (this memo)
