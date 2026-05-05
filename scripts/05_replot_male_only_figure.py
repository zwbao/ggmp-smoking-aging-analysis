#!/usr/bin/env python3
"""Re-render supp_male_only_sensitivity.{png,pdf} with the MaAsLin2-based
right panel.

Background
----------
The original right panel in `supp_male_only_sensitivity.png` plotted male-only
OTU smoking coefficients computed via Python OLS on arcsin-sqrt(rel_abundance)
against the full-sample MaAsLin2 LOG-transform coefficients from
`smk_status_sig_res.tsv`. The two coefficient sets live in different transform
spaces, so the y=x diagonal had no meaningful interpretation and points
appeared to lie far below it. `run_maaslin2_male_only_otu.R` re-fits the
male-only model in the published MaAsLin2 LOG framework so the two
coefficient sets are directly comparable; this script drops those new
coefficients into the figure.

Inputs
------
- `outputs/male_only_module_effects_by_age.tsv` (existing; left panel,
  unchanged)
- `outputs/male_only_overlap_otu_concordance_maaslin2.tsv` (new MaAsLin2
  male-only male coefs aligned with full-sample LOG coefs)

Outputs
-------
- `outputs/supp_male_only_sensitivity.png`
- `outputs/supp_male_only_sensitivity.pdf`

The buggy old PNG/PDF are renamed to *_arcsinsqrt_buggy.{png,pdf} for diagnostic
record, before being overwritten.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid", context="talk")


def main() -> None:
    # Run from the repository root: outputs/ holds intermediate TSVs.
    out_dir = Path("outputs")

    quartile_path = out_dir / "male_only_module_effects_by_age.tsv"
    maaslin_path = out_dir / "male_only_overlap_otu_concordance_maaslin2.tsv"
    fig_png = out_dir / "supp_male_only_sensitivity.png"
    fig_pdf = out_dir / "supp_male_only_sensitivity.pdf"

    # ------------------------------------------------------------------
    # Backup buggy figures (PNG + PDF) before we overwrite them.
    # ------------------------------------------------------------------
    backup_png = out_dir / "supp_male_only_sensitivity_arcsinsqrt_buggy.png"
    backup_pdf = out_dir / "supp_male_only_sensitivity_arcsinsqrt_buggy.pdf"
    if fig_png.exists() and not backup_png.exists():
        shutil.copy2(fig_png, backup_png)
        print(f"Backed up {fig_png.name} -> {backup_png.name}")
    if fig_pdf.exists() and not backup_pdf.exists():
        shutil.copy2(fig_pdf, backup_pdf)
        print(f"Backed up {fig_pdf.name} -> {backup_pdf.name}")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    quartile_df = pd.read_csv(quartile_path, sep="\t")
    # The quartile TSV stores age_categ ordered Quantile 1..4
    age_order = ["Quantile 1", "Quantile 2", "Quantile 3", "Quantile 4"]
    quartile_df["age_categ"] = pd.Categorical(
        quartile_df["age_categ"], categories=age_order, ordered=True
    )
    quartile_df = quartile_df.sort_values("age_categ").reset_index(drop=True)

    maaslin_df = pd.read_csv(maaslin_path, sep="\t")

    # ------------------------------------------------------------------
    # Build figure: 1x2, 14x6 — same layout as the buggy version
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel — quartile-level forest (unchanged)
    palette = ["#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"]
    axes[0].axvline(0, color="#9aa0a6", linestyle="--", linewidth=1)
    for i, row in quartile_df.iterrows():
        axes[0].errorbar(
            row["coef"],
            row["age_categ"],
            xerr=[[row["coef"] - row["ci_low"]], [row["ci_high"] - row["coef"]]],
            fmt="o",
            color=palette[i],
            markersize=10,
            capsize=4,
        )
    axes[0].set_xlabel("Adjusted effect of everyday smoking on pro-aging module score")
    axes[0].set_ylabel("Male age quartile")
    axes[0].set_title("Male-only sensitivity by age quartile")

    # Right panel — MaAsLin2 male-only vs full-sample MaAsLin2 LOG coefs
    scatter = maaslin_df.copy()
    same_dir = scatter["same_direction_as_full_maaslin2"].astype(bool)
    axes[1].scatter(
        scatter["smoking_coef"],
        scatter["male_smoking_coef_maaslin2"],
        c=np.where(same_dir, "#2a9d8f", "#d1495b"),
        alpha=0.8,
        edgecolor="white",
        linewidth=0.4,
    )
    lim = (
        np.nanmax(
            np.abs(
                scatter[["smoking_coef", "male_smoking_coef_maaslin2"]].to_numpy()
            )
        )
        * 1.1
    )
    axes[1].plot([-lim, lim], [-lim, lim], linestyle="--", color="#9aa0a6", linewidth=1)
    axes[1].axhline(0, color="#d0d7de", linewidth=1)
    axes[1].axvline(0, color="#d0d7de", linewidth=1)
    axes[1].set_xlim(-lim, lim)
    axes[1].set_ylim(-lim, lim)
    axes[1].set_xlabel("Full-sample MaAsLin2 smoking coefficient (LOG)")
    axes[1].set_ylabel("Male-only MaAsLin2 smoking coefficient (LOG)")
    axes[1].set_title("Concordance across concordant overlap OTUs")

    fig.tight_layout()
    fig.savefig(fig_png, dpi=300)
    fig.savefig(fig_pdf)
    plt.close(fig)
    print(f"Wrote: {fig_png}")
    print(f"Wrote: {fig_pdf}")

    # Console summary
    n = len(scatter)
    n_same = int(same_dir.sum())
    n_nominal = int((scatter["male_pval_maaslin2"] < 0.05).sum())
    n_fdr = int((scatter["male_qval_maaslin2"] < 0.05).sum())
    print(
        f"Right-panel summary: same-direction {n_same}/{n}, "
        f"nominal P<0.05 {n_nominal}/{n}, FDR q<0.05 {n_fdr}/{n}"
    )


if __name__ == "__main__":
    main()
