#!/usr/bin/env python3
"""Refit the male-only per-quartile pro-aging module regression with the
*correct* sample filter, fixing the np.where('nan') string-encoding bug in
`run_revision_analysis.py:clean_metadata` (line 60), which leaks 523
former_smoker / non-daily smoker males into the analysis as `smoke_everyday=0`.

We reproduce the module construction logic of
`run_revision_analysis.py:male_only_sensitivity` (lines 273-300) verbatim:
  - Pull the 40 direction-concordant overlap OTUs (sign(smoking_coef) ==
    sign(age_coef)).
  - Take their counts, divide by sample column-sums to get relative
    abundance, apply arcsin(sqrt(.)).
  - Orient each OTU by sign(age_coef).
  - z-score each OTU column on the male-subset, then average across the
    40 OTUs to get pro_aging_score.
  - Restrict to males with non-missing age, BMI, Bristol, age_categ.

Then we fit the per-quartile module-vs-smoking regression two ways:

    Approach A — exposure-as-binary, drop former+non-daily males:
        pro_aging_score ~ smoke_everyday + age + bmi + C(bristol_cat) + C(county)
        on males with smk_status in {everyday, never_smoker} only.

    Approach B — exposure-as-4-level-factor, all males kept:
        pro_aging_score ~ C(smk_status, Treatment(reference="never_smoker"))
                        + age + bmi + C(bristol_cat) + C(county)
        on all males. Extracts the smk_status[T.everyday] coefficient.
        This is the spec used in the published full-sample MaAsLin2 model
        and in `run_maaslin2_male_only_otu.R`.

Both approaches are written to separate TSVs. Approach B is also used to
re-render `outputs/supp_male_only_sensitivity.{png,pdf}` (overwriting the
previous fixed-right-panel-only render from `replot_male_only_figure.py`).
The pre-existing buggy backup `supp_male_only_sensitivity_arcsinsqrt_buggy.{png,pdf}`
is preserved.

Outputs
-------
- outputs/male_only_module_effects_by_age_corrected_everyday_vs_never.tsv  (Approach A)
- outputs/male_only_module_effects_by_age_corrected_4level.tsv             (Approach B)
- outputs/supp_male_only_sensitivity.{png,pdf}                             (re-rendered, Approach B left panel)
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import seaborn as sns
import statsmodels.api as sm
from biom import load_table
from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection


sns.set_theme(style="whitegrid", context="talk")


# ---------------------------------------------------------------------------
# I/O helpers (mirrors run_revision_analysis.py)
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--biom-path",
        default="data/GGMP7009_even10k.biom",
        help="Path to the public GGMP7009_even10k.biom file.",
    )
    parser.add_argument(
        "--metadata-path",
        default="data/GPf_metadata.tsv",
        help="Path to the project-specific GPf metadata table.",
    )
    parser.add_argument(
        "--overlap-path",
        default="outputs/otu_overlap_from_original_results.tsv",
        help="OTU overlap table produced by 01_revision_analysis.py.",
    )
    parser.add_argument(
        "--maaslin-male-otu-path",
        default="outputs/male_only_overlap_otu_concordance_maaslin2.tsv",
        help="Male-only MaAsLin2 OTU concordance produced by 02b_maaslin2_male_only_otu.R.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory for all outputs.",
    )
    return parser.parse_args()


def clean_metadata_corrected(path: Path) -> pd.DataFrame:
    """Same as run_revision_analysis.clean_metadata, but with the smoke_binary
    bug fixed (use pd.Series-based ifelse so the NaN third branch is real NaN,
    not the literal string 'nan')."""
    meta = pd.read_csv(path, sep="\t", low_memory=False).copy()
    meta["ID"] = meta["ID"].astype(str)
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["bmi"] = pd.to_numeric(meta["anthrop_BMI"], errors="coerce")
    meta["bristol"] = pd.to_numeric(meta["Bristol_stool_type"], errors="coerce")
    meta["bristol_cat"] = meta["bristol"].round().astype("Int64").astype(str)
    meta["gender"] = meta["gender"].astype(str)
    meta["county"] = meta["county_level_code"].astype(str)
    meta["age_categ"] = pd.Categorical(
        meta["age_categ"],
        categories=["Quantile 1", "Quantile 2", "Quantile 3", "Quantile 4"],
        ordered=True,
    )
    meta["smk_status"] = meta["smk_status"].astype(str)
    # CORRECT smoke_binary construction: pd.Series-based, so NaN stays NaN
    sb = pd.Series(pd.NA, index=meta.index, dtype="object")
    sb[meta["smk_status"].eq("everyday")] = "everyday"
    sb[meta["smk_status"].eq("never_smoker")] = "never_smoker"
    meta["smoke_binary"] = sb
    meta["smoke_everyday"] = (meta["smoke_binary"] == "everyday").astype(float)
    return meta


def subset_biom(biom_path: Path, sample_ids: Iterable[str]):
    table = load_table(str(biom_path))
    sample_ids = list(map(str, sample_ids))
    return table.filter(sample_ids, axis="sample", inplace=False)


def feature_table(table, feature_ids: list[str]) -> pd.DataFrame:
    sub = table.filter(feature_ids, axis="observation", inplace=False)
    data = sub.matrix_data.toarray()
    return pd.DataFrame(
        data,
        index=sub.ids(axis="observation"),
        columns=sub.ids(axis="sample"),
    )


def relative_abundance(counts: pd.DataFrame) -> pd.DataFrame:
    col_sums = counts.sum(axis=0)
    return counts.divide(col_sums, axis=1)


def arcsin_sqrt_transform(rel: pd.DataFrame) -> pd.DataFrame:
    return np.arcsin(np.sqrt(rel.clip(lower=0)))


# ---------------------------------------------------------------------------
# Module construction — copy-paste from run_revision_analysis.male_only_sensitivity
# ---------------------------------------------------------------------------

def build_male_module(
    table,
    meta: pd.DataFrame,
    overlap: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (male_df_with_score, concordant_overlap)."""
    # Use ALL males with non-missing covariates (smoke_binary may be NA — kept
    # for Approach B). This matches run_revision_analysis line 279-286 EXCEPT
    # we drop the smoke_binary.notna() filter because Approach B uses all 4
    # smk_status levels.
    male = meta[
        meta["gender"].eq("m")
        & meta["smk_status"].isin(
            ["everyday", "never_smoker", "former_smoker", "not_everyday"]
        )
        & meta["age"].notna()
        & meta["bmi"].notna()
        & meta["bristol"].notna()
        & meta["age_categ"].notna()
    ].copy()

    concordant = overlap[overlap["direction_concordant"]].copy()
    concordant["orientation"] = np.sign(concordant["age_coef"])
    feature_ids = concordant["feature"].astype(str).tolist()
    counts = feature_table(table, feature_ids)
    rel = relative_abundance(counts).reindex(columns=male["ID"])
    transformed = arcsin_sqrt_transform(rel).T  # samples x OTUs
    oriented = transformed.copy()
    orient_map = concordant.set_index("feature")["orientation"].to_dict()
    for feature, sign in orient_map.items():
        oriented[str(feature)] = oriented[str(feature)] * float(sign)
    score = oriented.apply(lambda col: (col - col.mean()) / col.std(ddof=0))
    male = male.reset_index(drop=True)
    male["pro_aging_score"] = score.mean(axis=1).values
    return male, concordant


# ---------------------------------------------------------------------------
# Per-quartile OLS for both approaches
# ---------------------------------------------------------------------------

def fit_one_quartile_approachA(sub: pd.DataFrame) -> dict:
    """Approach A: drop former + non-daily, smoke_everyday is binary."""
    sub_e = sub[sub["smk_status"].isin(["everyday", "never_smoker"])].copy()
    if sub_e.shape[0] < 10:
        return None
    # Drop categorical levels that have only one unique value (statsmodels chokes)
    formula = (
        "pro_aging_score ~ smoke_everyday + age + bmi"
        " + C(bristol_cat) + C(county)"
    )
    try:
        y, X = patsy.dmatrices(formula, data=sub_e, return_type="dataframe")
    except patsy.PatsyError as e:
        return {"_error": str(e)}
    res = sm.OLS(y, X).fit()
    coef_name = "smoke_everyday"
    coef = float(res.params[coef_name])
    se = float(res.bse[coef_name])
    pval = float(res.pvalues[coef_name])
    ci = res.conf_int(alpha=0.05).loc[coef_name]
    return {
        "coef": coef,
        "stderr": se,
        "pval": pval,
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
        "n": int(X.shape[0]),
        "n_everyday": int((sub_e["smk_status"] == "everyday").sum()),
        "n_never": int((sub_e["smk_status"] == "never_smoker").sum()),
    }


def fit_one_quartile_approachB(sub: pd.DataFrame) -> dict:
    """Approach B: 4-level smk_status with reference=never_smoker."""
    if sub.shape[0] < 10:
        return None
    formula = (
        'pro_aging_score ~ C(smk_status, Treatment(reference="never_smoker"))'
        " + age + bmi + C(bristol_cat) + C(county)"
    )
    try:
        y, X = patsy.dmatrices(formula, data=sub, return_type="dataframe")
    except patsy.PatsyError as e:
        return {"_error": str(e)}
    res = sm.OLS(y, X).fit()
    # The everyday coefficient column name in patsy:
    candidates = [
        c for c in X.columns
        if "smk_status" in c and "everyday" in c
    ]
    if not candidates:
        return {"_error": "no everyday coefficient column"}
    coef_name = candidates[0]
    coef = float(res.params[coef_name])
    se = float(res.bse[coef_name])
    pval = float(res.pvalues[coef_name])
    ci = res.conf_int(alpha=0.05).loc[coef_name]
    counts = sub["smk_status"].value_counts()
    return {
        "coef": coef,
        "stderr": se,
        "pval": pval,
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
        "coef_name": coef_name,
        "n": int(X.shape[0]),
        "n_everyday": int(counts.get("everyday", 0)),
        "n_never": int(counts.get("never_smoker", 0)),
        "n_former": int(counts.get("former_smoker", 0)),
        "n_notdaily": int(counts.get("not_everyday", 0)),
    }


def run_quartile_panel(
    male: pd.DataFrame, fit_fn, label: str
) -> pd.DataFrame:
    rows = []
    for quartile in male["age_categ"].cat.categories:
        sub = male[male["age_categ"] == quartile].copy()
        res = fit_fn(sub)
        if res is None:
            print(f"[{label}] Quartile {quartile}: skipped (n too small)")
            continue
        if "_error" in res:
            print(f"[{label}] Quartile {quartile}: ERROR {res['_error']}")
            continue
        res["age_categ"] = quartile
        rows.append(res)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    _, qvals = fdrcorrection(df["pval"].values)
    df["qval"] = qvals
    return df


# ---------------------------------------------------------------------------
# Figure rendering — both panels
# ---------------------------------------------------------------------------

def render_supp_figure(
    quartile_df: pd.DataFrame,
    maaslin_df: pd.DataFrame,
    out_dir: Path,
) -> None:
    fig_png = out_dir / "supp_male_only_sensitivity.png"
    fig_pdf = out_dir / "supp_male_only_sensitivity.pdf"
    backup_png = out_dir / "supp_male_only_sensitivity_arcsinsqrt_buggy.png"
    backup_pdf = out_dir / "supp_male_only_sensitivity_arcsinsqrt_buggy.pdf"
    # Don't overwrite the buggy backup if it exists
    if fig_png.exists() and not backup_png.exists():
        shutil.copy2(fig_png, backup_png)
        print(f"Backed up old figure to {backup_png.name}")
    if fig_pdf.exists() and not backup_pdf.exists():
        shutil.copy2(fig_pdf, backup_pdf)
        print(f"Backed up old figure to {backup_pdf.name}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Left panel: corrected quartile-level forest ----
    palette = ["#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"]
    qplot = quartile_df.copy()
    age_order = ["Quantile 1", "Quantile 2", "Quantile 3", "Quantile 4"]
    qplot["age_categ"] = pd.Categorical(
        qplot["age_categ"], categories=age_order, ordered=True
    )
    qplot = qplot.sort_values("age_categ").reset_index(drop=True)
    axes[0].axvline(0, color="#9aa0a6", linestyle="--", linewidth=1)
    for i, row in qplot.iterrows():
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

    # ---- Right panel: MaAsLin2 male-only vs full-sample LOG (already correct) ----
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
    print(f"Wrote {fig_png}")
    print(f"Wrote {fig_pdf}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading metadata from {args.metadata_path}")
    meta = clean_metadata_corrected(Path(args.metadata_path))
    n_pre = (meta["gender"].eq("m") & meta["smk_status"].notna()).sum()
    print(f"Total males with non-empty smk_status: {n_pre}")
    print(f"  smk_status breakdown (males):\n{meta.loc[meta['gender'].eq('m'), 'smk_status'].value_counts()}")

    print(f"Loading BIOM and overlap")
    overlap = pd.read_csv(args.overlap_path, sep="\t")
    overlap["feature"] = overlap["feature"].astype(str)
    table = subset_biom(Path(args.biom_path), meta["ID"])

    male, concordant = build_male_module(table, meta, overlap)
    print(f"\nMale subset (kept smk_status in {{everyday, never, former, notdaily}} + complete covariates): n = {male.shape[0]}")
    print(f"  smk_status breakdown:\n{male['smk_status'].value_counts()}")

    # =====================================================================
    # Approach A: smoke_everyday vs never_smoker, drop former + notdaily
    # =====================================================================
    print("\n=== Approach A: smoke_everyday (binary), drop former + not_everyday ===")
    male_A = male[male["smk_status"].isin(["everyday", "never_smoker"])].copy()
    print(f"Males kept for Approach A: n = {male_A.shape[0]}")
    print(f"  per-quartile counts:\n{male_A.groupby('age_categ', observed=True).size()}")

    df_A = run_quartile_panel(male_A, fit_one_quartile_approachA, "Approach A")
    df_A["approach"] = "A_everyday_vs_never"
    cols_A = [
        "age_categ", "approach", "coef", "stderr", "pval", "qval",
        "ci_low", "ci_high", "n", "n_everyday", "n_never",
    ]
    df_A = df_A[cols_A]
    out_A = out_dir / "male_only_module_effects_by_age_corrected_everyday_vs_never.tsv"
    df_A.to_csv(out_A, sep="\t", index=False)
    print(f"\nApproach A results -> {out_A}")
    print(df_A.to_string(index=False))

    # =====================================================================
    # Approach B: 4-level smk_status, reference=never_smoker
    # =====================================================================
    print("\n=== Approach B: C(smk_status, Treatment(ref=never_smoker)), all males ===")
    print(f"Males kept for Approach B: n = {male.shape[0]}")
    print(f"  per-quartile counts:\n{male.groupby('age_categ', observed=True).size()}")

    df_B = run_quartile_panel(male, fit_one_quartile_approachB, "Approach B")
    df_B["approach"] = "B_4level_factor"
    cols_B = [
        "age_categ", "approach", "coef", "stderr", "pval", "qval",
        "ci_low", "ci_high", "n", "n_everyday", "n_never",
        "n_former", "n_notdaily",
    ]
    df_B = df_B[cols_B]
    out_B = out_dir / "male_only_module_effects_by_age_corrected_4level.tsv"
    df_B.to_csv(out_B, sep="\t", index=False)
    print(f"\nApproach B results -> {out_B}")
    print(df_B.to_string(index=False))

    # =====================================================================
    # Re-render figure (Approach B left panel + MaAsLin2 right panel)
    # =====================================================================
    maaslin_path = Path(args.maaslin_male_otu_path)
    if not maaslin_path.exists():
        print(
            f"\nWARNING: {maaslin_path} not found; skipping figure re-render."
            " Run run_maaslin2_male_only_otu.R first."
        )
        return
    maaslin_df = pd.read_csv(maaslin_path, sep="\t")
    print("\n=== Re-rendering supp_male_only_sensitivity figure (Approach B left, MaAsLin2 right) ===")
    render_supp_figure(df_B, maaslin_df, out_dir)


if __name__ == "__main__":
    main()
