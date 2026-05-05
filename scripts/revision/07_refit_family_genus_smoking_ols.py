#!/usr/bin/env python3
"""Refit the family-level and genus-level OLS smoking models with the
*correct* sample filter, fixing the np.where('nan') string-encoding bug
that propagated through `01_revision_analysis.run_taxonomic_robustness`
(line 400 of that script).

Produces the manuscript-ready data for Supp Tables 12, 14, 15, 17 and
Supp Figs 8-9 (family- and genus-level smoking + shared scatter plots),
superseding the corresponding sub-pipelines of `01_revision_analysis.py`.

The bug:
    smoking_meta = meta[meta["smoke_binary"].notna() & ...]
relies on `smoke_binary`, which `clean_metadata` builds via
    np.where(eq("everyday"), "everyday",
             np.where(eq("never_smoker"), "never_smoker", np.nan))
The inner `np.where` returns an `object` array of strings, so the third
branch becomes the literal string "nan" rather than NaN. Then
`.notna()` is True for ALL rows, and former_smoker / not_everyday rows
slip through with `smoke_everyday=0`, contaminating the never-smoker
baseline.

Fix:
    smoking_meta = meta[
        meta["smk_status"].isin(["everyday", "never_smoker"])
        & meta["age"].notna() & meta["bmi"].notna() & meta["bristol"].notna()
    ]

This restricts to the intended n=5,926 (was 6,496 buggy).

This refit script:
  - reproduces the family- and genus-level OLS smoking model from
    `run_revision_analysis.run_taxonomic_robustness` lines 393-470 verbatim,
    with the corrected filter;
  - rejoins the (CLEAN) age-side columns from the existing
    `*_age_results.tsv` to produce a corrected `*_shared_results.tsv`;
  - re-renders the family / genus overlap scatter figure with the
    corrected smoking coefficients;
  - preserves the buggy outputs as `*_buggy_smoke_binary.{tsv,png,pdf}`.
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
from biom import load_table
from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection


sns.set_theme(style="whitegrid", context="talk")


# ---------------------------------------------------------------------------
# I/O helpers (copy-paste from run_revision_analysis.py)
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
    meta["smk_status"] = meta["smk_status"].astype(str)
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


def taxonomy_label(md: dict, level: str) -> str:
    level_map = {"family": 4, "genus": 5}
    idx = level_map[level]
    tax = list(md["taxonomy"])
    while len(tax) < 7:
        tax.append("")
    vals = [x.split("__", 1)[1] if "__" in x else x for x in tax]
    cur = vals[idx].strip()
    if cur:
        return cur
    for back in range(idx - 1, -1, -1):
        prev = vals[back].strip()
        if prev:
            return f"Unclassified_{prev}"
    return "Unclassified"


def collapse_table(table, level: str) -> pd.DataFrame:
    collapsed = table.collapse(
        lambda _id, md: taxonomy_label(md, level),
        axis="observation",
        norm=False,
    )
    data = collapsed.matrix_data.toarray()
    return pd.DataFrame(
        data,
        index=collapsed.ids(axis="observation"),
        columns=collapsed.ids(axis="sample"),
    )


def relative_abundance(counts: pd.DataFrame) -> pd.DataFrame:
    col_sums = counts.sum(axis=0)
    return counts.divide(col_sums, axis=1)


def arcsin_sqrt_transform(rel: pd.DataFrame) -> pd.DataFrame:
    return np.arcsin(np.sqrt(rel.clip(lower=0)))


def build_design(data: pd.DataFrame, formula: str) -> pd.DataFrame:
    return patsy.dmatrix(formula, data=data, return_type="dataframe")


def matrix_ols(y: pd.DataFrame, x: pd.DataFrame, coef_name: str) -> pd.DataFrame:
    x_mat = np.asarray(x, dtype=float)
    y_mat = np.asarray(y, dtype=float)
    xtx_inv = np.linalg.inv(x_mat.T @ x_mat)
    beta = xtx_inv @ x_mat.T @ y_mat
    resid = y_mat - x_mat @ beta
    n, p = x_mat.shape
    df_resid = n - p
    sse = np.sum(resid**2, axis=0)
    sigma2 = sse / df_resid
    coef_idx = list(x.columns).index(coef_name)
    se = np.sqrt(xtx_inv[coef_idx, coef_idx] * sigma2)
    coef = beta[coef_idx, :]
    t_stat = coef / se
    pval = 2 * t.sf(np.abs(t_stat), df=df_resid)
    ci_low = coef - t.ppf(0.975, df_resid) * se
    ci_high = coef + t.ppf(0.975, df_resid) * se
    return pd.DataFrame(
        {
            "feature": y.columns,
            "coef": coef,
            "stderr": se,
            "pval": pval,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "n": n,
        }
    )


def bh_qvalues(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    _, qvals = fdrcorrection(df["pval"].values)
    df["qval"] = qvals
    return df


def prevalence_filter(rel: pd.DataFrame, threshold: float = 0.1) -> pd.DataFrame:
    prev = (rel > 0).mean(axis=1)
    return rel.loc[prev >= threshold]


# ---------------------------------------------------------------------------
# The corrected refit (mirrors run_taxonomic_robustness, smoking-side only)
# ---------------------------------------------------------------------------

def refit_smoking_for_level(
    table,
    meta: pd.DataFrame,
    level: str,
    output_dir: Path,
    age_results_path: Path,
) -> dict:
    counts = collapse_table(table, level)
    ids = meta["ID"].tolist()
    counts = counts.loc[:, ids]
    rel = relative_abundance(counts)

    # CORRECTED filter: smk_status.isin(["everyday", "never_smoker"])
    # was: meta["smoke_binary"].notna()  (which let "nan" strings through)
    smoking_meta = meta[
        meta["smk_status"].isin(["everyday", "never_smoker"])
        & meta["age"].notna()
        & meta["bmi"].notna()
        & meta["bristol"].notna()
    ].copy()
    smoking_ids = smoking_meta["ID"].tolist()
    smoking_rel = prevalence_filter(rel.loc[:, smoking_ids]).T
    smoking_y = arcsin_sqrt_transform(smoking_rel.T).T
    smoking_x = build_design(
        smoking_meta,
        "1 + smoke_everyday + age + bmi + C(gender) + C(bristol_cat) + C(county)",
    )
    smoking_res = bh_qvalues(matrix_ols(smoking_y, smoking_x, "smoke_everyday"))
    smoking_res = smoking_res.sort_values(["qval", "pval"]).reset_index(drop=True)

    # Save corrected smoking results
    out_smoking = output_dir / f"{level}_smoking_results_corrected.tsv"
    smoking_res.to_csv(out_smoking, sep="\t", index=False)
    print(f"  wrote {out_smoking} (n features = {len(smoking_res)}, n samples = {smoking_meta.shape[0]})")

    # Rejoin age-side columns from the existing CLEAN age TSV
    age_res = pd.read_csv(age_results_path, sep="\t")
    shared = smoking_res.merge(age_res, on="feature", suffixes=("_smoking", "_age"))
    shared["direction_concordant"] = (
        np.sign(shared["coef_smoking"]) == np.sign(shared["coef_age"])
    )
    shared["both_fdr_sig"] = (shared["qval_smoking"] < 0.05) & (shared["qval_age"] < 0.05)
    out_shared = output_dir / f"{level}_shared_results_corrected.tsv"
    shared.to_csv(out_shared, sep="\t", index=False)
    print(f"  wrote {out_shared}")

    # Re-render the overlap scatter
    fig, ax = plt.subplots(figsize=(7, 6))
    plot_df = shared[shared["both_fdr_sig"]].copy()
    if not plot_df.empty:
        ax.scatter(
            plot_df["coef_smoking"],
            plot_df["coef_age"],
            c=np.where(plot_df["direction_concordant"], "#2a9d8f", "#d1495b"),
            alpha=0.85,
            edgecolor="white",
            linewidth=0.4,
        )
        top = plot_df.assign(
            rank=plot_df["qval_smoking"] + plot_df["qval_age"]
        ).sort_values("rank").head(10)
        for _, row in top.iterrows():
            ax.text(row["coef_smoking"], row["coef_age"], row["feature"], fontsize=9)
    ax.axhline(0, color="#d0d7de", linewidth=1)
    ax.axvline(0, color="#d0d7de", linewidth=1)
    ax.set_xlabel("Smoking coefficient (everyday vs never)")
    ax.set_ylabel("Age coefficient (within never smokers)")
    ax.set_title(f"{level.capitalize()}-level overlap")
    fig.tight_layout()

    # Backup originals before overwriting
    for ext in ("png", "pdf"):
        orig = output_dir / f"supp_{level}_overlap_scatter.{ext}"
        bak = output_dir / f"supp_{level}_overlap_scatter_buggy_smoke_binary.{ext}"
        if orig.exists() and not bak.exists():
            shutil.copy2(orig, bak)
            print(f"  backed up {orig.name} -> {bak.name}")

    fig.savefig(output_dir / f"supp_{level}_overlap_scatter.png", dpi=300)
    fig.savefig(output_dir / f"supp_{level}_overlap_scatter.pdf")
    plt.close(fig)
    print(f"  re-rendered supp_{level}_overlap_scatter.{{png,pdf}}")

    return {
        "level": level,
        "n_samples_smoking": int(smoking_meta.shape[0]),
        "smoking_tested_features": int(smoking_res.shape[0]),
        "smoking_fdr_sig": int((smoking_res["qval"] < 0.05).sum()),
        "shared_both_fdr_sig": int(shared["both_fdr_sig"].sum()),
        "shared_concordant": int(
            (shared["both_fdr_sig"] & shared["direction_concordant"]).sum()
        ),
        "smoking_results": smoking_res,
        "shared_results": shared,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)

    # First, back up the buggy TSV outputs (one-time)
    for level in ("family", "genus"):
        for suffix in ("smoking_results.tsv", "shared_results.tsv"):
            orig = output_dir / f"{level}_{suffix}"
            bak = output_dir / f"{level}_{suffix.replace('.tsv', '_buggy_smoke_binary.tsv')}"
            if orig.exists() and not bak.exists():
                shutil.copy2(orig, bak)
                print(f"backed up {orig.name} -> {bak.name}")

    meta = clean_metadata_corrected(Path(args.metadata_path))
    table = subset_biom(Path(args.biom_path), meta["ID"])

    summaries = []
    for level in ("family", "genus"):
        print(f"--- {level} ---")
        age_path = output_dir / f"{level}_age_results.tsv"
        summary = refit_smoking_for_level(
            table, meta, level, output_dir, age_path
        )
        summaries.append(summary)

    print("\n=== Summary ===")
    for s in summaries:
        print(
            f"{s['level']}: n_samples_smoking={s['n_samples_smoking']}, "
            f"features={s['smoking_tested_features']}, "
            f"smoking_fdr_sig={s['smoking_fdr_sig']}, "
            f"shared_both_fdr_sig={s['shared_both_fdr_sig']}, "
            f"shared_concordant={s['shared_concordant']}"
        )


if __name__ == "__main__":
    main()
