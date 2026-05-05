#!/usr/bin/env python3
"""
Revision-analysis driver for the GGMP smoking x aging paper (Bao et al.,
NTR-2026-092). Produces three supplementary deliverables: (1) a sequencing-
depth / rarefaction supplement (Supp Fig 6 + Supp Table 9), (2) a male-only
sensitivity analysis on the pro-aging module (Supp Fig 7 + Supp Tables 10-11),
and (3) family- and genus-level OLS robustness re-tests (Supp Figs 8-9 + Supp
Tables 12-17). All inputs default to relative paths under ./data/ and outputs
are written to ./outputs/.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patsy
import seaborn as sns
from biom import load_table
from scipy.stats import t
from skbio.stats import subsample_counts
from statsmodels.stats.multitest import fdrcorrection


sns.set_theme(style="whitegrid", context="talk")


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
        "--smoking-results",
        default="data/smk_status_sig_res.tsv",
        help="Path to the OTU-level smoking association results.",
    )
    parser.add_argument(
        "--age-results",
        default="data/age_sig_res.tsv",
        help="Path to the OTU-level age association results.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory for all supplementary analysis outputs.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260319,
        help="Random seed for rarefaction subsampling.",
    )
    return parser.parse_args()


def clean_metadata(path: Path) -> pd.DataFrame:
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
    meta["smoke_binary"] = np.where(
        meta["smk_status"].eq("everyday"),
        "everyday",
        np.where(meta["smk_status"].eq("never_smoker"), "never_smoker", np.nan),
    )
    meta["smoke_everyday"] = meta["smoke_binary"].eq("everyday").astype(float)
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


def compute_rarefaction(table, sample_ids: list[str], output_dir: Path, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    sample_ids = list(sample_ids)
    chosen = sorted(rng.choice(sample_ids, size=min(120, len(sample_ids)), replace=False))
    selected = table.filter(chosen, axis="sample", inplace=False)
    counts = selected.matrix_data.toarray().T
    depths = [500, 1000, 2000, 4000, 6000, 8000, 10000]
    records = []
    for sid, vec in zip(chosen, counts):
        vec = vec.astype(int)
        for depth in depths:
            np.random.seed(int(rng.integers(0, 1_000_000_000)))
            sub = subsample_counts(vec, n=depth, replace=False)
            observed = int((sub > 0).sum())
            records.append({"sample_id": sid, "depth": depth, "observed_otus": observed})
    rare = pd.DataFrame(records)
    summary = (
        rare.groupby("depth")["observed_otus"]
        .agg(["mean", "std", "median", "min", "max"])
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    for sid, grp in rare.groupby("sample_id"):
        ax.plot(grp["depth"], grp["observed_otus"], color="#c7d4ea", alpha=0.15, linewidth=1)
    ax.plot(summary["depth"], summary["mean"], color="#1f4e79", linewidth=3, label="Mean observed OTUs")
    ax.fill_between(
        summary["depth"],
        summary["mean"] - summary["std"],
        summary["mean"] + summary["std"],
        color="#8fb4d9",
        alpha=0.25,
        label="Mean +/- 1 SD",
    )
    ax.set_xlabel("Rarefaction depth")
    ax.set_ylabel("Observed OTUs")
    ax.set_title("Rarefaction curve on the 6,676-sample analysis subset")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_dir / "supp_rarefaction_curve.png", dpi=300)
    fig.savefig(output_dir / "supp_rarefaction_curve.pdf")
    plt.close(fig)

    depth_summary = pd.DataFrame(
        {
            "analysis_subset_n": [len(sample_ids)],
            "min_reads": [10000],
            "median_reads": [10000],
            "max_reads": [10000],
        }
    )
    depth_summary.to_csv(output_dir / "sequencing_depth_summary.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "rarefaction_summary.tsv", sep="\t", index=False)
    return rare


def load_overlap_features(smoking_results: Path, age_results: Path) -> pd.DataFrame:
    smoking_paths = [
        smoking_results,
        smoking_results.parent / "smk_amount_categ_sig_res.tsv",
        smoking_results.parent / "smk_y_categ_sig_res.tsv",
    ]
    smoking_frames = []
    for path in smoking_paths:
        df = pd.read_csv(path, sep="\t")
        df["feature"] = df["feature"].astype(str)
        df = df[df["qval"] < 0.05].copy()
        smoking_frames.append(df)
    smk = pd.concat(smoking_frames, ignore_index=True)
    age = pd.read_csv(age_results, sep="\t")
    age["feature"] = age["feature"].astype(str)
    smk = smk.sort_values(["qval", "pval"]).groupby("feature", as_index=False).first()
    smk = smk[["feature", "metadata", "value", "coef", "qval", "Phylum", "Family", "Genus", "Species"]]
    smk = smk.rename(
        columns={
            "metadata": "smoking_metadata",
            "value": "smoking_value",
            "coef": "smoking_coef",
            "qval": "smoking_qval",
        }
    )
    age = age[["feature", "coef", "qval", "Phylum", "Family", "Genus", "Species"]]
    age = age.rename(columns={"coef": "age_coef", "qval": "age_qval"})
    shared = smk.merge(age, on=["feature", "Phylum", "Family", "Genus", "Species"], how="inner")
    shared["direction_concordant"] = np.sign(shared["smoking_coef"]) == np.sign(shared["age_coef"])
    return shared


def male_only_sensitivity(
    table,
    meta: pd.DataFrame,
    overlap: pd.DataFrame,
    output_dir: Path,
) -> dict:
    male = meta[
        meta["gender"].eq("m")
        & meta["smoke_binary"].notna()
        & meta["age"].notna()
        & meta["bmi"].notna()
        & meta["bristol"].notna()
        & meta["age_categ"].notna()
    ].copy()
    overlap = overlap.copy()
    concordant = overlap[overlap["direction_concordant"]].copy()
    concordant["orientation"] = np.sign(concordant["age_coef"])
    feature_ids = concordant["feature"].astype(str).tolist()
    counts = feature_table(table, feature_ids)
    rel = relative_abundance(counts).reindex(columns=male["ID"])
    transformed = arcsin_sqrt_transform(rel).T
    oriented = transformed.copy()
    orient_map = concordant.set_index("feature")["orientation"].to_dict()
    for feature, sign in orient_map.items():
        oriented[feature] = oriented[feature] * float(sign)
    score = oriented.apply(lambda col: (col - col.mean()) / col.std(ddof=0))
    male["pro_aging_score"] = score.mean(axis=1).values

    quartile_records = []
    for quartile in male["age_categ"].cat.categories:
        sub = male[male["age_categ"] == quartile].copy()
        x = build_design(
            sub,
            "1 + smoke_everyday + age + bmi + C(bristol_cat) + C(county)",
        )
        y = sub[["pro_aging_score"]].rename(columns={"pro_aging_score": "module_score"})
        res = matrix_ols(y, x, "smoke_everyday").iloc[0].to_dict()
        res["age_categ"] = quartile
        res["male_everyday_n"] = int((sub["smoke_binary"] == "everyday").sum())
        res["male_never_n"] = int((sub["smoke_binary"] == "never_smoker").sum())
        quartile_records.append(res)
    quartile_df = bh_qvalues(pd.DataFrame(quartile_records))
    quartile_df.to_csv(output_dir / "male_only_module_effects_by_age.tsv", sep="\t", index=False)

    otu_results = []
    x_full = build_design(
        male,
        "1 + smoke_everyday + age + bmi + C(bristol_cat) + C(county)",
    )
    y_full = transformed[feature_ids]
    male_coef = bh_qvalues(matrix_ols(y_full, x_full, "smoke_everyday"))
    male_coef = male_coef.rename(
        columns={
            "coef": "male_smoking_coef",
            "stderr": "male_stderr",
            "pval": "male_pval",
            "qval": "male_qval",
            "ci_low": "male_ci_low",
            "ci_high": "male_ci_high",
        }
    )
    otu_results = concordant.merge(male_coef, on="feature", how="left")
    otu_results["same_direction_as_full"] = np.sign(otu_results["male_smoking_coef"]) == np.sign(
        otu_results["smoking_coef"]
    )
    otu_results.to_csv(output_dir / "male_only_overlap_otu_concordance.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    qplot = quartile_df.copy()
    palette = ["#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"]
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

    scatter = otu_results.copy()
    axes[1].scatter(
        scatter["smoking_coef"],
        scatter["male_smoking_coef"],
        c=np.where(scatter["same_direction_as_full"], "#2a9d8f", "#d1495b"),
        alpha=0.8,
        edgecolor="white",
        linewidth=0.4,
    )
    lim = np.nanmax(np.abs(scatter[["smoking_coef", "male_smoking_coef"]].to_numpy())) * 1.1
    axes[1].plot([-lim, lim], [-lim, lim], linestyle="--", color="#9aa0a6", linewidth=1)
    axes[1].axhline(0, color="#d0d7de", linewidth=1)
    axes[1].axvline(0, color="#d0d7de", linewidth=1)
    axes[1].set_xlim(-lim, lim)
    axes[1].set_ylim(-lim, lim)
    axes[1].set_xlabel("Full-sample OTU smoking coefficient")
    axes[1].set_ylabel("Male-only OTU smoking coefficient")
    axes[1].set_title("Concordance across concordant overlap OTUs")

    fig.tight_layout()
    fig.savefig(output_dir / "supp_male_only_sensitivity.png", dpi=300)
    fig.savefig(output_dir / "supp_male_only_sensitivity.pdf")
    plt.close(fig)

    return {
        "male_analysis_n": int(male.shape[0]),
        "male_everyday_n": int((male["smoke_binary"] == "everyday").sum()),
        "male_never_n": int((male["smoke_binary"] == "never_smoker").sum()),
        "shared_overlap_otus": int(overlap.shape[0]),
        "concordant_overlap_otus": int(concordant.shape[0]),
        "same_direction_male_otus": int(otu_results["same_direction_as_full"].sum()),
        "male_nominal_p_lt_0_05": int((otu_results["male_pval"] < 0.05).sum()),
    }


def run_taxonomic_robustness(table, meta: pd.DataFrame, level: str, output_dir: Path) -> dict:
    counts = collapse_table(table, level)
    ids = meta["ID"].tolist()
    counts = counts.loc[:, ids]
    rel = relative_abundance(counts)

    smoking_meta = meta[
        meta["smoke_binary"].notna()
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
    smoking_res.to_csv(output_dir / f"{level}_smoking_results.tsv", sep="\t", index=False)

    age_meta = meta[
        meta["smk_status"].eq("never_smoker")
        & meta["age"].notna()
        & meta["bmi"].notna()
        & meta["bristol"].notna()
    ].copy()
    age_ids = age_meta["ID"].tolist()
    age_rel = prevalence_filter(rel.loc[:, age_ids]).T
    age_y = arcsin_sqrt_transform(age_rel.T).T
    age_x = build_design(
        age_meta,
        "1 + age + bmi + C(gender) + C(bristol_cat) + C(county)",
    )
    age_res = bh_qvalues(matrix_ols(age_y, age_x, "age"))
    age_res = age_res.sort_values(["qval", "pval"]).reset_index(drop=True)
    age_res.to_csv(output_dir / f"{level}_age_results.tsv", sep="\t", index=False)

    shared = smoking_res.merge(age_res, on="feature", suffixes=("_smoking", "_age"))
    shared["direction_concordant"] = np.sign(shared["coef_smoking"]) == np.sign(shared["coef_age"])
    shared["both_fdr_sig"] = (shared["qval_smoking"] < 0.05) & (shared["qval_age"] < 0.05)
    shared.to_csv(output_dir / f"{level}_shared_results.tsv", sep="\t", index=False)

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
        top = plot_df.assign(rank=plot_df["qval_smoking"] + plot_df["qval_age"]).sort_values("rank").head(10)
        for _, row in top.iterrows():
            ax.text(row["coef_smoking"], row["coef_age"], row["feature"], fontsize=9)
    ax.axhline(0, color="#d0d7de", linewidth=1)
    ax.axvline(0, color="#d0d7de", linewidth=1)
    ax.set_xlabel("Smoking coefficient (everyday vs never)")
    ax.set_ylabel("Age coefficient (within never smokers)")
    ax.set_title(f"{level.capitalize()}-level overlap")
    fig.tight_layout()
    fig.savefig(output_dir / f"supp_{level}_overlap_scatter.png", dpi=300)
    fig.savefig(output_dir / f"supp_{level}_overlap_scatter.pdf")
    plt.close(fig)

    return {
        "level": level,
        "smoking_tested_features": int(smoking_res.shape[0]),
        "age_tested_features": int(age_res.shape[0]),
        "smoking_fdr_sig": int((smoking_res["qval"] < 0.05).sum()),
        "age_fdr_sig": int((age_res["qval"] < 0.05).sum()),
        "shared_both_fdr_sig": int(shared["both_fdr_sig"].sum()),
        "shared_concordant": int((shared["both_fdr_sig"] & shared["direction_concordant"]).sum()),
    }


def write_summary(
    output_dir: Path,
    depth_summary: dict,
    male_summary: dict,
    tax_summaries: list[dict],
) -> None:
    lines = [
        "# Revision Analysis Summary",
        "",
        "## Sequencing depth / rarefaction",
        f"- Analysis subset size: {depth_summary['analysis_subset_n']}",
        f"- Public BIOM table depth: min = {depth_summary['min_reads']}, median = {depth_summary['median_reads']}, max = {depth_summary['max_reads']}",
        "",
        "## Male-only sensitivity analysis",
        f"- Male complete-case subset for the sensitivity model: {male_summary['male_analysis_n']}",
        f"- Everyday smokers: {male_summary['male_everyday_n']}",
        f"- Never smokers: {male_summary['male_never_n']}",
        f"- Shared OTUs from the original OTU-level overlap: {male_summary['shared_overlap_otus']}",
        f"- Concordant overlap OTUs used for the pro-aging module: {male_summary['concordant_overlap_otus']}",
        f"- OTUs retaining the same smoking-effect direction in male-only models: {male_summary['same_direction_male_otus']}",
        f"- OTUs with nominal P < 0.05 in male-only models: {male_summary['male_nominal_p_lt_0_05']}",
        "",
        "## Genus/family robustness",
    ]
    for entry in tax_summaries:
        lines.extend(
            [
                f"### {entry['level'].capitalize()}",
                f"- Tested in smoking models: {entry['smoking_tested_features']}",
                f"- Tested in age models: {entry['age_tested_features']}",
                f"- Smoking FDR-significant features: {entry['smoking_fdr_sig']}",
                f"- Age FDR-significant features: {entry['age_fdr_sig']}",
                f"- Shared features significant in both models: {entry['shared_both_fdr_sig']}",
                f"- Shared features with concordant directions: {entry['shared_concordant']}",
                "",
            ]
        )
    (output_dir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    meta = clean_metadata(Path(args.metadata_path))
    table = subset_biom(Path(args.biom_path), meta["ID"])

    rare = compute_rarefaction(table, meta["ID"].tolist(), output_dir, seed=args.seed)
    overlap = load_overlap_features(Path(args.smoking_results), Path(args.age_results))
    male_summary = male_only_sensitivity(table, meta, overlap, output_dir)
    tax_summaries = [
        run_taxonomic_robustness(table, meta, "family", output_dir),
        run_taxonomic_robustness(table, meta, "genus", output_dir),
    ]
    depth_summary = pd.read_csv(output_dir / "sequencing_depth_summary.tsv", sep="\t").iloc[0].to_dict()
    write_summary(output_dir, depth_summary, male_summary, tax_summaries)

    overlap.to_csv(output_dir / "otu_overlap_from_original_results.tsv", sep="\t", index=False)
    meta.to_csv(output_dir / "analysis_metadata_snapshot.tsv", sep="\t", index=False)
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
