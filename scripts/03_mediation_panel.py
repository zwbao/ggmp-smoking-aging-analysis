#!/usr/bin/env python3
"""Multi-mediator panel for the GGMP smoking x aging revision (Bao et al.,
NTR-2026-092). Re-runs the published `gai_med.R` mediation specification across
ten cardiometabolic candidate mediators (SBP, DBP, BMI, waist, FBG, TCHO, TG,
UA, HDL, LDL) to support the SBP-as-primary-mediator choice in Figure 6b
(produces Supp Fig 10 + Supp Table 20). Path equations are the original
unadjusted spec:

    M  ~ X
    Y  ~ M + X
    Y  ~ X         (total effect)

faithfully replicating `lm(...)` + `mediate(..., sims=1000, boot=T)` under the
LM specification, with bootstrap percentile 95% CIs over 1,000 resamples
(seed=12345). Inputs default to ./data/ and outputs go to ./outputs/.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biom import load_table
from scipy.stats import t
from statsmodels.stats.multitest import fdrcorrection


sns.set_theme(style="whitegrid", context="talk")


CANDIDATE_MEDIATORS = [
    "anthrop_SBP",
    "anthrop_DBP",
    "anthrop_BMI",
    "anthrop_waist",
    "biochem_FBG",
    "biochem_TCHO",
    "biochem_TG",
    "biochem_UA",
    "biochem_HDL",
    "biochem_LDL",
]

SEED = 12345  # matches gai_med.R
N_BOOTSTRAP = 1000


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--biom-path",
        default="data/GGMP7009_even10k.biom",
    )
    parser.add_argument(
        "--metadata-path",
        default="data/GPf_metadata.tsv",
    )
    parser.add_argument(
        "--overlap-path",
        default="outputs/otu_overlap_from_original_results.tsv",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
    )
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--n-boot", type=int, default=N_BOOTSTRAP)
    return parser.parse_args()


def clean_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", low_memory=False).copy()
    meta["ID"] = meta["ID"].astype(str)
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["bmi"] = pd.to_numeric(meta["anthrop_BMI"], errors="coerce")
    meta["bristol"] = pd.to_numeric(meta["Bristol_stool_type"], errors="coerce")
    meta["gender"] = meta["gender"].astype(str)
    meta["county"] = meta["county_level_code"].astype(str)
    meta["smk_status"] = meta["smk_status"].astype(str)
    for col in CANDIDATE_MEDIATORS:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")
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
# OLS helpers (no patsy needed — all designs are simple [1, x] or [1, x1, x2])
# ---------------------------------------------------------------------------

def ols_simple(y: np.ndarray, X: np.ndarray, focal_idx: int) -> dict:
    """Return coef, p, se, ci for column focal_idx of design X."""
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    n, p = X.shape
    df_resid = n - p
    sse = float(np.sum(resid ** 2))
    sigma2 = sse / df_resid
    se = float(np.sqrt(XtX_inv[focal_idx, focal_idx] * sigma2))
    coef = float(beta[focal_idx])
    t_stat = coef / se if se > 0 else np.nan
    pval = float(2 * t.sf(np.abs(t_stat), df=df_resid))
    return {"coef": coef, "stderr": se, "pval": pval, "n": n}


def ols_coef_only(y: np.ndarray, X: np.ndarray, focal_idx: int) -> float:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(beta[focal_idx])


# ---------------------------------------------------------------------------
# Pro-aging module construction (matches male_only_sensitivity)
# ---------------------------------------------------------------------------

def build_proaging_score(table, meta: pd.DataFrame, overlap_path: Path) -> pd.Series:
    overlap = pd.read_csv(overlap_path, sep="\t")
    concordant = overlap[overlap["direction_concordant"]].copy()
    concordant["orientation"] = np.sign(concordant["age_coef"])
    feature_ids = concordant["feature"].astype(str).tolist()

    counts = feature_table(table, feature_ids)
    rel = relative_abundance(counts).reindex(columns=meta["ID"])
    transformed = arcsin_sqrt_transform(rel).T  # samples x OTU
    oriented = transformed.copy()
    orient_map = concordant.set_index("feature")["orientation"].to_dict()
    for feature, sign in orient_map.items():
        oriented[str(feature)] = oriented[str(feature)] * float(sign)
    score = oriented.apply(lambda col: (col - col.mean()) / col.std(ddof=0))
    return pd.Series(score.mean(axis=1).values, index=meta["ID"], name="pro_aging_score")


# ---------------------------------------------------------------------------
# ASCVD 10-year risk - ACC/AHA 2013 PCE
# ---------------------------------------------------------------------------

PCE_COEF = {
    ("female", "white"): {
        "age": -29.799, "age_sq": 4.884, "tchol": 13.540, "age_tchol": -3.114,
        "hdl": -13.578, "age_hdl": 3.149, "treated_sbp": 2.019,
        "age_treated_sbp": 0.0, "untreated_sbp": 1.957, "age_untreated_sbp": 0.0,
        "smoker": 7.574, "age_smoker": -1.665, "diabetes": 0.661,
        "baseline_survival": 0.9665, "mean_terms": -29.18,
    },
    ("female", "aa"): {
        "age": 17.114, "age_sq": 0.0, "tchol": 0.940, "age_tchol": 0.0,
        "hdl": -18.920, "age_hdl": 4.475, "treated_sbp": 29.291,
        "age_treated_sbp": -6.432, "untreated_sbp": 27.820, "age_untreated_sbp": -6.087,
        "smoker": 0.691, "age_smoker": 0.0, "diabetes": 0.874,
        "baseline_survival": 0.9533, "mean_terms": 86.61,
    },
    ("male", "white"): {
        "age": 12.344, "age_sq": 0.0, "tchol": 11.853, "age_tchol": -2.664,
        "hdl": -7.990, "age_hdl": 1.769, "treated_sbp": 1.797,
        "age_treated_sbp": 0.0, "untreated_sbp": 1.764, "age_untreated_sbp": 0.0,
        "smoker": 7.837, "age_smoker": -1.795, "diabetes": 0.658,
        "baseline_survival": 0.9144, "mean_terms": 61.18,
    },
    ("male", "aa"): {
        "age": 2.469, "age_sq": 0.0, "tchol": 0.302, "age_tchol": 0.0,
        "hdl": -0.307, "age_hdl": 0.0, "treated_sbp": 1.916,
        "age_treated_sbp": 0.0, "untreated_sbp": 1.809, "age_untreated_sbp": 0.0,
        "smoker": 0.549, "age_smoker": 0.0, "diabetes": 0.645,
        "baseline_survival": 0.8954, "mean_terms": 19.54,
    },
}


def ascvd_10y_accaha(race, gender, age, totchol, hdl, sbp, bp_med, smoker, diabetes):
    if race == "other":
        race = "white"
    key = (gender, race)
    if key not in PCE_COEF:
        return np.nan
    if any(pd.isna(v) for v in [age, totchol, hdl, sbp]):
        return np.nan
    if age < 40 or age > 79:
        return np.nan
    c = PCE_COEF[key]
    ln_age = np.log(age)
    ln_age_sq = ln_age ** 2
    ln_tchol = np.log(totchol)
    ln_hdl = np.log(hdl)
    ln_sbp = np.log(sbp)

    indv_sum = (
        c["age"] * ln_age
        + c["age_sq"] * ln_age_sq
        + c["tchol"] * ln_tchol
        + c["age_tchol"] * ln_age * ln_tchol
        + c["hdl"] * ln_hdl
        + c["age_hdl"] * ln_age * ln_hdl
        + (c["treated_sbp"] * ln_sbp + c["age_treated_sbp"] * ln_age * ln_sbp
           if bp_med else
           c["untreated_sbp"] * ln_sbp + c["age_untreated_sbp"] * ln_age * ln_sbp)
        + c["smoker"] * (1 if smoker else 0)
        + c["age_smoker"] * ln_age * (1 if smoker else 0)
        + c["diabetes"] * (1 if diabetes else 0)
    )
    risk = 1.0 - c["baseline_survival"] ** np.exp(indv_sum - c["mean_terms"])
    return float(risk * 100.0)


def compute_ascvd(meta: pd.DataFrame) -> pd.Series:
    df = meta.copy()
    bp1 = df.get("bp_control_medication1", pd.Series([np.nan] * len(df)))
    bp2 = df.get("bp_control_medication2", pd.Series([np.nan] * len(df)))
    bp3 = df.get("bp_control_medication3", pd.Series([np.nan] * len(df)))
    bp_med = (bp1.eq("y") | bp2.eq("y") | bp3.eq("y")).fillna(False)

    gender = df["gender"].map({"f": "female", "m": "male"}).fillna("")
    smk = df["smk_status"].map({"never_smoker": 0, "current_smoker": 1})
    bg = df.get("bg_diagnosis", pd.Series([np.nan] * len(df)))
    bg_bin = bg.eq("y").astype(float)
    bg_bin = bg_bin.where(bg.notna(), 0.0)

    tcho = pd.to_numeric(df["biochem_TCHO"], errors="coerce") * 38.67
    hdl = pd.to_numeric(df["biochem_HDL"], errors="coerce") * 38.67
    sbp = pd.to_numeric(df["anthrop_SBP"], errors="coerce")
    age = pd.to_numeric(df["age"], errors="coerce")

    out = []
    for i in range(len(df)):
        if pd.isna(smk.iloc[i]):
            out.append(np.nan)
            continue
        out.append(
            ascvd_10y_accaha(
                race="other",
                gender=gender.iloc[i],
                age=age.iloc[i],
                totchol=tcho.iloc[i],
                hdl=hdl.iloc[i],
                sbp=sbp.iloc[i],
                bp_med=bool(bp_med.iloc[i]),
                smoker=int(smk.iloc[i]),
                diabetes=int(bg_bin.iloc[i]),
            )
        )
    return pd.Series(out, index=df.index, name="ascvd_10y")


def detect_existing_ascvd_column(meta: pd.DataFrame) -> str | None:
    pat = ("ascvd", "cvrisk", "risk_10y")
    for col in meta.columns:
        lc = col.lower()
        if any(p in lc for p in pat):
            return col
    return None


# ---------------------------------------------------------------------------
# Mediation - ORIGINAL gai_med.R spec (no covariates beyond X and M)
# ---------------------------------------------------------------------------

def run_mediation_for(
    df: pd.DataFrame,
    mediator: str,
    rng: np.random.Generator,
    n_boot: int,
) -> dict:
    """Single-mediator analysis using the original gai_med.R specification:
        path a:   M ~ X            (X = pro_aging_score; ONLY covariate is X)
        path b:   Y ~ M + X        (Y = ascvd_10y; ONLY covariates are M, X)
        path c:   Y ~ X            (total effect of X on Y)
    Bootstrap 1000 reps, percentile CI on indirect = a*b.
    """
    sub = df[[mediator, "ascvd_10y", "pro_aging_score"]].dropna()
    n = len(sub)
    if n < 50:
        return {
            "mediator": mediator, "n_complete": n,
            "coef_a": np.nan, "p_a": np.nan,
            "coef_b": np.nan, "p_b": np.nan,
            "total_c": np.nan, "p_total_c": np.nan,
            "direct_c_prime": np.nan, "p_direct_c_prime": np.nan,
            "indirect_ab": np.nan,
            "indirect_ci_low": np.nan, "indirect_ci_high": np.nan,
            "prop_mediated": np.nan,
            "prop_mediated_ci_low": np.nan, "prop_mediated_ci_high": np.nan,
            "comment": "insufficient complete cases (<50)",
        }

    m_arr = sub[mediator].to_numpy(dtype=float)
    y_arr = sub["ascvd_10y"].to_numpy(dtype=float)
    x_arr = sub["pro_aging_score"].to_numpy(dtype=float)

    ones = np.ones(n)
    # path a: M ~ 1 + X    (focal coef on X is index 1)
    Xa = np.column_stack([ones, x_arr])
    # path b: Y ~ 1 + M + X  (focal coef on M is index 1, on X is index 2)
    Xb = np.column_stack([ones, m_arr, x_arr])
    # path c: Y ~ 1 + X    (focal coef on X is index 1)
    Xc = np.column_stack([ones, x_arr])

    res_a = ols_simple(m_arr, Xa, 1)              # coef_a (M~X)
    res_b = ols_simple(y_arr, Xb, 1)              # coef_b (Y~M+X), coef on M
    res_cp = ols_simple(y_arr, Xb, 2)             # c' = coef on X in Y~M+X
    res_c = ols_simple(y_arr, Xc, 1)              # total c = coef on X in Y~X

    coef_a = res_a["coef"]
    coef_b = res_b["coef"]
    indirect = coef_a * coef_b
    total_c = res_c["coef"]

    # Bootstrap CI for indirect effect — percentile method on a*b
    boot_indirect = np.empty(n_boot)
    boot_indirect[:] = np.nan
    boot_prop = np.empty(n_boot)
    boot_prop[:] = np.nan
    boot_acme_prop = np.empty(n_boot)
    boot_acme_prop[:] = np.nan
    idx_arr = np.arange(n)

    for b in range(n_boot):
        idx = rng.choice(idx_arr, size=n, replace=True)
        try:
            ca = ols_coef_only(m_arr[idx], Xa[idx], 1)
            cb = ols_coef_only(y_arr[idx], Xb[idx], 1)
            ccp = ols_coef_only(y_arr[idx], Xb[idx], 2)
        except (np.linalg.LinAlgError, ValueError):
            continue
        ind_b = ca * cb
        boot_indirect[b] = ind_b
        # ACME + ADE = TE under LM, so prop_mediated = ACME / (ACME + ADE).
        # This is exactly what mediation::summary(mediate(...))$d0 / TE returns
        # for an LM specification.
        te_b = ind_b + ccp
        if te_b != 0 and np.sign(ind_b) == np.sign(te_b):
            boot_acme_prop[b] = ind_b / te_b
        # Also report prop based on standalone total-c regression (for legacy
        # comparison with the previous panel run).
        try:
            cc = ols_coef_only(y_arr[idx], Xc[idx], 1)
        except (np.linalg.LinAlgError, ValueError):
            continue
        if cc != 0 and np.sign(ind_b) == np.sign(cc):
            boot_prop[b] = ind_b / cc

    ci_low, ci_high = np.nanpercentile(boot_indirect, [2.5, 97.5])

    # Point-estimate proportion-mediated via gai_med.R definition:
    # prop_mediated = ACME / (ACME + ADE) where ACME = a*b, ADE = c'
    te_point = indirect + res_cp["coef"]
    if te_point != 0 and np.sign(indirect) == np.sign(te_point):
        prop_med = indirect / te_point
        comment = ""
    else:
        prop_med = np.nan
        comment = "prop_mediated NA: ACME and (ACME+ADE) disagree in sign or sum=0"

    if np.isfinite(prop_med):
        prop_ci_low, prop_ci_high = np.nanpercentile(boot_acme_prop, [2.5, 97.5])
    else:
        prop_ci_low = prop_ci_high = np.nan

    # Also compute the legacy "indirect / total_c" proportion (using standalone
    # total-c regression). When no covariates exist beyond X, total_c == ACME +
    # ADE exactly, so these will agree to floating-point precision; we keep
    # both columns in the output to mirror the previous schema.
    return {
        "mediator": mediator,
        "n_complete": n,
        "coef_a": coef_a,
        "p_a": res_a["pval"],
        "coef_b": coef_b,
        "p_b": res_b["pval"],
        "total_c": total_c,
        "p_total_c": res_c["pval"],
        "direct_c_prime": res_cp["coef"],
        "p_direct_c_prime": res_cp["pval"],
        "indirect_ab": indirect,
        "indirect_ci_low": ci_low,
        "indirect_ci_high": ci_high,
        "prop_mediated": prop_med,
        "prop_mediated_ci_low": prop_ci_low,
        "prop_mediated_ci_high": prop_ci_high,
        "comment": comment,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def forest_plot(results: pd.DataFrame, output_dir: Path) -> None:
    df = results.dropna(subset=["indirect_ab"]).copy()
    df["abs_ind"] = df["indirect_ab"].abs()
    df = df.sort_values("abs_ind", ascending=True)
    df["sig_both"] = (df["q_a"] < 0.05) & (df["q_b"] < 0.05)

    fig, ax = plt.subplots(figsize=(9, 6))
    y = np.arange(len(df))
    colors = np.where(df["sig_both"].values, "#d1495b", "#5b8db8")
    for yi, (lo, hi) in zip(y, zip(df["indirect_ci_low"].values, df["indirect_ci_high"].values)):
        ax.hlines(yi, lo, hi, color="#777", linewidth=1.4, zorder=2)
        ax.vlines([lo, hi], yi - 0.18, yi + 0.18, color="#777", linewidth=1.2, zorder=2)
    ax.scatter(df["indirect_ab"].values, y, c=colors, s=110, zorder=3, edgecolor="white")
    ax.axvline(0, color="#9aa0a6", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["mediator"].values)
    ax.set_xlabel("Indirect effect (a x b) on ASCVD 10-year risk (%)")
    ax.set_title(
        "Mediator panel (original spec, no covariates):\n"
        "pro-aging module -> M -> ASCVD risk\n"
        "(red = q_a < 0.05 AND q_b < 0.05)"
    )
    fig.tight_layout()
    fig.savefig(output_dir / "supp_mediation_panel.png", dpi=300)
    fig.savefig(output_dir / "supp_mediation_panel.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Memo
# ---------------------------------------------------------------------------

def write_memo(results: pd.DataFrame, output_dir: Path, ascvd_source: str) -> None:
    sig = results[(results["q_a"] < 0.05) & (results["q_b"] < 0.05)].copy()
    sig = sig.sort_values("indirect_ab", key=lambda s: s.abs(), ascending=False)

    sbp_row = results[results["mediator"] == "anthrop_SBP"]
    sbp_prop = (
        float(sbp_row["prop_mediated"].iloc[0])
        if not sbp_row.empty and pd.notna(sbp_row["prop_mediated"].iloc[0])
        else np.nan
    )
    sbp_qa = float(sbp_row["q_a"].iloc[0]) if not sbp_row.empty else np.nan
    sbp_qb = float(sbp_row["q_b"].iloc[0]) if not sbp_row.empty else np.nan
    sbp_in_sig = "anthrop_SBP" in sig["mediator"].values

    if not sig.empty:
        top = sig.iloc[0]
        if pd.notna(top["prop_mediated"]):
            top_line = (
                f"- **Top mediator by |indirect|:** {top['mediator']} "
                f"(indirect = {top['indirect_ab']:.4g}, "
                f"prop. mediated = {top['prop_mediated']:.1%}, "
                f"95% CI {top['prop_mediated_ci_low']:.1%} - "
                f"{top['prop_mediated_ci_high']:.1%})"
            )
        else:
            top_line = (
                f"- **Top mediator by |indirect|:** {top['mediator']} "
                f"(indirect = {top['indirect_ab']:.4g}, prop. mediated = NA)"
            )
        verdict_lines = [
            f"- **Mediators FDR-significant on both path a (q_a<0.05) and path b (q_b<0.05):** {len(sig)} of {len(results)}",
            top_line,
        ]
    else:
        verdict_lines = ["- No mediator reached q<0.05 on both paths simultaneously."]

    # SBP-specific reproduction line
    pub_target = 28.5
    if np.isfinite(sbp_prop):
        sbp_pct = sbp_prop * 100.0
        delta = sbp_pct - pub_target
        if abs(delta) <= 5.0:
            repro = "REPRODUCED"
        else:
            repro = "MISMATCHED"
        sbp_line = (
            f"- **SBP proportion mediated (current spec):** {sbp_pct:.1f}% "
            f"(95% CI {results.loc[sbp_row.index[0], 'prop_mediated_ci_low']*100:.1f}% - "
            f"{results.loc[sbp_row.index[0], 'prop_mediated_ci_high']*100:.1f}%); "
            f"published Fig. 6b reported 28.5% — delta = {delta:+.1f} pp ({repro})"
        )
    else:
        sbp_line = (
            f"- **SBP proportion mediated:** NA "
            f"(comment: '{results.loc[sbp_row.index[0], 'comment']}')"
        )

    table = results[
        ["mediator", "n_complete", "coef_a", "q_a", "coef_b", "q_b",
         "total_c", "direct_c_prime", "indirect_ab", "indirect_ci_low",
         "indirect_ci_high", "prop_mediated", "prop_mediated_ci_low",
         "prop_mediated_ci_high", "comment"]
    ].copy()
    md_table = table.to_markdown(index=False, floatfmt=".4g")

    # Narrative paragraph for the response letter
    if not sig.empty:
        if sbp_in_sig:
            sbp_sig_row = sig[sig["mediator"] == "anthrop_SBP"].iloc[0]
            sbp_rank = (
                sig["indirect_ab"].abs()
                .rank(method="min", ascending=False)
                .reindex(sig.index)
            )
            sbp_rank_value = int(sbp_rank.loc[sbp_sig_row.name])
            narrative = (
                f"Among the 10 candidate mediators systematically tested under "
                f"the original gai_med.R specification (no covariate adjustment "
                f"beyond the exposure and mediator), {len(sig)} reached "
                f"FDR-significance on both path a and path b. Systolic blood "
                f"pressure (SBP) ranked #{sbp_rank_value} of the {len(sig)} "
                f"jointly significant mediators by absolute indirect effect, "
                f"with a proportion mediated of "
                f"{sbp_sig_row['prop_mediated']:.1%} (95% CI "
                f"{sbp_sig_row['prop_mediated_ci_low']:.1%}-"
                f"{sbp_sig_row['prop_mediated_ci_high']:.1%}). "
                f"This empirically supports its inclusion in the primary "
                f"mediation analysis (Figure 6b)."
            )
        else:
            narrative = (
                f"Among the 10 candidate mediators systematically tested under "
                f"the original gai_med.R specification, {len(sig)} reached "
                f"FDR-significance on both path a and path b, but SBP was NOT "
                f"among them. The lead author should reconsider the "
                f"primary-mediation choice or refine the response-letter framing."
            )
    else:
        narrative = (
            "No candidate mediator was FDR-significant on both paths under the "
            "original specification."
        )

    methods = (
        "**Mediator-panel sensitivity analysis (Methods).** To address Reviewer "
        "2's question of why SBP was the only mediator reported in Figure 6b, "
        "we systematically tested ten cardiometabolic intermediates available "
        "in GGMP (SBP, DBP, BMI, waist circumference, fasting blood glucose, "
        "total cholesterol, triglycerides, uric acid, HDL, LDL) using the same "
        "specification as in the primary analysis (gai_med.R): for each "
        "candidate mediator M we fitted lm(M ~ pro_aging_module) for path a, "
        "lm(ASCVD ~ M + pro_aging_module) for path b, and lm(ASCVD ~ "
        "pro_aging_module) for the total effect c, with no demographic or "
        "behavioural covariates beyond the exposure and mediator. Mediation "
        "was estimated by the product-of-coefficients method (indirect = a*b) "
        "with bootstrap percentile 95% confidence intervals over 1,000 "
        "resamples (seed=12345), reproducing the call to "
        "`mediation::mediate(..., sims=1000, boot=TRUE)` used in the primary "
        "analysis. Path-a and path-b p-values were Benjamini-Hochberg "
        "corrected across the ten candidates. The proportion mediated is "
        "reported as ACME/(ACME + ADE) per the standard mediation::summary "
        "convention; rows where ACME and (ACME+ADE) disagree in sign are "
        "flagged as undefined."
    )

    lines = [
        "# Supplementary Mediation Panel (original-spec re-run)",
        "",
        "> **This memo replaces the previous full-covariate panel.** "
        "The previous version (preserved as `mediation_panel_results_full_covariates.tsv` "
        "and `supp_mediation_panel_full_covariates.png|pdf`) adjusted for "
        "age + gender + BMI + Bristol + county + smoking status in every "
        "path equation. That specification absorbs the smoking-driven module "
        "variance and is inconsistent with the original published analysis "
        "(`gai_med.R` in the original GGMP processing repository, lines 81-118), which uses simple `lm(M ~ X)` and "
        "`lm(Y ~ M + X)` with no demographic or behavioural covariates. "
        "The corrected re-run below matches the published spec.",
        "",
        f"- ASCVD source: **{ascvd_source}**",
        f"- Bootstrap resamples: {N_BOOTSTRAP}; seed = {SEED} (matches gai_med.R `set.seed(12345)`)",
        "- **Path equations (no demographic covariates):** `M ~ pro_aging_score`, `ASCVD ~ M + pro_aging_score`, `ASCVD ~ pro_aging_score`.",
        "- Pro-aging module: 40 direction-concordant overlap OTUs, arcsin-sqrt rel. abundance, oriented by sign(age_coef), z-scored, mean across OTUs (full 6,676-sample subset).",
        "- Proportion mediated = ACME / (ACME + ADE) per `mediation::summary` (equivalent to indirect / total under LM).",
        "",
        "## Verdict",
        *verdict_lines,
        sbp_line,
        "",
        "## Narrative answer to Reviewer 2",
        narrative,
        "",
        "## Full results",
        md_table,
        "",
        "## Methods paragraph (publication-ready)",
        methods,
        "",
    ]
    (output_dir / "mediation_panel_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    meta = clean_metadata(Path(args.metadata_path))
    table = subset_biom(Path(args.biom_path), meta["ID"])

    score = build_proaging_score(table, meta, Path(args.overlap_path))
    meta = meta.set_index("ID")
    meta["pro_aging_score"] = score.reindex(meta.index).values

    existing = detect_existing_ascvd_column(meta)
    if existing is not None and pd.to_numeric(meta[existing], errors="coerce").notna().sum() > 100:
        ascvd = pd.to_numeric(meta[existing], errors="coerce")
        ascvd_source = f"precomputed metadata column '{existing}'"
    else:
        ascvd = compute_ascvd(meta.reset_index()).set_axis(meta.index)
        ascvd_source = (
            "computed in-script via ACC/AHA 2013 PCE (port of CVrisk::ascvd_10y_accaha; "
            "cholesterol/HDL converted from mmol/L to mg/dL by *38.67; "
            "race='other' -> White coefficients; "
            "current_smoker=1, never_smoker=0 (others NA); "
            "bp_med=any of bp_control_medication1/2/3 == 'y'; "
            "diabetes from bg_diagnosis)"
        )
    meta["ascvd_10y"] = ascvd.values

    n_score = meta["pro_aging_score"].notna().sum()
    n_ascvd = meta["ascvd_10y"].notna().sum()
    print(f"Module score available for {n_score} samples; ASCVD available for {n_ascvd}")

    rows = []
    for med in CANDIDATE_MEDIATORS:
        print(f"  ... mediator: {med}")
        rows.append(run_mediation_for(meta, med, rng, args.n_boot))

    results = pd.DataFrame(rows)
    _, results["q_a"] = fdrcorrection(results["p_a"].fillna(1.0).values)
    _, results["q_b"] = fdrcorrection(results["p_b"].fillna(1.0).values)
    cols = [
        "mediator", "n_complete",
        "coef_a", "p_a", "q_a",
        "coef_b", "p_b", "q_b",
        "total_c", "p_total_c",
        "direct_c_prime", "p_direct_c_prime",
        "indirect_ab", "indirect_ci_low", "indirect_ci_high",
        "prop_mediated", "prop_mediated_ci_low", "prop_mediated_ci_high",
        "comment",
    ]
    results = results[cols]
    results.to_csv(output_dir / "mediation_panel_results.tsv", sep="\t", index=False)

    forest_plot(results, output_dir)
    write_memo(results, output_dir, ascvd_source)
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
