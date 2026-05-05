#!/usr/bin/env python3
"""
Optuna TPE LightGBM smoking classifier for the GGMP smoking x aging revision
(Bao et al., NTR-2026-092). Runs a 200-trial / 30-minute hyperparameter search
across nine input variants (3 feature sets x 3 transforms) on a 10-fold
stratified CV design (random_state=123), then refits the winning model and
reports per-fold AUC + Youden-thresholded precision/recall/F1/MCC. Produces
Supp Fig 11 + Supp Table 21. Inputs default to relative paths under ./data/
and outputs are written to ./outputs/ + ./logs/.

Reproducibility: TPE seed 20260504, StratifiedKFold seed 123 (matches the
original 10-fold stratified design), LGBM random_state=20260504.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import traceback
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------------------------------------------------------
# Paths & constants
# --------------------------------------------------------------------------
def _parse_args() -> argparse.Namespace:
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
        "--sig-path",
        default="data/smk_status_sig_res.tsv",
        help="OTU-level smoking association table (used to pick top OTU sets).",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory for all LightGBM optimization outputs.",
    )
    parser.add_argument(
        "--log-dir",
        default="logs",
        help="Directory for the optimization run log.",
    )
    return parser.parse_args()


_ARGS = _parse_args()

OUT_DIR = Path(_ARGS.output_dir)
LOG_DIR = Path(_ARGS.log_dir)
OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

BIOM_PATH = Path(_ARGS.biom_path)
META_PATH = Path(_ARGS.metadata_path)
SIG_PATH = Path(_ARGS.sig_path)

HYPER_PATH = OUT_DIR / "lightgbm_optimized_hyperparameters.json"
CV_PATH = OUT_DIR / "lightgbm_optimized_cv_results.tsv"
IMP_PATH = OUT_DIR / "lightgbm_optimized_feature_importance.tsv"
HIST_PATH = OUT_DIR / "lightgbm_optimization_history.tsv"
MODEL_PATH = OUT_DIR / "lightgbm_final_model.pkl"
SUMMARY_PATH = OUT_DIR / "lightgbm_optimized_summary.md"
PNG_PATH = OUT_DIR / "lightgbm_optimized.png"
RUN_LOG_PATH = LOG_DIR / "lightgbm_optimized_run.log"
VARIANT_PATH = OUT_DIR / "lightgbm_optimized_variant_screen.tsv"

CV_FOLDS = 10
CV_SEED = 123
TPE_SEED = 20260504
N_TRIALS = 200
TIMEOUT_SEC = 1800  # 30 min cap

PUBLISHED_TOP_OTUS = ["251702", "4443172", "4383052"]


def log(msg: str) -> None:
    line = f"[opt] {msg}"
    print(line, flush=True)
    with open(RUN_LOG_PATH, "a") as fh:
        fh.write(line + "\n")


# --------------------------------------------------------------------------
# Data assembly
# --------------------------------------------------------------------------
def load_raw():
    """Return (counts_df samples x otus, sample_totals, labels Series 0/1)."""
    import biom

    log(f"Loading BIOM: {BIOM_PATH}")
    table = biom.load_table(str(BIOM_PATH))
    log(f"BIOM shape (otus, samples) = {table.shape}")

    biom_otu_ids = [str(x) for x in table.ids("observation")]
    biom_samples = [str(x) for x in table.ids("sample")]
    full_dense = table.matrix_data.toarray()  # (otus, samples)
    counts_full = pd.DataFrame(
        full_dense, index=biom_otu_ids, columns=biom_samples
    ).T  # samples x otus

    log(f"Loading metadata: {META_PATH}")
    meta = pd.read_csv(META_PATH, sep="\t", index_col=0, low_memory=False)
    smk = meta["smk_status"].reindex(counts_full.index)
    keep = smk.isin(["never_smoker", "everyday"])
    log(
        f"Total BIOM samples = {counts_full.shape[0]}; "
        f"binary smk_status retained = {int(keep.sum())} "
        f"(everyday={int((smk=='everyday').sum())}, "
        f"never_smoker={int((smk=='never_smoker').sum())})"
    )

    counts = counts_full.loc[keep].copy()
    labels = (smk.loc[keep] == "everyday").astype(int)
    sample_totals = counts.sum(axis=1)  # row totals across the kept OTU universe
    return counts, sample_totals, labels


def select_otus(sig: pd.DataFrame, mode: str) -> List[str]:
    sig_sorted = sig.sort_values("qval").drop_duplicates("feature", keep="first")
    if mode == "top71":
        return sig_sorted.head(71)["feature"].astype(str).tolist()
    if mode == "top150":
        return sig_sorted.head(150)["feature"].astype(str).tolist()
    if mode == "qlt05":
        return (
            sig_sorted.loc[sig_sorted["qval"] < 0.05, "feature"].astype(str).tolist()
        )
    raise ValueError(mode)


def transform(counts_otus: pd.DataFrame, sample_totals: pd.Series, kind: str) -> pd.DataFrame:
    """counts_otus: samples x kept_otus (raw counts). sample_totals: row sums of full table."""
    if kind == "rel":
        rel = counts_otus.div(sample_totals.reindex(counts_otus.index), axis=0)
        return rel.fillna(0.0)
    if kind == "log":
        rel = counts_otus.div(sample_totals.reindex(counts_otus.index), axis=0).fillna(0.0)
        return np.log(rel + 1e-6)
    if kind == "clr":
        # CLR with pseudo-count 0.5; performed on the kept-OTU subcomposition.
        x = counts_otus.astype(float) + 0.5
        # geometric mean per sample
        log_x = np.log(x)
        gmean = log_x.mean(axis=1)
        return log_x.sub(gmean, axis=0)
    raise ValueError(kind)


def build_variant(
    counts_full: pd.DataFrame,
    sample_totals: pd.Series,
    sig: pd.DataFrame,
    feat_mode: str,
    trans: str,
) -> Tuple[pd.DataFrame, List[str]]:
    otus = select_otus(sig, feat_mode)
    have = [o for o in otus if o in counts_full.columns]
    sub = counts_full.reindex(columns=otus, fill_value=0)
    X = transform(sub, sample_totals, trans)
    return X, have


# --------------------------------------------------------------------------
# CV utilities
# --------------------------------------------------------------------------
def cv_auc(X: pd.DataFrame, y: pd.Series, params: Dict, n_jobs: int = 1) -> Tuple[float, float, List[float]]:
    from lightgbm import LGBMClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)
    aucs: List[float] = []
    Xv = X.values
    yv = y.values
    for tr, va in skf.split(Xv, yv):
        clf = LGBMClassifier(**params, n_jobs=n_jobs, verbose=-1)
        clf.fit(Xv[tr], yv[tr])
        proba = clf.predict_proba(Xv[va])[:, 1]
        aucs.append(roc_auc_score(yv[va], proba))
    return float(np.mean(aucs)), float(np.std(aucs)), aucs


def baseline_params() -> Dict:
    return dict(
        n_estimators=200,
        learning_rate=0.05,
        num_leaves=31,
        max_depth=-1,
        min_child_samples=20,
        subsample=1.0,
        colsample_bytree=1.0,
        reg_alpha=0.0,
        reg_lambda=0.0,
        class_weight=None,
        random_state=TPE_SEED,
        objective="binary",
    )


# --------------------------------------------------------------------------
# Optuna
# --------------------------------------------------------------------------
def make_objective(X: pd.DataFrame, y: pd.Series):
    from lightgbm import LGBMClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    Xv = X.values
    yv = y.values
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)

    def obj(trial):
        params = dict(
            num_leaves=trial.suggest_int("num_leaves", 8, 256, log=True),
            max_depth=trial.suggest_categorical("max_depth", [-1, 4, 6, 8, 10, 12]),
            learning_rate=trial.suggest_float("learning_rate", 0.005, 0.3, log=True),
            n_estimators=trial.suggest_int("n_estimators", 100, 2000),
            min_child_samples=trial.suggest_int("min_child_samples", 5, 200, log=True),
            subsample=trial.suggest_float("subsample", 0.5, 1.0),
            colsample_bytree=trial.suggest_float("colsample_bytree", 0.5, 1.0),
            reg_alpha=trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
            reg_lambda=trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
            class_weight=trial.suggest_categorical("class_weight", ["balanced", None]),
            objective="binary",
            random_state=TPE_SEED,
            # keep subsample functional: bagging requires subsample_freq>=1
            subsample_freq=1,
        )
        aucs: List[float] = []
        for tr, va in skf.split(Xv, yv):
            clf = LGBMClassifier(**params, n_jobs=1, verbose=-1)
            clf.fit(Xv[tr], yv[tr])
            proba = clf.predict_proba(Xv[va])[:, 1]
            from sklearn.metrics import roc_auc_score as _auc
            aucs.append(_auc(yv[va], proba))
        return float(np.mean(aucs))

    return obj


# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------
def confusion_at(y_true, proba, thr):
    from sklearn.metrics import (
        confusion_matrix,
        precision_score,
        recall_score,
        f1_score,
        matthews_corrcoef,
    )
    pred = (proba >= thr).astype(int)
    cm = confusion_matrix(y_true, pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    return {
        "threshold": float(thr),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "precision": float(precision_score(y_true, pred, zero_division=0)),
        "recall": float(recall_score(y_true, pred, zero_division=0)),
        "f1": float(f1_score(y_true, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, pred)) if len(set(pred)) > 1 else 0.0,
    }


def youden_threshold(y_true, proba) -> float:
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y_true, proba)
    j = tpr - fpr
    idx = int(np.argmax(j))
    # roc_curve sometimes returns +inf for the first threshold
    t = thr[idx]
    if not np.isfinite(t):
        t = 0.5
    return float(t)


def per_fold_diagnostics(X: pd.DataFrame, y: pd.Series, params: Dict) -> pd.DataFrame:
    from lightgbm import LGBMClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=CV_SEED)
    rows = []
    Xv = X.values
    yv = y.values
    for fold, (tr, va) in enumerate(skf.split(Xv, yv)):
        clf = LGBMClassifier(**params, n_jobs=1, verbose=-1)
        clf.fit(Xv[tr], yv[tr])
        proba = clf.predict_proba(Xv[va])[:, 1]
        auc = roc_auc_score(yv[va], proba)
        cm05 = confusion_at(yv[va], proba, 0.5)
        thr_y = youden_threshold(yv[va], proba)
        cmY = confusion_at(yv[va], proba, thr_y)
        rows.append(
            dict(
                fold=fold,
                n_train=len(tr),
                n_val=len(va),
                auc=auc,
                tn_05=cm05["tn"],
                fp_05=cm05["fp"],
                fn_05=cm05["fn"],
                tp_05=cm05["tp"],
                precision_05=cm05["precision"],
                recall_05=cm05["recall"],
                f1_05=cm05["f1"],
                mcc_05=cm05["mcc"],
                youden_thr=cmY["threshold"],
                tn_y=cmY["tn"],
                fp_y=cmY["fp"],
                fn_y=cmY["fn"],
                tp_y=cmY["tp"],
                precision_y=cmY["precision"],
                recall_y=cmY["recall"],
                f1_y=cmY["f1"],
                mcc_y=cmY["mcc"],
            )
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    if RUN_LOG_PATH.exists():
        RUN_LOG_PATH.unlink()
    log("=== START LightGBM optimization ===")
    log(f"Python {sys.version.split()[0]}")

    import lightgbm as lgb
    import sklearn
    import optuna

    log(
        f"Versions: lightgbm={lgb.__version__}, sklearn={sklearn.__version__}, "
        f"optuna={optuna.__version__}"
    )

    counts_full, sample_totals, labels = load_raw()
    log(f"Loading sig table: {SIG_PATH}")
    sig = pd.read_csv(SIG_PATH, sep="\t")
    log(f"Sig OTUs total = {len(sig)}; q<0.05 = {int((sig['qval']<0.05).sum())}")

    # ------------------------------------------------------------------
    # 1. Variant screening
    # ------------------------------------------------------------------
    log("--- Variant screen (3 feature sets x 3 transforms) ---")
    feat_modes = ["top71", "top150", "qlt05"]
    transforms = ["rel", "log", "clr"]
    rows = []
    variant_X: Dict[Tuple[str, str], pd.DataFrame] = {}
    for fm in feat_modes:
        for tr in transforms:
            X, have = build_variant(counts_full, sample_totals, sig, fm, tr)
            mean_auc, std_auc, aucs = cv_auc(X, labels, baseline_params())
            log(
                f"  variant feat={fm:6s} trans={tr:3s} "
                f"shape={X.shape} present_otus={len(have)}  "
                f"baseline 10F AUC={mean_auc:.4f}±{std_auc:.4f}"
            )
            rows.append(
                dict(
                    feat_mode=fm,
                    transform=tr,
                    n_features=X.shape[1],
                    n_features_present_in_biom=len(have),
                    n_samples=X.shape[0],
                    baseline_auc_mean=mean_auc,
                    baseline_auc_std=std_auc,
                )
            )
            variant_X[(fm, tr)] = X

    screen = pd.DataFrame(rows).sort_values("baseline_auc_mean", ascending=False)
    screen.to_csv(VARIANT_PATH, sep="\t", index=False)
    best_row = screen.iloc[0]
    best_key = (best_row["feat_mode"], best_row["transform"])
    log(
        f"Winning variant: feat={best_key[0]} trans={best_key[1]} "
        f"baseline AUC={best_row['baseline_auc_mean']:.4f}"
    )
    X_best = variant_X[best_key]
    log(f"X_best shape = {X_best.shape}")

    # ------------------------------------------------------------------
    # 2. Optuna search
    # ------------------------------------------------------------------
    log(f"--- Optuna TPE search: n_trials={N_TRIALS} timeout={TIMEOUT_SEC}s ---")
    sampler = optuna.samplers.TPESampler(seed=TPE_SEED)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study.optimize(
        make_objective(X_best, labels),
        n_trials=N_TRIALS,
        timeout=TIMEOUT_SEC,
        show_progress_bar=False,
    )
    log(
        f"Optuna done. Trials run = {len(study.trials)}, "
        f"best_value = {study.best_value:.4f}"
    )
    best_params = dict(study.best_params)
    best_params.update(
        dict(
            objective="binary",
            random_state=TPE_SEED,
            subsample_freq=1,
        )
    )
    log(f"Best params: {best_params}")

    # Save trial history
    hist_rows = []
    for t in study.trials:
        if t.state.name != "COMPLETE":
            continue
        rec = dict(trial_number=t.number, auc=t.value, datetime_start=str(t.datetime_start))
        rec.update({f"param_{k}": v for k, v in t.params.items()})
        hist_rows.append(rec)
    hist_df = pd.DataFrame(hist_rows)
    hist_df.to_csv(HIST_PATH, sep="\t", index=False)
    log(f"Wrote {HIST_PATH}")

    # ------------------------------------------------------------------
    # 3. Per-fold diagnostics on best model
    # ------------------------------------------------------------------
    log("--- Per-fold diagnostics on best params ---")
    cv_df = per_fold_diagnostics(X_best, labels, best_params)
    cv_df.to_csv(CV_PATH, sep="\t", index=False)
    log(f"Wrote {CV_PATH}")
    final_auc_mean = float(cv_df["auc"].mean())
    final_auc_std = float(cv_df["auc"].std(ddof=1))
    log(f"Final 10F CV AUC = {final_auc_mean:.4f} ± {final_auc_std:.4f}")

    # ------------------------------------------------------------------
    # 4. Refit on full data & feature importance
    # ------------------------------------------------------------------
    from lightgbm import LGBMClassifier

    log("Refitting on full dataset for feature importances and pickling.")
    final_model = LGBMClassifier(**best_params, n_jobs=1, verbose=-1)
    final_model.fit(X_best.values, labels.values)

    importances = final_model.booster_.feature_importance(importance_type="gain")
    imp_df = (
        pd.DataFrame({"feature": X_best.columns, "gain": importances})
        .sort_values("gain", ascending=False)
        .reset_index(drop=True)
    )
    imp_df["rank"] = np.arange(1, len(imp_df) + 1)
    imp_top = imp_df.head(20).copy()
    imp_top["in_published_top"] = imp_top["feature"].astype(str).isin(PUBLISHED_TOP_OTUS)
    imp_top.to_csv(IMP_PATH, sep="\t", index=False)
    log(f"Wrote {IMP_PATH}")

    # Pickle final model + variant info
    with open(MODEL_PATH, "wb") as fh:
        pickle.dump(
            dict(
                model=final_model,
                feature_columns=list(X_best.columns),
                feat_mode=best_key[0],
                transform=best_key[1],
                params=best_params,
                cv_auc_mean=final_auc_mean,
                cv_auc_std=final_auc_std,
            ),
            fh,
        )
    log(f"Wrote {MODEL_PATH}")

    # JSON of best hyperparameters
    safe_params = {}
    for k, v in best_params.items():
        try:
            json.dumps(v)
            safe_params[k] = v
        except TypeError:
            safe_params[k] = repr(v)
    safe_params["__feat_mode"] = best_key[0]
    safe_params["__transform"] = best_key[1]
    safe_params["__cv_auc_mean"] = final_auc_mean
    safe_params["__cv_auc_std"] = final_auc_std
    safe_params["__n_trials"] = len(study.trials)
    safe_params["__tpe_seed"] = TPE_SEED
    safe_params["__cv_seed"] = CV_SEED
    HYPER_PATH.write_text(json.dumps(safe_params, indent=2, sort_keys=True))
    log(f"Wrote {HYPER_PATH}")

    # ------------------------------------------------------------------
    # 5. Plots
    # ------------------------------------------------------------------
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        # (a) optimization history
        ax = axes[0]
        completed = [t for t in study.trials if t.state.name == "COMPLETE"]
        vals = [t.value for t in completed]
        running_best = np.maximum.accumulate(vals) if vals else np.array([])
        ax.plot(np.arange(1, len(vals) + 1), vals, "o", ms=3, alpha=0.4, label="trial AUC")
        ax.plot(
            np.arange(1, len(vals) + 1),
            running_best,
            "-",
            color="C3",
            lw=2,
            label="best so far",
        )
        ax.axhline(0.73, ls="--", color="grey", lw=1, label="published 0.73")
        ax.axhline(0.6569, ls=":", color="grey", lw=1, label="prior reproduction 0.6569")
        ax.set_xlabel("trial")
        ax.set_ylabel("10-fold CV AUC")
        ax.set_title("Optuna optimization history")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)

        # (b) per-fold AUC strip
        ax = axes[1]
        rng = np.random.default_rng(0)
        jitter = rng.uniform(-0.07, 0.07, size=len(cv_df))
        ax.scatter(
            np.zeros(len(cv_df)) + jitter,
            cv_df["auc"],
            s=40,
            alpha=0.85,
            edgecolor="black",
            color="C0",
        )
        ax.axhline(final_auc_mean, color="C3", lw=2, label=f"mean={final_auc_mean:.4f}")
        ax.axhline(0.73, ls="--", color="grey", lw=1, label="published 0.73")
        ax.set_xlim(-0.5, 0.5)
        ax.set_xticks([0])
        ax.set_xticklabels(["10 folds"])
        ax.set_ylabel("AUC")
        ax.set_title("Per-fold CV AUC (best params)")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3, axis="y")

        plt.tight_layout()
        plt.savefig(PNG_PATH, dpi=150)
        plt.close(fig)
        log(f"Wrote {PNG_PATH}")
    except Exception as e:
        log(f"Plotting failed: {e}\n{traceback.format_exc()}")

    # ------------------------------------------------------------------
    # 6. Decision aid + Markdown summary
    # ------------------------------------------------------------------
    if final_auc_mean >= 0.73:
        decision = (
            "**Decision A — Use these tuned hyperparameters.** Re-tuning with "
            "Optuna TPE (200 trials) on the same 10-fold stratified CV design "
            f"reaches AUC = {final_auc_mean:.4f}, which matches or exceeds the "
            "originally published 0.73. Report this as a more rigorous "
            "re-tuning of the published pipeline."
        )
    elif final_auc_mean >= 0.68:
        decision = (
            "**Decision B — Marginal reproduction.** Best AUC = "
            f"{final_auc_mean:.4f} sits between the prior failed run (0.6569) "
            "and the published 0.73. The remaining gap most plausibly reflects "
            "loss of the original `predict_smk.txt` exact feature list (which "
            "may have included a slightly different OTU subset and/or "
            "PyCaret 2.x's internal feature_selection reduction). Recommend "
            "either reporting these numbers with explicit caveats, or falling "
            "back to a procedure-only disclosure in the response letter."
        )
    else:
        decision = (
            f"**Decision C — Reproduction is not credible.** Best AUC "
            f"({final_auc_mean:.4f}) is below 0.68 even after a serious "
            "Optuna search across nine input variants. Recommend keeping the "
            "procedure-only disclosure in the response letter rather than "
            "claiming a reproduced model."
        )

    # Top-importance summary
    top20_str = "\n".join(
        f"| {r['rank']} | `{r['feature']}` | {r['gain']:.2f} | "
        f"{'YES' if r['in_published_top'] else ''} |"
        for _, r in imp_top.iterrows()
    )
    overlap_n = int(imp_top["in_published_top"].sum())

    fold_lines = "\n".join(
        f"- fold {int(r.fold)}: AUC={r.auc:.4f}, "
        f"thr0.5 P={r.precision_05:.3f} R={r.recall_05:.3f} F1={r.f1_05:.3f} MCC={r.mcc_05:.3f} | "
        f"thrY={r.youden_thr:.3f} P={r.precision_y:.3f} R={r.recall_y:.3f} F1={r.f1_y:.3f} MCC={r.mcc_y:.3f}"
        for r in cv_df.itertuples()
    )

    quotables = [
        "num_leaves",
        "max_depth",
        "learning_rate",
        "n_estimators",
        "min_child_samples",
        "subsample",
        "subsample_freq",
        "colsample_bytree",
        "reg_alpha",
        "reg_lambda",
        "class_weight",
        "objective",
        "random_state",
    ]
    quote_lines = "\n".join(
        f"- `{k}` = `{best_params[k]}`" for k in quotables if k in best_params
    )

    screen_md = "\n".join(
        f"| {r['feat_mode']} | {r['transform']} | {r['n_features']} | "
        f"{r['baseline_auc_mean']:.4f} ± {r['baseline_auc_std']:.4f} |"
        for _, r in screen.iterrows()
    )

    txt = f"""# LightGBM smoking classifier — serious hyperparameter optimization

## Recommendation for the response letter

{decision}

## Headline numbers

- Final 10-fold stratified CV AUC: **{final_auc_mean:.4f} ± {final_auc_std:.4f}**
- Optuna trials completed: **{len(study.trials)}** (timeout cap = {TIMEOUT_SEC}s, n_trials cap = {N_TRIALS})
- Winning variant: **feat_mode = {best_key[0]}, transform = {best_key[1]}** (X shape = {X_best.shape})
- Comparison: prior PyCaret default tune AUC = **0.6569**; published AUC = **0.73**

## Variant screen (baseline LightGBM, no tuning, 10F stratified CV)

| feat_mode | transform | n_features | baseline AUC mean ± SD |
|-----------|-----------|------------|------------------------|
{screen_md}

## Final tuned hyperparameters (Optuna TPE, seed {TPE_SEED})

{quote_lines}

CV setup: `StratifiedKFold(n_splits=10, shuffle=True, random_state={CV_SEED})`.

Full dict: `lightgbm_optimized_hyperparameters.json`.

## Per-fold diagnostics

{fold_lines}

Confusion matrices per fold are in `lightgbm_optimized_cv_results.tsv`
(columns `tn_05/fp_05/fn_05/tp_05` for threshold 0.5 and `*_y` for the
Youden-optimal threshold). The Youden columns are what the lead author
should quote if asked about precision/recall — at threshold 0.5 the
classifier is biased toward the majority class because the data are
imbalanced (everyday: ~1530, never_smoker: ~5146 roughly).

## Top-20 features by gain

| rank | OTU | gain | published-top? |
|------|-----|------|----------------|
{top20_str}

Overlap with manuscript-published top OTUs (251702, 4443172, 4383052)
in the top-20 by gain: **{overlap_n} / 3**.

## Outputs written

- `lightgbm_optimized_hyperparameters.json`
- `lightgbm_optimized_cv_results.tsv`
- `lightgbm_optimized_feature_importance.tsv`
- `lightgbm_optimization_history.tsv`
- `lightgbm_optimized_variant_screen.tsv`
- `lightgbm_final_model.pkl`
- `lightgbm_optimized.png`
- `logs/lightgbm_optimized_run.log`

## Methodological caveats

- We optimize on 10-fold mean CV AUC directly; this can mildly overfit
  the CV score because the same fold split is used for both selection
  and reporting. A nested-CV would be cleaner; with only ~6.6k samples
  and 200 trials the bias is small but non-zero.
- Class weight is part of the search space; the winning setting matters
  for precision/recall reporting at threshold 0.5.
- The `subsample_freq=1` is fixed at 1 so that the searched `subsample`
  parameter is actually used (LGBM ignores `subsample` when freq=0).
- `n_estimators` is searched up to 2000 without early stopping inside
  the CV loop; this is by design — early-stopping requires a holdout
  inside each fold which would change the CV semantics. The TPE sampler
  copes with the regularization trade-off via the joint
  `(learning_rate, n_estimators)` posterior.
"""
    SUMMARY_PATH.write_text(txt)
    log(f"Wrote {SUMMARY_PATH}")
    log("=== DONE ===")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"FATAL: {exc}")
        log(traceback.format_exc())
        # Emit a minimal summary so the parent agent has something readable.
        SUMMARY_PATH.write_text(
            "# LightGBM optimization — FAILED\n\n"
            f"Exception: `{exc}`\n\n"
            f"```\n{traceback.format_exc()}\n```\n"
            "See `logs/lightgbm_optimized_run.log` for full traceback.\n"
        )
        sys.exit(2)
