#!/usr/bin/env python3
"""
09_soce_risk_score.py  –  SOCE Composite Risk Score (TCGA-LIHC)
================================================================
Computes a composite SOCE risk score by combining the expression of
STIM1, TRPC6, TRPC1, and ORAI1 using PCA-based weights, then performs
Kaplan-Meier survival analysis on high vs. low risk patients.

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal (lihc_tcga)
  Clinical data              : cBioPortal

Method
------
  1. Log2-transform RSEM values
  2. Z-score normalise each gene across patients
  3. Compute PC1 score as composite risk score
  4. Dichotomise at median (high vs. low)
  5. KM survival analysis + log-rank test

Outputs
-------
  soce_risk_score_patients.csv     – Per-patient risk scores and group
  soce_risk_score_summary.csv      – Group summary statistics
  soce_risk_score_km.png / .svg
  soce_risk_score_km_v2.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn lifelines scikit-learn requests
  python soce_risk_score.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

OUTDIR     = os.path.dirname(os.path.abspath(__file__))
SOCE_GENES = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]


def fetch_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load expression and clinical data (uses sibling caches if available)."""
    # Expression
    expr_cache = None
    for p in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
    ):
        if os.path.exists(p):
            df = pd.read_csv(p, index_col=0)
            if all(g in df.columns for g in SOCE_GENES):
                expr_cache = df[SOCE_GENES]
                print(f"  [expr cache] {p}")
                break

    if expr_cache is None:
        print("  Downloading expression from cBioPortal …")
        base    = "https://www.cbioportal.org/api"
        profile = "lihc_tcga_rna_seq_v2_mrna"
        r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
        r.raise_for_status()
        sample_ids = r.json()
        r = requests.get(f"{base}/genes?geneIds={','.join(SOCE_GENES)}", timeout=60)
        r.raise_for_status()
        gene_map = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
        data = {}
        for symbol, entrez in gene_map.items():
            payload = {"entrezGeneId": entrez, "molecularProfileId": profile, "sampleIds": sample_ids}
            r = requests.post(f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
                              json=payload, timeout=120)
            for row in r.json():
                s = row["sampleId"]
                if s not in data: data[s] = {}
                data[s][symbol] = row["value"]
        expr_cache = pd.DataFrame.from_dict(data, orient="index")[SOCE_GENES]

    # Clinical
    clin_cache = None
    for p in (
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_clinical_cache.csv"),
        os.path.join(OUTDIR, "..", "07_cox_regression", "_clinical_cache.csv"),
    ):
        if os.path.exists(p):
            clin_cache = pd.read_csv(p, index_col=0)
            print(f"  [clin cache] {p}")
            break

    if clin_cache is None:
        print("  Downloading clinical data …")
        base = "https://www.cbioportal.org/api"
        r = requests.get(
            f"{base}/studies/lihc_tcga/clinical-data?clinicalDataType=SAMPLE",
            timeout=120)
        r.raise_for_status()
        clin_cache = pd.DataFrame(r.json()).pivot(
            index="sampleId", columns="clinicalAttributeId", values="value")

    return expr_cache, clin_cache


def compute_risk_score(expr: pd.DataFrame) -> pd.DataFrame:
    """Z-score → PCA(1 component) → risk score per patient."""
    X = expr.dropna()
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=1)
    scores = pca.fit_transform(X_scaled).flatten()
    loadings = dict(zip(SOCE_GENES, pca.components_[0]))
    explained = pca.explained_variance_ratio_[0]

    print(f"  PC1 explained variance: {explained:.1%}")
    print(f"  Loadings: { {k: round(v,3) for k,v in loadings.items()} }")

    df_score = pd.DataFrame({
        "risk_score": scores,
        "group":      pd.cut(scores, bins=[-np.inf, np.median(scores), np.inf],
                             labels=["Low", "High"]),
    }, index=X.index)
    df_score = df_score.join(X)
    return df_score, loadings, explained


def km_analysis(df_score: pd.DataFrame, clin: pd.DataFrame):
    """Kaplan-Meier survival for High vs. Low risk groups."""
    os_time_col  = next((c for c in clin.columns if "OS_MONTHS"  in c.upper()), None)
    os_event_col = next((c for c in clin.columns if "OS_STATUS"  in c.upper()), None)
    if os_time_col is None:
        print("  ⚠  No OS columns found. Skipping KM."); return

    merged = df_score.join(clin[[os_time_col, os_event_col]], how="inner")
    merged["os_time"]  = pd.to_numeric(merged[os_time_col],  errors="coerce")
    merged["os_event"] = (
        merged[os_event_col].astype(str).str.extract(r"(\d+)")[0]
        .astype(float).fillna(0).astype(int)
    )
    merged = merged.dropna(subset=["os_time", "os_event"])

    grp_h = merged[merged["group"] == "High"]
    grp_l = merged[merged["group"] == "Low"]

    lr = logrank_test(grp_h["os_time"], grp_l["os_time"],
                      grp_h["os_event"], grp_l["os_event"])

    # Plot KM
    for version, palette in [("km", ["#c0392b", "#2980b9"]),
                              ("km_v2", ["#8e44ad", "#16a085"])]:
        fig, ax = plt.subplots(figsize=(9, 6))
        kmf_h = KaplanMeierFitter()
        kmf_l = KaplanMeierFitter()
        kmf_h.fit(grp_h["os_time"], grp_h["os_event"],
                  label=f"High Risk (n={len(grp_h)})")
        kmf_l.fit(grp_l["os_time"], grp_l["os_event"],
                  label=f"Low Risk  (n={len(grp_l)})")
        kmf_h.plot_survival_function(ax=ax, ci_show=True, color=palette[0])
        kmf_l.plot_survival_function(ax=ax, ci_show=True, color=palette[1])

        p_str = f"p = {lr.p_value:.4f}" if lr.p_value >= 0.001 else f"p = {lr.p_value:.2e}"
        ax.text(0.98, 0.98, p_str, transform=ax.transAxes,
                ha="right", va="top", fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.9))
        ax.set_xlabel("Months", fontsize=12)
        ax.set_ylabel("Overall Survival Probability", fontsize=12)
        ax.set_title(
            "SOCE Composite Risk Score — Kaplan-Meier Survival\n"
            "(TCGA-LIHC, PC1 of STIM1 + TRPC6 + TRPC1 + ORAI1)",
            fontsize=12, fontweight="bold"
        )
        plt.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(OUTDIR, f"soce_risk_score_{version}.{ext}"),
                        dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  Saved soce_risk_score_{version}.png / .svg")

    print(f"  Log-rank p = {lr.p_value:.4f}")
    return lr.p_value


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 09 — SOCE Composite Risk Score")
    print("=" * 60)

    print("\n[1/3] Loading data …")
    expr, clin = fetch_data()

    print("\n[2/3] Computing PCA-based risk score …")
    df_score, loadings, var_exp = compute_risk_score(expr)

    # Save patient-level data
    df_score.to_csv(os.path.join(OUTDIR, "soce_risk_score_patients.csv"))
    print("  Saved soce_risk_score_patients.csv")

    summary = df_score.groupby("group")["risk_score"].agg(["count", "mean", "std"])
    summary.index.name = "group"
    summary.to_csv(os.path.join(OUTDIR, "soce_risk_score_summary.csv"))
    print("  Saved soce_risk_score_summary.csv")

    print("\n[3/3] KM survival analysis …")
    km_analysis(df_score, clin)

    print("\n  Done. All outputs saved to:", OUTDIR)
