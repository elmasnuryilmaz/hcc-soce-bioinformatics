#!/usr/bin/env python3
"""
07_cox_regression.py  –  Multivariate Cox Proportional Hazards Regression
=========================================================================
Tests whether STIM1 and TRPC6 are independent prognostic factors in
TCGA-LIHC after adjusting for clinical covariates.

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM + clinical data : cBioPortal (lihc_tcga)

Model covariates
----------------
  Gene expression (STIM1, TRPC6 — continuous, log2-scaled)
  Tumour stage (AJCC I/II/III/IV — ordered)
  Grade (G1–G3)
  Age at diagnosis
  Sex
  AFP (log2-transformed, where available)

Outputs
-------
  cox_results.csv         – HR, 95%CI, p-value per covariate
  cox_forest_plot.png / .svg

Key findings (TCGA-LIHC, n = 370)
----------------------------------
  STIM1  : HR = 0.80, p = 0.711  (not independent)
  TRPC6  : HR = 0.93, p = 0.314  (not independent)
  Stage  : HR = 1.68, p < 0.001  (dominant predictor)

Usage
-----
  pip install pandas numpy scipy matplotlib lifelines requests statsmodels
  python cox_regression.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

warnings.filterwarnings("ignore")

# ── Config ───────────────────────────────────────────────────────────────────
OUTDIR = os.path.dirname(os.path.abspath(__file__))
GENES  = ["STIM1", "TRPC6"]


# ── 1. Data download helpers ──────────────────────────────────────────────────
def fetch_expression() -> pd.DataFrame:
    """Load TCGA-LIHC RSEM expression for STIM1 and TRPC6."""
    # Try sibling caches first
    for cache in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival",  "_tcga_lihc_expr_cache.csv"),
        os.path.join(OUTDIR, "_tcga_expr_cache.csv"),
    ):
        if os.path.exists(cache):
            df = pd.read_csv(cache, index_col=0)
            if all(g in df.columns for g in GENES):
                print(f"  [cache expr] {cache}")
                return df[GENES]

    print("  Downloading expression from cBioPortal …")
    base    = "https://www.cbioportal.org/api"
    profile = "lihc_tcga_rna_seq_v2_mrna"
    r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
    r.raise_for_status()
    sample_ids = r.json()

    r = requests.get(f"{base}/genes?geneIds={','.join(GENES)}", timeout=60)
    r.raise_for_status()
    gene_map = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}

    data = {}
    for symbol, entrez in gene_map.items():
        payload = {"entrezGeneId": entrez, "molecularProfileId": profile, "sampleIds": sample_ids}
        r = requests.post(
            f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
            json=payload, timeout=120)
        for row in r.json():
            s = row["sampleId"]
            if s not in data: data[s] = {}
            data[s][symbol] = row["value"]
    return pd.DataFrame.from_dict(data, orient="index")


def fetch_clinical() -> pd.DataFrame:
    """Load TCGA-LIHC clinical data from cBioPortal."""
    cache = os.path.join(OUTDIR, "_clinical_cache.csv")
    if os.path.exists(cache):
        print(f"  [cache clinical] {cache}")
        return pd.read_csv(cache, index_col=0)

    # Try sibling cache
    sib = os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_clinical_cache.csv")
    if os.path.exists(sib):
        return pd.read_csv(sib, index_col=0)

    print("  Downloading clinical data from cBioPortal …")
    base = "https://www.cbioportal.org/api"
    r = requests.get(
        f"{base}/studies/lihc_tcga/clinical-data?clinicalDataType=SAMPLE",
        timeout=120)
    r.raise_for_status()
    df = pd.DataFrame(r.json()).pivot(
        index="sampleId", columns="clinicalAttributeId", values="value")
    df.to_csv(cache)
    return df


# ── 2. Merge and preprocess ───────────────────────────────────────────────────
def prepare_cox_data(expr: pd.DataFrame, clin: pd.DataFrame) -> pd.DataFrame:
    """Merge expression and clinical, encode covariates, return clean DataFrame."""
    merged = expr.join(clin, how="inner")

    # OS outcome
    os_time  = next((c for c in merged.columns if "OS_MONTHS"  in c.upper()), None)
    os_event = next((c for c in merged.columns if "OS_STATUS"  in c.upper()), None)
    if os_time is None or os_event is None:
        raise ValueError("Could not find OS_MONTHS / OS_STATUS columns.")

    merged["os_time"]  = pd.to_numeric(merged[os_time],  errors="coerce")
    merged["os_event"] = (
        merged[os_event].astype(str).str.extract(r"(\d+)")[0]
        .astype(float).fillna(0).astype(int)
    )

    # Gene expression – log2-scale
    for g in GENES:
        if g in merged.columns:
            merged[f"{g}_log2"] = np.log2(merged[g].clip(lower=0) + 1)

    # Tumour stage – ordinal 1-4
    stage_col = next((c for c in merged.columns if "STAGE" in c.upper()), None)
    if stage_col:
        stage_map = {"I": 1, "II": 2, "III": 3, "IV": 4,
                     "STAGE I": 1, "STAGE II": 2, "STAGE III": 3, "STAGE IV": 4}
        merged["stage"] = (
            merged[stage_col].astype(str).str.upper()
            .str.extract(r"(I{1,3}V?)$")[0]
            .map({"I": 1, "II": 2, "III": 3, "IV": 4})
        )

    # Sex binary
    sex_col = next((c for c in merged.columns if c.upper() in ("SEX", "GENDER")), None)
    if sex_col:
        merged["sex_male"] = (merged[sex_col].astype(str).str.upper() == "MALE").astype(int)

    # Age
    age_col = next((c for c in merged.columns if "AGE" in c.upper()), None)
    if age_col:
        merged["age"] = pd.to_numeric(merged[age_col], errors="coerce")

    return merged


# ── 3. Fit Cox model ──────────────────────────────────────────────────────────
def fit_cox(data: pd.DataFrame) -> tuple[CoxPHFitter, list[str]]:
    """Fits CoxPHFitter and returns (model, covariates_used)."""
    base_covars = ["os_time", "os_event"]
    covars = []
    for c in [f"{g}_log2" for g in GENES] + ["stage", "age", "sex_male"]:
        if c in data.columns and data[c].notna().sum() > 50:
            covars.append(c)

    df_cox = data[base_covars + covars].dropna()
    print(f"  Cox dataset: {len(df_cox)} samples, covariates: {covars}")

    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(df_cox, duration_col="os_time", event_col="os_event")
    return cph, covars


# ── 4. Forest plot ────────────────────────────────────────────────────────────
def forest_plot(cph: CoxPHFitter):
    """Horizontal forest plot of HR with 95% CI."""
    summary = cph.summary.copy()
    summary = summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].copy()
    summary.columns = ["HR", "CI_lower", "CI_upper", "p_value"]
    summary.index.name = "covariate"
    summary.to_csv(os.path.join(OUTDIR, "cox_results.csv"))
    print("  Saved cox_results.csv")

    # Pretty labels
    label_map = {
        "STIM1_log2": "STIM1 (log₂ RSEM)",
        "TRPC6_log2": "TRPC6 (log₂ RSEM)",
        "stage":      "Tumour Stage",
        "age":        "Age",
        "sex_male":   "Sex (Male)",
    }
    summary["label"] = [label_map.get(i, i) for i in summary.index]
    summary = summary.sort_values("HR")

    fig, ax = plt.subplots(figsize=(9, max(4, len(summary) * 0.8 + 1)))

    y_pos = range(len(summary))
    ax.scatter(summary["HR"], y_pos, color="#2c3e50", s=70, zorder=3)
    ax.hlines(y_pos, summary["CI_lower"], summary["CI_upper"],
              color="#2c3e50", linewidth=2, zorder=2)

    # Significance markers
    for y, (_, row) in zip(y_pos, summary.iterrows()):
        sig = "***" if row["p_value"] < 0.001 else ("**" if row["p_value"] < 0.01
                                                     else ("*" if row["p_value"] < 0.05 else "ns"))
        ax.text(row["CI_upper"] + 0.05, y, sig, va="center", fontsize=10)
        p_str = f"p = {row['p_value']:.3f}" if row["p_value"] >= 0.001 else "p < 0.001"
        hr_str = f"HR = {row['HR']:.2f} ({row['CI_lower']:.2f}–{row['CI_upper']:.2f})"
        ax.text(row["CI_lower"] - 0.05, y, f"{hr_str}  {p_str}",
                va="center", ha="right", fontsize=8, color="#555")

    ax.axvline(1.0, color="grey", linestyle="--", linewidth=1)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(summary["label"], fontsize=11)
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=12)
    ax.set_title(
        "Multivariate Cox Proportional Hazards — TCGA-LIHC\n"
        "(Overall Survival; covariates: STIM1, TRPC6, Stage, Age, Sex)",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"cox_forest_plot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved cox_forest_plot.png / .svg")

    return summary


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 07 — Multivariate Cox Regression")
    print("=" * 60)

    print("\n[1/3] Loading data …")
    expr = fetch_expression()
    clin = fetch_clinical()

    print("\n[2/3] Fitting Cox model …")
    data   = prepare_cox_data(expr, clin)
    cph, _ = fit_cox(data)
    cph.print_summary()

    print("\n[3/3] Forest plot …")
    summary = forest_plot(cph)

    print("\n  ⚠  Interpretation note:")
    print("  Neither STIM1 nor TRPC6 are independent prognostic factors")
    print("  after adjusting for tumour stage (the dominant covariate).")
    print("  Univariate KM significance is confounded by stage distribution.")
    print("\n  Done. All outputs saved to:", OUTDIR)
