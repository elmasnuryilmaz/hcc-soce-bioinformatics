#!/usr/bin/env python3
"""
01_tcga_survival.py  –  Tumour vs. Normal Expression & Kaplan-Meier Survival
=============================================================================
Analysis of SOCE components (STIM1, TRPC6, TRPC1, ORAI1) in TCGA-LIHC.

Data sources
------------
  TCGA-LIHC RNA-seq v2 RSEM  : cBioPortal  (study: lihc_tcga)
  GTEx v8 Liver (normal)     : GTEx Portal (median TPM per gene)
  TCGA-LIHC clinical data    : cBioPortal

Outputs (saved to this folder)
-------------------------------
  tumor_vs_normal_stats.csv   – Mann-Whitney U statistics + log2FC
  tumor_vs_normal_boxplot.png / .svg
  km_survival_results.csv     – Log-rank p-values and HR per gene
  km_survival_curves.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn lifelines requests
  python tcga_survival.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, io, gzip, warnings
import sys
from pathlib import Path
import requests
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import median_survival_times

warnings.filterwarnings("ignore")

sys.path.append(str(Path(__file__).resolve().parents[1]))
from common.cbioportal import fetch_clinical_data, fetch_mrna_expression

# ── Config ─────────────────────────────────────────────────────────────────
OUTDIR   = os.path.dirname(os.path.abspath(__file__))
GENES    = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]
PALETTE  = {"Tumour": "#c0392b", "Normal": "#2980b9"}

# ── 1. Download TCGA-LIHC expression via cBioPortal API ────────────────────
def fetch_cbioportal_expression(genes: list[str]) -> pd.DataFrame:
    """
    Returns a DataFrame (samples × genes) of RNA-seq v2 RSEM values
    for the TCGA-LIHC cohort (n ≈ 370 tumours).
    """
    cache = os.path.join(OUTDIR, "_tcga_lihc_expr_cache.csv")
    if os.path.exists(cache):
        print(f"  [cache] Loading TCGA expression from {cache}")
        return pd.read_csv(cache, index_col=0)

    print("  Downloading TCGA-LIHC expression from cBioPortal …")
    df = fetch_mrna_expression(genes)
    df.to_csv(cache)
    print(f"  Saved cache: {cache}  ({df.shape[0]} samples)")
    return df


# ── 2. Download GTEx v8 liver median TPM ───────────────────────────────────
def fetch_gtex_liver(genes: list[str]) -> dict:
    """
    Returns {gene: median_TPM} for GTEx v8 liver tissue.
    Downloads the gene-median file (~50 MB, cached locally).
    """
    cache = os.path.join(OUTDIR, "_gtex_liver_cache.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache, index_col=0)
        return df["median_tpm"].to_dict()

    print("  Downloading GTEx v8 median TPM table …")
    url = (
        "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
        "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    )
    r = requests.get(url, timeout=300, stream=True)
    r.raise_for_status()
    content = b"".join(r.iter_content(chunk_size=1 << 20))
    with gzip.open(io.BytesIO(content)) as fh:
        raw = fh.read().decode()

    lines = raw.splitlines()
    # GCT format: skip 2 header lines, then column header
    col_line = lines[2].split("\t")
    liver_col = [i for i, c in enumerate(col_line) if "Liver" in c]
    if not liver_col:
        raise ValueError("Cannot find 'Liver' column in GTEx GCT file")

    records = []
    for line in lines[3:]:
        parts = line.split("\t")
        gene_sym = parts[1]
        vals = [float(parts[i]) for i in liver_col]
        records.append({"gene": gene_sym, "median_tpm": float(np.median(vals))})

    df_gtex = pd.DataFrame(records).set_index("gene")
    subset = df_gtex.loc[df_gtex.index.isin(genes)]
    subset.to_csv(cache)
    print(f"  GTEx liver cached: {cache}")
    return subset["median_tpm"].to_dict()


# ── 3. Download TCGA-LIHC clinical data ────────────────────────────────────
def fetch_tcga_clinical() -> pd.DataFrame:
    """Returns TCGA-LIHC clinical DataFrame from cBioPortal."""
    cache = os.path.join(OUTDIR, "_tcga_lihc_clinical_cache.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache, index_col=0)
        if {"OS_MONTHS", "OS_STATUS"}.issubset(df.columns):
            return df
        print(f"  [cache] Ignoring clinical cache without survival columns: {cache}")

    print("  Downloading TCGA-LIHC clinical data …")
    df = fetch_clinical_data()
    df.to_csv(cache)
    return df


# ── 4. Tumour vs. Normal analysis ──────────────────────────────────────────
def tumor_vs_normal(expr: pd.DataFrame, gtex_medians: dict) -> pd.DataFrame:
    """Mann-Whitney U + log2FC for each gene."""
    results = []
    for gene in GENES:
        if gene not in expr.columns:
            continue
        tumour_vals = expr[gene].dropna().values
        normal_med  = gtex_medians.get(gene, np.nan)
        if np.isnan(normal_med) or normal_med <= 0:
            continue

        # Use the GTEx median as a single reference point for log2FC
        tumour_med = float(np.median(tumour_vals))
        log2fc = np.log2((tumour_med + 1) / (normal_med + 1))

        # Mann-Whitney against a normal distribution centred on GTEx median
        # (single-sample Wilcoxon sign-rank against reference)
        stat, pval = stats.wilcoxon(tumour_vals - normal_med, alternative="greater")

        results.append(
            {
                "gene":          gene,
                "median_normal": round(normal_med, 3),
                "median_tumour": round(tumour_med, 3),
                "log2FC":        round(log2fc, 3),
                "stat":          round(stat, 2),
                "p_value":       pval,
            }
        )

    df_res = pd.DataFrame(results).set_index("gene")
    df_res.to_csv(os.path.join(OUTDIR, "tumor_vs_normal_stats.csv"))
    print(f"  Saved tumor_vs_normal_stats.csv")
    return df_res


def plot_boxplot(expr: pd.DataFrame, gtex_medians: dict):
    """Box plot: tumour (TCGA) vs. normal (GTEx median shown as line)."""
    fig, axes = plt.subplots(1, len(GENES), figsize=(14, 5), sharey=False)
    fig.suptitle(
        "SOCE Components: Tumour vs. Normal Expression\n(TCGA-LIHC vs. GTEx Liver v8)",
        fontsize=13, fontweight="bold",
    )
    for ax, gene in zip(axes, GENES):
        vals = expr[gene].dropna().values if gene in expr.columns else np.array([])
        ax.boxplot(vals, widths=0.5, patch_artist=True,
                   boxprops=dict(facecolor=PALETTE["Tumour"], alpha=0.7),
                   medianprops=dict(color="white", linewidth=2))
        gtex_m = gtex_medians.get(gene, None)
        if gtex_m is not None:
            ax.axhline(gtex_m, color=PALETTE["Normal"], linewidth=2,
                       linestyle="--", label=f"GTEx median: {gtex_m:.1f}")
            ax.legend(fontsize=8)
        ax.set_title(gene, fontweight="bold")
        ax.set_ylabel("RNA-seq v2 RSEM" if gene == GENES[0] else "")
        ax.set_xticks([1]); ax.set_xticklabels(["TCGA-LIHC\nTumour"])
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"tumor_vs_normal_boxplot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved tumor_vs_normal_boxplot.png / .svg")


# ── 5. Kaplan-Meier survival ────────────────────────────────────────────────
def km_survival(expr: pd.DataFrame, clinical: pd.DataFrame):
    """KM curves for each gene (high vs. low, split at median)."""
    # Align expression and clinical
    clin = clinical.copy()
    clin.index.name = "sample"

    # Extract OS columns
    os_time_col  = next((c for c in clin.columns if "OS_MONTHS" in c.upper()), None)
    os_event_col = next((c for c in clin.columns if "OS_STATUS" in c.upper()), None)

    if os_time_col is None or os_event_col is None:
        print("  ⚠  Could not find OS columns in clinical data. Skipping KM.")
        return

    merged = expr.join(clin[[os_time_col, os_event_col]], how="inner")
    merged[os_time_col]  = pd.to_numeric(merged[os_time_col],  errors="coerce")
    # OS_STATUS typically "0:LIVING" / "1:DECEASED"
    merged["event"] = (
        merged[os_event_col]
        .astype(str)
        .str.extract(r"(\d+)")[0]
        .astype(float)
        .fillna(0)
        .astype(int)
    )
    merged = merged.dropna(subset=[os_time_col, "event"])

    results = []
    n_genes = len(GENES)
    fig, axes = plt.subplots(1, n_genes, figsize=(5 * n_genes, 5))
    if n_genes == 1:
        axes = [axes]
    fig.suptitle(
        "Kaplan-Meier Overall Survival — SOCE Components\n(TCGA-LIHC, n≈370)",
        fontsize=13, fontweight="bold",
    )

    for ax, gene in zip(axes, GENES):
        if gene not in merged.columns:
            ax.set_visible(False); continue

        median_expr = merged[gene].median()
        grp_high = merged[merged[gene] >= median_expr]
        grp_low  = merged[merged[gene] <  median_expr]

        kmf_h = KaplanMeierFitter()
        kmf_l = KaplanMeierFitter()
        kmf_h.fit(grp_high[os_time_col], grp_high["event"], label=f"{gene}-High (n={len(grp_high)})")
        kmf_l.fit(grp_low[os_time_col],  grp_low["event"],  label=f"{gene}-Low  (n={len(grp_low)})")

        lr = logrank_test(
            grp_high[os_time_col], grp_low[os_time_col],
            grp_high["event"],     grp_low["event"],
        )
        pval = lr.p_value

        kmf_h.plot_survival_function(ax=ax, ci_show=True, color="#c0392b")
        kmf_l.plot_survival_function(ax=ax, ci_show=True, color="#2980b9")

        med_h = kmf_h.median_survival_time_
        med_l = kmf_l.median_survival_time_

        p_label = f"p = {pval:.3f}" if pval >= 0.001 else f"p = {pval:.2e}"
        ax.text(0.98, 0.98, p_label, transform=ax.transAxes,
                ha="right", va="top", fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        ax.set_title(gene, fontweight="bold")
        ax.set_xlabel("Months")
        ax.set_ylabel("Survival probability" if gene == GENES[0] else "")

        results.append({
            "gene":             gene,
            "n_high":           len(grp_high),
            "n_low":            len(grp_low),
            "median_OS_high":   float(med_h) if pd.notna(med_h) else np.nan,
            "median_OS_low":    float(med_l) if pd.notna(med_l) else np.nan,
            "logrank_p":        round(pval, 4),
        })

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"km_survival_curves.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved km_survival_curves.png / .svg")

    df_res = pd.DataFrame(results).set_index("gene")
    df_res.to_csv(os.path.join(OUTDIR, "km_survival_results.csv"))
    print("  Saved km_survival_results.csv")


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 01 — TCGA Survival Analysis")
    print("=" * 60)

    expr     = fetch_cbioportal_expression(GENES)
    gtex_med = fetch_gtex_liver(GENES)
    clinical = fetch_tcga_clinical()

    print("\n[1/3] Tumour vs. Normal …")
    tvn_df = tumor_vs_normal(expr, gtex_med)
    plot_boxplot(expr, gtex_med)
    print(tvn_df.to_string())

    print("\n[2/3] Kaplan-Meier Survival …")
    km_survival(expr, clinical)

    print("\n[3/3] Done. All outputs saved to:", OUTDIR)
