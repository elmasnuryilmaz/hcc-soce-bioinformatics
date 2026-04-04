#!/usr/bin/env python3
"""
05_coexpression.py  –  Spearman Co-expression Analysis (TCGA-LIHC)
===================================================================
Computes Spearman correlations between SOCE components and key
EMT / NF-κB partner genes in TCGA-LIHC (n ≈ 370 tumours).

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal (lihc_tcga)

Outputs
-------
  coexpression_results.csv    – ρ, p-value, FDR for each pair
  correlation_heatmap.png / .svg
  scatter_plots.png / .svg    – Individual scatter plots for top pairs

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn requests statsmodels
  python coexpression.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# ── Config ───────────────────────────────────────────────────────────────────
OUTDIR       = os.path.dirname(os.path.abspath(__file__))

SOCE_GENES   = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]

PARTNER_GENES = [
    # EMT TFs
    "VIM", "CDH1", "CDH2", "FN1",
    "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
    # NF-κB / IL-6
    "RELA", "RELB", "NFKB1", "IL6", "IL6R", "STAT3",
    # Drug resistance
    "ABCC2", "ABCB1", "ABCG2",
    # Other
    "TNF", "MAP3K7", "TGFB1", "MMP2", "MMP9",
]

ALL_GENES = SOCE_GENES + PARTNER_GENES


# ── 1. Download TCGA-LIHC expression ─────────────────────────────────────────
def fetch_expression(genes: list[str]) -> pd.DataFrame:
    """Returns (samples × genes) DataFrame from cBioPortal."""
    cache = os.path.join(OUTDIR, "_tcga_expr_cache.csv")
    if os.path.exists(cache):
        print(f"  [cache] Loading from {cache}")
        return pd.read_csv(cache, index_col=0)

    print("  Downloading TCGA-LIHC expression from cBioPortal …")
    base    = "https://www.cbioportal.org/api"
    profile = "lihc_tcga_rna_seq_v2_mrna"

    # Get sample IDs
    r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
    r.raise_for_status()
    sample_ids = r.json()

    # Get Entrez IDs
    r = requests.get(f"{base}/genes?geneIds={','.join(genes)}", timeout=60)
    r.raise_for_status()
    gene_map = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}

    all_data = {}
    for symbol, entrez in gene_map.items():
        payload = {
            "entrezGeneId":       entrez,
            "molecularProfileId": profile,
            "sampleIds":          sample_ids,
        }
        r = requests.post(
            f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
            json=payload, timeout=120,
        )
        if r.ok:
            for row in r.json():
                s = row["sampleId"]
                if s not in all_data:
                    all_data[s] = {}
                all_data[s][symbol] = row["value"]

    df = pd.DataFrame.from_dict(all_data, orient="index")
    df.to_csv(cache)
    print(f"  Saved expression cache: {df.shape}")
    return df


# ── 2. Spearman correlation for all SOCE × partner pairs ─────────────────────
def compute_correlations(df: pd.DataFrame) -> pd.DataFrame:
    """All pairwise Spearman r between SOCE genes and partner genes."""
    records = []
    for soce in SOCE_GENES:
        if soce not in df.columns:
            continue
        for partner in PARTNER_GENES:
            if partner not in df.columns:
                continue
            x = df[soce].dropna()
            y = df[partner].dropna()
            common = x.index.intersection(y.index)
            if len(common) < 30:
                continue
            rho, pval = stats.spearmanr(x[common], y[common])
            records.append({
                "SOCE_gene":   soce,
                "partner":     partner,
                "rho":         round(rho, 4),
                "p_value":     pval,
                "n":           len(common),
            })

    df_corr = pd.DataFrame(records)
    _, fdr, _, _ = multipletests(df_corr["p_value"], method="fdr_bh")
    df_corr["FDR"] = fdr.round(4)
    df_corr = df_corr.sort_values("rho", ascending=False)
    return df_corr


# ── 3. Heatmap ────────────────────────────────────────────────────────────────
def plot_heatmap(df_corr: pd.DataFrame, expr: pd.DataFrame):
    """Heatmap of Spearman ρ (SOCE × partners)."""
    # Build pivot matrix
    pivot = df_corr.pivot(index="SOCE_gene", columns="partner", values="rho")
    pivot = pivot.reindex(index=[g for g in SOCE_GENES if g in pivot.index])

    # Sort columns by absolute max correlation
    col_order = pivot.abs().max(axis=0).sort_values(ascending=False).index
    pivot = pivot[col_order]

    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 0.7), 5))
    sns.heatmap(
        pivot, ax=ax,
        cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
        annot=True, fmt=".2f", annot_kws={"size": 8},
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Spearman ρ"},
    )
    ax.set_title(
        "SOCE Component Co-expression Heatmap\n(TCGA-LIHC, n≈370, Spearman ρ)",
        fontsize=13, fontweight="bold"
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"correlation_heatmap.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved correlation_heatmap.png / .svg")


# ── 4. Scatter plots for key pairs ────────────────────────────────────────────
def plot_scatters(df_corr: pd.DataFrame, expr: pd.DataFrame):
    """Top-6 scatter plots (highest absolute ρ)."""
    top = df_corr.nlargest(6, "rho")
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    for ax, (_, row) in zip(axes, top.iterrows()):
        soce, partner = row["SOCE_gene"], row["partner"]
        if soce not in expr.columns or partner not in expr.columns:
            ax.set_visible(False); continue
        common = expr[soce].dropna().index.intersection(expr[partner].dropna().index)
        x = expr.loc[common, soce]
        y = expr.loc[common, partner]
        ax.scatter(x, y, alpha=0.35, s=15, c="#2980b9", edgecolors="none")
        # Regression line
        m, b = np.polyfit(x, y, 1)
        xline = np.linspace(x.min(), x.max(), 100)
        ax.plot(xline, m * xline + b, c="#c0392b", linewidth=2)
        ax.set_xlabel(f"{soce} (RSEM)", fontsize=10)
        ax.set_ylabel(f"{partner} (RSEM)", fontsize=10)
        ax.set_title(
            f"{soce} vs {partner}\nρ = {row['rho']:.3f}, FDR = {row['FDR']:.3f}",
            fontsize=10, fontweight="bold"
        )

    for ax in axes[len(top):]:
        ax.set_visible(False)

    fig.suptitle(
        "Top Co-expression Pairs — SOCE Components (TCGA-LIHC)",
        fontsize=13, fontweight="bold"
    )
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"scatter_plots.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved scatter_plots.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 05 — Spearman Co-expression Analysis")
    print("=" * 60)

    expr = fetch_expression(ALL_GENES)
    print(f"  Expression matrix: {expr.shape}")

    print("\n[1/3] Computing Spearman correlations …")
    df_corr = compute_correlations(expr)
    df_corr.to_csv(os.path.join(OUTDIR, "coexpression_results.csv"), index=False)
    print(f"  {len(df_corr)} pairs computed. Saved coexpression_results.csv")
    print("\n  Top correlations:")
    print(df_corr.head(10).to_string(index=False))

    print("\n[2/3] Heatmap …")
    plot_heatmap(df_corr, expr)

    print("\n[3/3] Scatter plots …")
    plot_scatters(df_corr, expr)

    print("\n  Done. All outputs saved to:", OUTDIR)
