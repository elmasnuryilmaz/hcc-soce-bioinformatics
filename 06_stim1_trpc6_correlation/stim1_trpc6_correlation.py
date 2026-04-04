#!/usr/bin/env python3
"""
06_stim1_trpc6_correlation.py  –  STIM1 vs TRPC6 Direct Correlation
=====================================================================
Scatter plot and Spearman/Pearson correlation between STIM1 and TRPC6
expression in TCGA-LIHC tumours, plus stratification analysis.

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal (lihc_tcga)

Outputs
-------
  stim1_trpc6_analysis.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn requests
  python stim1_trpc6_correlation.py

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
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR = os.path.dirname(os.path.abspath(__file__))
GENES  = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]


def fetch_expression(genes: list[str]) -> pd.DataFrame:
    cache = os.path.join(OUTDIR, "_tcga_expr_cache.csv")
    # Re-use cache from module 05 if present
    cache2 = os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv")
    for c in (cache, cache2):
        if os.path.exists(c):
            df = pd.read_csv(c, index_col=0)
            if all(g in df.columns for g in genes):
                print(f"  [cache] {c}")
                return df[genes]

    print("  Downloading TCGA-LIHC expression …")
    base    = "https://www.cbioportal.org/api"
    profile = "lihc_tcga_rna_seq_v2_mrna"
    r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
    r.raise_for_status()
    sample_ids = r.json()

    r = requests.get(f"{base}/genes?geneIds={','.join(genes)}", timeout=60)
    r.raise_for_status()
    gene_map = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}

    data = {}
    for symbol, entrez in gene_map.items():
        payload = {"entrezGeneId": entrez, "molecularProfileId": profile, "sampleIds": sample_ids}
        r = requests.post(f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
                          json=payload, timeout=120)
        for row in r.json():
            s = row["sampleId"]
            if s not in data:
                data[s] = {}
            data[s][symbol] = row["value"]

    df = pd.DataFrame.from_dict(data, orient="index")
    df.to_csv(cache)
    return df


def make_figure(df: pd.DataFrame):
    """
    2×2 panel:
      [A] STIM1 vs TRPC6 scatter (coloured by expression density)
      [B] STIM1 vs TRPC1 scatter
      [C] STIM1 vs ORAI1 scatter
      [D] Quadrant analysis: dual-high / dual-low / discordant
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    fig.suptitle(
        "STIM1 Co-expression with SOCE Partners in TCGA-LIHC (n≈370)",
        fontsize=14, fontweight="bold"
    )

    pairs = [("STIM1", "TRPC6"), ("STIM1", "TRPC1"), ("STIM1", "ORAI1")]
    for ax, (gx, gy) in zip(axes.flatten()[:3], pairs):
        common = df[gx].dropna().index.intersection(df[gy].dropna().index)
        x = df.loc[common, gx]
        y = df.loc[common, gy]

        # Kernel density scatter
        from scipy.stats import gaussian_kde
        xy = np.vstack([x, y])
        kde = gaussian_kde(xy)(xy)

        scatter = ax.scatter(x, y, c=kde, s=18, cmap="viridis", alpha=0.7, linewidths=0)
        plt.colorbar(scatter, ax=ax, label="Density")

        rho, p_rho = stats.spearmanr(x, y)
        r,   p_r   = stats.pearsonr(x, y)
        m, b = np.polyfit(x, y, 1)
        xline = np.linspace(x.min(), x.max(), 100)
        ax.plot(xline, m * xline + b, c="#e74c3c", linewidth=2)

        ax.set_xlabel(f"{gx} (RSEM)", fontsize=11)
        ax.set_ylabel(f"{gy} (RSEM)", fontsize=11)
        pval_str = f"p = {p_rho:.3f}" if p_rho >= 0.001 else f"p = {p_rho:.2e}"
        ax.set_title(f"{gx} vs {gy}\nSpearman ρ = {rho:.3f}, {pval_str}", fontsize=11)

    # Panel D: Quadrant analysis (STIM1 high/low × TRPC6 high/low)
    ax = axes[1, 1]
    common = df["STIM1"].dropna().index.intersection(df["TRPC6"].dropna().index)
    x = df.loc[common, "STIM1"]
    y = df.loc[common, "TRPC6"]
    med_x, med_y = x.median(), y.median()

    colours_q = np.where(
        (x >= med_x) & (y >= med_y), "#c0392b",   # dual-high
        np.where(
            (x <  med_x) & (y <  med_y), "#2980b9",  # dual-low
            "#95a5a6"                               # discordant
        )
    )
    ax.scatter(x, y, c=colours_q, s=20, alpha=0.7, linewidths=0)
    ax.axvline(med_x, color="grey", linestyle="--", linewidth=1)
    ax.axhline(med_y, color="grey", linestyle="--", linewidth=1)

    n_dh = ((x >= med_x) & (y >= med_y)).sum()
    n_dl = ((x <  med_x) & (y <  med_y)).sum()
    n_di = len(x) - n_dh - n_dl

    ax.text(x.max() * 0.7, y.max() * 0.95,
            f"Dual-high\n(n={n_dh})", color="#c0392b", fontsize=9, ha="center")
    ax.text(x.min() * 1.1, y.min() + (y.max() - y.min()) * 0.05,
            f"Dual-low\n(n={n_dl})", color="#2980b9", fontsize=9)

    ax.set_xlabel("STIM1 (RSEM)", fontsize=11)
    ax.set_ylabel("TRPC6 (RSEM)", fontsize=11)
    ax.set_title("STIM1 × TRPC6 Quadrant Analysis\n(split at median)", fontsize=11)

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"stim1_trpc6_analysis.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved stim1_trpc6_analysis.png / .svg")


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 06 — STIM1 vs TRPC6 Correlation")
    print("=" * 60)

    df = fetch_expression(GENES)
    print(f"  {df.shape[0]} samples, {df.shape[1]} genes")
    make_figure(df)

    # Print key correlations
    for g in ["TRPC6", "TRPC1", "ORAI1"]:
        if g in df.columns:
            common = df["STIM1"].dropna().index.intersection(df[g].dropna().index)
            rho, p = stats.spearmanr(df.loc[common, "STIM1"], df.loc[common, g])
            print(f"  STIM1 vs {g}: ρ = {rho:.3f}, p = {p:.4f}")

    print("\n  Done. All outputs saved to:", OUTDIR)
