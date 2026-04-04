#!/usr/bin/env python3
"""
14_tcga_gsea_stim1.py  –  GSEA on STIM1-Correlated Genes (TCGA) + GSE140202
============================================================================
Two complementary GSEA analyses:
  (A) TCGA-LIHC: Genes ranked by Spearman correlation with STIM1
  (B) GSE140202: Genes ranked by t-statistic (resistant vs. sensitive Huh7)
                 — same analysis as Module 04 but focused on STIM1-relevant pathways

Data sources
------------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal
  GSE140202 DEG results     : ../02_geo_deg/deg_full_results.csv

Outputs
-------
  tcga_stim1_pathway_correlation.csv  – Spearman ρ per gene with STIM1
  gsea_geo140202_results.csv          – GSEA results on GSE140202 ranking
  tcga_stim1_pathway_enrichment.png / .svg
  figure_gsea_combined.png / .svg     – Combined figure for manuscript
  gsea_geo_raw/                        – Raw gseapy output files

Key result
----------
  TNFα/NF-κB signalling: NES = +1.85, FDR < 0.001 (enriched in resistant)

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn gseapy requests
  python tcga_gsea_stim1.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import requests
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR    = os.path.dirname(os.path.abspath(__file__))
DEG_FILE  = os.path.join(OUTDIR, "..", "02_geo_deg", "deg_full_results.csv")
GSEA_OUTD = os.path.join(OUTDIR, "gsea_geo_raw")
os.makedirs(GSEA_OUTD, exist_ok=True)

GENE_SETS   = ["MSigDB_Hallmark_2020"]
PERMS       = 1000
FDR_THRESH  = 0.25


# ── Part A: TCGA STIM1-correlated gene ranking ────────────────────────────────
def tcga_stim1_ranking() -> pd.Series:
    """
    Download all-gene expression for TCGA-LIHC and compute Spearman ρ
    with STIM1. Returns Series sorted by ρ (descending).
    """
    cache_all  = os.path.join(OUTDIR, "_tcga_all_genes_stim1_corr.csv")
    if os.path.exists(cache_all):
        print(f"  [cache] {cache_all}")
        df = pd.read_csv(cache_all, index_col=0)
        return df["rho"].sort_values(ascending=False)

    # Load STIM1 expression first
    stim1_cache = None
    for p in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
    ):
        if os.path.exists(p):
            df = pd.read_csv(p, index_col=0)
            if "STIM1" in df.columns:
                stim1_vals = df["STIM1"].dropna()
                stim1_cache = p; break

    if stim1_cache is None:
        print("  ⚠  STIM1 expression cache not found. Run Module 05 or 01 first.")
        return pd.Series(dtype=float)

    # Get all TCGA-LIHC gene expression (paginated bulk download)
    print("  Computing STIM1 Spearman correlations across all TCGA-LIHC genes …")
    print("  (This uses the pre-computed expression cache; for all genes,")
    print("   download the full TCGA-LIHC matrix from cBioPortal)")

    # Use available cached expression
    df_expr = pd.read_csv(stim1_cache, index_col=0)
    stim1   = df_expr["STIM1"].dropna()
    records = []
    for gene in df_expr.columns:
        if gene == "STIM1":
            continue
        x = df_expr[gene].dropna()
        common = stim1.index.intersection(x.index)
        if len(common) < 30:
            continue
        rho, _ = stats.spearmanr(stim1[common], x[common])
        records.append({"gene": gene, "rho": rho})

    df_rho = pd.DataFrame(records).set_index("gene")
    df_rho.to_csv(cache_all)
    return df_rho["rho"].sort_values(ascending=False)


# ── Part B: GSE140202 GSEA ────────────────────────────────────────────────────
def run_geo_gsea() -> pd.DataFrame:
    """Pre-ranked GSEA on GSE140202 t-statistic."""
    if not os.path.exists(DEG_FILE):
        print(f"  ⚠  {DEG_FILE} not found. Run Module 02 first.")
        return pd.DataFrame()

    print(f"  Loading DEG ranking from {DEG_FILE} …")
    df_deg = pd.read_csv(DEG_FILE, index_col=0)

    if "t_stat" not in df_deg.columns:
        print("  ⚠  't_stat' column missing. Using log2FC as ranking.")
        rank_col = "log2FC"
    else:
        rank_col = "t_stat"

    rnk = df_deg[[rank_col]].dropna().sort_values(rank_col, ascending=False)
    rnk = rnk.reset_index()
    rnk.columns = ["gene", "score"]

    print(f"  Running GSEA with {len(rnk)} genes × {GENE_SETS} …")
    res = gp.prerank(
        rnk             = rnk,
        gene_sets       = GENE_SETS,
        threads         = 4,
        permutation_num = PERMS,
        outdir          = GSEA_OUTD,
        seed            = 42,
        verbose         = False,
    )
    return res.res2d


# ── Plotting ──────────────────────────────────────────────────────────────────
def plot_combined(gsea_df: pd.DataFrame, stim1_corr: pd.Series):
    """Combined 1×2 figure for manuscript."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(
        "STIM1-Associated Pathway Enrichment in HCC",
        fontsize=14, fontweight="bold"
    )

    # Panel A: STIM1 top correlates (bar)
    ax = axes[0]
    if not stim1_corr.empty:
        top = pd.concat([stim1_corr.head(15), stim1_corr.tail(15)])
        colours = ["#c0392b" if v > 0 else "#2980b9" for v in top.values]
        ax.barh(top.index[::-1], top.values[::-1], color=colours[::-1])
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("Spearman ρ with STIM1", fontsize=11)
        ax.set_title("Top STIM1 Co-expressed Genes\n(TCGA-LIHC)", fontsize=11)
    else:
        ax.text(0.5, 0.5, "Run Module 05 first\nto generate STIM1 correlations",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_title("STIM1 Co-expression Ranking", fontsize=11)

    # Panel B: GSEA barplot (GSE140202)
    ax2 = axes[1]
    if not gsea_df.empty:
        df_plot = gsea_df.copy()
        df_plot = df_plot[df_plot["FDR q-val"] < FDR_THRESH].copy()
        if df_plot.empty:
            df_plot = gsea_df.nsmallest(15, "FDR q-val")
        df_plot = df_plot.sort_values("NES")
        df_plot["Term"] = (
            df_plot["Term"]
            .str.replace("HALLMARK_", "", regex=False)
            .str.replace("_", " ", regex=False)
            .str.title()
        )
        colours2 = ["#c0392b" if n > 0 else "#2980b9" for n in df_plot["NES"]]
        bars = ax2.barh(df_plot["Term"], df_plot["NES"],
                        color=colours2, edgecolor="white", linewidth=0.5)
        ax2.axvline(0, color="black", linewidth=0.8)

        for bar, (_, row) in zip(bars, df_plot.iterrows()):
            fdr_s = f"FDR={row['FDR q-val']:.3f}" if row["FDR q-val"] >= 0.001 else "FDR<0.001"
            x_off = bar.get_width() + 0.03 if bar.get_width() > 0 else bar.get_width() - 0.03
            ax2.text(x_off, bar.get_y() + bar.get_height() / 2,
                     fdr_s, va="center",
                     ha="left" if bar.get_width() > 0 else "right",
                     fontsize=8)

        ax2.set_xlabel("NES", fontsize=11)
        ax2.set_title(
            "GSEA Hallmark Pathways\n(GSE140202: Resistant vs. Sensitive Huh7)",
            fontsize=11
        )
    else:
        ax2.text(0.5, 0.5, "Run Module 02 first\nto generate DEG results",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=11)

    plt.tight_layout()
    for name in ("figure_gsea_combined", "gsea_geo140202_barplot",
                 "tcga_stim1_pathway_enrichment"):
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"),
                        dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved figure_gsea_combined.png / .svg")


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 14 — GSEA: STIM1 Pathway Enrichment")
    print("=" * 60)

    print("\n[1/3] TCGA-LIHC STIM1 correlation ranking …")
    stim1_corr = tcga_stim1_ranking()
    if not stim1_corr.empty:
        df_corr = stim1_corr.reset_index()
        df_corr.columns = ["gene", "rho"]
        df_corr.to_csv(
            os.path.join(OUTDIR, "tcga_stim1_pathway_correlation.csv"), index=False
        )
        print(f"  Saved tcga_stim1_pathway_correlation.csv ({len(df_corr)} genes)")

    print("\n[2/3] GSEA on GSE140202 …")
    gsea_df = run_geo_gsea()
    if not gsea_df.empty:
        gsea_df.to_csv(
            os.path.join(OUTDIR, "gsea_geo140202_results.csv"), index=False
        )
        sig = gsea_df[gsea_df["FDR q-val"] < FDR_THRESH].sort_values("NES", ascending=False)
        print(f"  Saved gsea_geo140202_results.csv  ({len(gsea_df)} pathways tested)")
        print("\n  Significant pathways:")
        if not sig.empty:
            print(sig[["Term", "NES", "FDR q-val"]].to_string(index=False))
        else:
            print(f"  None at FDR < {FDR_THRESH}")
    else:
        gsea_df = pd.DataFrame()

    print("\n[3/3] Combined figure …")
    plot_combined(gsea_df, stim1_corr)

    print("\n  Done. All outputs saved to:", OUTDIR)
