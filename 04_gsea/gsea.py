#!/usr/bin/env python3
"""
04_gsea.py  –  Gene Set Enrichment Analysis (GSE140202 DEG)
============================================================
Runs GSEA on sorafenib-resistance DEG rankings from Module 02,
testing Hallmark gene sets from MSigDB.

Inputs
------
  ../02_geo_deg/deg_full_results.csv  (pre-ranked by t-statistic)

Method
------
  gseapy.prerank() on t-statistic ranking
  Gene sets : MSigDB Hallmarks (h.all.v2023.2.Hs.symbols.gmt)
              downloaded automatically by gseapy

Outputs
-------
  gsea_barplot.png / .svg
  gsea_barplot_results.csv

Usage
-----
  pip install pandas gseapy matplotlib seaborn
  python gsea.py  (run from this folder; Module 02 must have been run first)

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# ── Config ───────────────────────────────────────────────────────────────────
OUTDIR     = os.path.dirname(os.path.abspath(__file__))
DEG_FILE   = os.path.join(OUTDIR, "..", "02_geo_deg", "deg_full_results.csv")
GENE_SETS  = ["MSigDB_Hallmark_2020"]   # gseapy built-in gene set library
THREADS    = 4
PERMUTATIONS = 1000
FDR_THRESH = 0.25    # standard GSEA threshold

# ── 1. Load pre-ranked gene list ──────────────────────────────────────────────
def load_ranking(deg_file: str) -> pd.Series:
    """
    Returns a Series (index = gene, values = t-stat) sorted descending.
    Uses a t-statistic ranking when available (more robust than log2FC alone).
    """
    df = pd.read_csv(deg_file, index_col=0)
    if "t_stat" in df.columns:
        ranking = df["t_stat"]
    elif "t" in df.columns:
        if "gene" in df.columns:
            ranking = df.set_index("gene")["t"]
        else:
            ranking = df["t"]
    else:
        raise ValueError("deg_full_results.csv must contain 't_stat' or 't' column.")
    ranking = ranking.dropna().sort_values(ascending=False)
    ranking = ranking[~ranking.index.duplicated(keep="first")]
    return ranking


# ── 2. Run pre-ranked GSEA ────────────────────────────────────────────────────
def run_gsea(ranking: pd.Series) -> pd.DataFrame:
    """Runs gseapy.prerank and returns results DataFrame."""
    rnk = ranking.reset_index()
    rnk.columns = ["gene", "score"]

    results = gp.prerank(
        rnk              = rnk,
        gene_sets        = GENE_SETS,
        threads          = THREADS,
        permutation_num  = PERMUTATIONS,
        outdir           = None,
        seed             = 42,
        verbose          = False,
    )
    return results.res2d


# ── 3. Plot barplot ────────────────────────────────────────────────────────────
def barplot(df: pd.DataFrame):
    """Horizontal barplot coloured by direction (NES > 0 = enriched in resistant)."""
    # Filter and sort
    df = df.dropna(subset=["NES", "FDR q-val"])
    df = df[df["FDR q-val"] < FDR_THRESH].copy()
    if df.empty:
        # Relax threshold if nothing significant
        df = df.nsmallest(20, "FDR q-val")

    df = df.sort_values("NES")
    df["Term"] = df["Term"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ", regex=False).str.title()

    colours = ["#c0392b" if n > 0 else "#2980b9" for n in df["NES"]]

    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.4 + 1)))
    bars = ax.barh(df["Term"], df["NES"], color=colours, edgecolor="white", linewidth=0.5)

    # FDR labels
    for bar, (_, row) in zip(bars, df.iterrows()):
        fdr_label = f"FDR={row['FDR q-val']:.3f}" if row["FDR q-val"] >= 0.001 else f"FDR<0.001"
        x_pos = bar.get_width() + 0.05 if bar.get_width() > 0 else bar.get_width() - 0.05
        ha = "left" if bar.get_width() > 0 else "right"
        ax.text(x_pos, bar.get_y() + bar.get_height() / 2, fdr_label,
                va="center", ha=ha, fontsize=8)

    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalised Enrichment Score (NES)", fontsize=11)
    ax.set_title(
        "GSEA Hallmark Pathways — Sorafenib-Resistant vs. Sensitive Huh7\n"
        f"(GSE140202, FDR < {FDR_THRESH})",
        fontsize=12, fontweight="bold"
    )
    from matplotlib.patches import Patch
    legend = [
        Patch(color="#c0392b", label="Enriched in Resistant"),
        Patch(color="#2980b9", label="Enriched in Sensitive"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=9)
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"gsea_barplot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved gsea_barplot.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 04 — GSEA on GSE140202 DEG Rankings")
    print("=" * 60)

    print(f"\n[1/3] Loading ranking from {DEG_FILE} …")
    if not os.path.exists(DEG_FILE):
        print(f"  ⚠  {DEG_FILE} not found. Run Module 02 first.")
        import sys; sys.exit(1)
    ranking = load_ranking(DEG_FILE)
    print(f"  Loaded {len(ranking)} genes.")

    print(f"\n[2/3] Running pre-ranked GSEA ({PERMUTATIONS} permutations) …")
    res_df = run_gsea(ranking)
    res_df.to_csv(os.path.join(OUTDIR, "gsea_barplot_results.csv"), index=False)
    print(f"  Tested {len(res_df)} gene sets. Saved gsea_barplot_results.csv")

    sig = res_df[res_df["FDR q-val"] < FDR_THRESH].sort_values("NES", ascending=False)
    print(f"\n  Significant pathways (FDR < {FDR_THRESH}):")
    if not sig.empty:
        print(sig[["Term", "NES", "FDR q-val", "ES"]].to_string(index=False))
    else:
        print("  None at this threshold; check raw results file.")

    print("\n[3/3] Barplot …")
    barplot(res_df)

    print("\n  Done. All outputs saved to:", OUTDIR)
