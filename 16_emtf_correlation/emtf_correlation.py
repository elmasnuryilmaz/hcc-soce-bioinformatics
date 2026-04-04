#!/usr/bin/env python3
"""
16_emtf_correlation.py  –  EMT Transcription Factor Correlation Bubble/Dot Plot
================================================================================
Visualises the Spearman correlations between SOCE components (STIM1, TRPC6,
TRPC1, ORAI1) and a comprehensive panel of EMT transcription factors as a
publication-quality bubble dot plot.

Bubble size    = |ρ|  (larger = stronger)
Bubble colour  = ρ    (red = positive, blue = negative)
Significance   = FDR < 0.05 highlighted with black border

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal  (lihc_tcga)

Outputs
-------
  emtf_correlation_data.csv          – Full ρ table with FDR
  emtf_bubble_dotplot.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn requests statsmodels
  python emtf_correlation.py

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
import matplotlib.colors as mcolors
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR     = os.path.dirname(os.path.abspath(__file__))
SOCE_GENES = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]

# Comprehensive EMT transcription factor panel
EMTF_GENES = [
    # Core EMT TFs
    "SNAI1", "SNAI2", "TWIST1", "TWIST2",
    "ZEB1",  "ZEB2",
    # Mesenchymal markers
    "VIM", "FN1", "CDH2",
    # Epithelial markers (expected negative ρ)
    "CDH1", "EPCAM", "KRT18", "KRT19",
    # Accessory TFs
    "FOXC1", "FOXC2", "PRRX1", "PRRX2",
    "TCF4",  "TCF3",
    # MMP / invasion
    "MMP2", "MMP9", "MMP14",
    # TGF-β pathway
    "TGFB1", "TGFB2", "SMAD2", "SMAD3", "SMAD7",
    # Additional
    "AXL", "GAS6", "PDGFRA",
]

ALL_GENES = list(dict.fromkeys(SOCE_GENES + EMTF_GENES))


def load_expression() -> pd.DataFrame:
    for cache in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "15_clean_correlation", "_tcga_expr_cache.csv"),
    ):
        if os.path.exists(cache):
            df = pd.read_csv(cache, index_col=0)
            avail = [g for g in ALL_GENES if g in df.columns]
            if len(avail) > len(SOCE_GENES):
                print(f"  [cache] {cache}  ({len(avail)} / {len(ALL_GENES)} genes)")
                return df[avail]

    print("  Downloading expression from cBioPortal …")
    base    = "https://www.cbioportal.org/api"
    profile = "lihc_tcga_rna_seq_v2_mrna"
    r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
    r.raise_for_status()
    sids = r.json()
    r = requests.get(f"{base}/genes?geneIds={','.join(ALL_GENES)}", timeout=60)
    r.raise_for_status()
    gmap = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
    data = {}
    for sym, eid in gmap.items():
        pl = {"entrezGeneId": eid, "molecularProfileId": profile, "sampleIds": sids}
        r = requests.post(f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
                          json=pl, timeout=120)
        for row in r.json():
            s = row["sampleId"]
            data.setdefault(s, {})[sym] = row["value"]
    return pd.DataFrame.from_dict(data, orient="index")


def compute_corr(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for soce in SOCE_GENES:
        if soce not in df.columns:
            continue
        for emtf in EMTF_GENES:
            if emtf not in df.columns:
                continue
            common = df[soce].dropna().index.intersection(df[emtf].dropna().index)
            if len(common) < 30:
                continue
            rho, p = stats.spearmanr(df.loc[common, soce], df.loc[common, emtf])
            records.append({"SOCE": soce, "EMTF": emtf,
                             "rho": round(rho, 4), "p_value": p, "n": len(common)})

    df_res = pd.DataFrame(records)
    if not df_res.empty:
        _, fdr, _, _ = multipletests(df_res["p_value"], method="fdr_bh")
        df_res["FDR"] = fdr.round(4)
    return df_res


def bubble_dotplot(df_corr: pd.DataFrame):
    """
    Bubble dot plot:
      x-axis  = EMT TF gene
      y-axis  = SOCE gene
      size    = |ρ| (scaled)
      colour  = ρ  (RdBu_r)
      border  = black if FDR < 0.05, grey otherwise
    """
    # Build pivot matrix for ordering
    pivot_rho = df_corr.pivot(index="SOCE", columns="EMTF", values="rho")
    pivot_fdr = df_corr.pivot(index="SOCE", columns="EMTF", values="FDR")

    # Order EMT TFs by mean |ρ| across SOCE genes
    avail_emtf = [g for g in EMTF_GENES if g in pivot_rho.columns]
    col_order   = (pivot_rho[avail_emtf].abs().mean(axis=0)
                   .sort_values(ascending=False).index.tolist())
    row_order   = [g for g in SOCE_GENES if g in pivot_rho.index]

    pivot_rho = pivot_rho.reindex(index=row_order, columns=col_order)
    pivot_fdr = pivot_fdr.reindex(index=row_order, columns=col_order)

    cmap   = plt.cm.RdBu_r
    norm   = mcolors.Normalize(vmin=-0.6, vmax=0.6)
    max_sz = 800   # max bubble area in points²

    fig, ax = plt.subplots(figsize=(max(14, len(col_order) * 0.55), 5))

    for yi, soce in enumerate(row_order):
        for xi, emtf in enumerate(col_order):
            rho = pivot_rho.loc[soce, emtf] if emtf in pivot_rho.columns else np.nan
            fdr = pivot_fdr.loc[soce, emtf] if emtf in pivot_fdr.columns else 1.0
            if np.isnan(rho):
                continue
            colour    = cmap(norm(rho))
            size      = (abs(rho) ** 1.5) * max_sz
            edgecol   = "black" if fdr < 0.05 else "#cccccc"
            linewidth = 1.5 if fdr < 0.05 else 0.5
            ax.scatter(xi, yi, s=size, c=[colour],
                       edgecolors=edgecol, linewidths=linewidth,
                       zorder=3, alpha=0.9)

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order, fontsize=11)

    # Grid
    ax.set_xlim(-0.6, len(col_order) - 0.4)
    ax.set_ylim(-0.6, len(row_order) - 0.4)
    ax.set_axisbelow(True)
    ax.grid(True, linewidth=0.4, color="#dddddd")

    # Colour bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, shrink=0.7)
    cbar.set_label("Spearman ρ", fontsize=10)

    # Size legend
    for rho_ex in [0.2, 0.4, 0.6]:
        sz = (rho_ex ** 1.5) * max_sz
        ax.scatter([], [], s=sz, c="grey", alpha=0.6, label=f"|ρ| = {rho_ex}")
    ax.legend(title="Bubble size", loc="lower right", fontsize=8,
              framealpha=0.9, title_fontsize=9)

    ax.set_title(
        "SOCE Component Correlation with EMT Transcription Factors\n"
        "(TCGA-LIHC, n ≈ 370, Spearman ρ; black border = FDR < 0.05)",
        fontsize=12, fontweight="bold"
    )

    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"emtf_bubble_dotplot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved emtf_bubble_dotplot.png / .svg")


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 16 — EMT-TF Correlation Bubble Dot Plot")
    print("=" * 60)

    print("\n[1/3] Loading expression …")
    df = load_expression()
    print(f"  {df.shape[0]} samples")

    print("\n[2/3] Computing Spearman correlations …")
    df_corr = compute_corr(df)
    df_corr.to_csv(os.path.join(OUTDIR, "emtf_correlation_data.csv"), index=False)
    print(f"  {len(df_corr)} gene pairs. Saved emtf_correlation_data.csv")

    print("\n  Top positive correlations:")
    print(df_corr.nlargest(10, "rho")[["SOCE", "EMTF", "rho", "FDR"]].to_string(index=False))

    print("\n[3/3] Bubble dot plot …")
    bubble_dotplot(df_corr)

    print("\n  Done. All outputs saved to:", OUTDIR)
