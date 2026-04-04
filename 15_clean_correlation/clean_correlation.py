#!/usr/bin/env python3
"""
15_clean_correlation.py  –  SOCE × EMT/NF-κB Clean Correlation Heatmap
=======================================================================
Publication-quality Spearman correlation heatmap between SOCE components
and a curated panel of EMT and NF-κB pathway genes in TCGA-LIHC.
This module produces the final cleaned version used in the manuscript figure.

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM : cBioPortal  (lihc_tcga)
  Coexpression results       : ../05_coexpression/coexpression_results.csv

Key results (manuscript)
------------------------
  TRPC6–VIM  : ρ = 0.48  (p < 0.001)
  TRPC6–ZEB2 : ρ = 0.45  (p < 0.001)
  STIM1–IL6R : ρ = 0.41  (p < 0.001)
  STIM1–RELA : ρ = 0.25  (p < 0.001)

Outputs
-------
  soce_emtf_corr_matrix.csv              – ρ matrix (SOCE × partners)
  soce_partner_correlations.csv          – Long-form table with p-values
  soce_correlation_heatmap_clean.png / .svg
  soce_correlation_heatmap_final.png / .svg
  soce_correlation_heatmap_v2.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn requests
  python clean_correlation.py

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
import matplotlib.patches as mpatches
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR     = os.path.dirname(os.path.abspath(__file__))
SOCE_GENES = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]

# Curated partner gene panel for the figure
EMT_GENES  = ["VIM", "CDH1", "CDH2", "FN1", "ZEB1", "ZEB2",
              "SNAI1", "SNAI2", "TWIST1", "TWIST2", "MMP2", "MMP9"]
NFKB_GENES = ["RELA", "RELB", "NFKB1", "TNF", "IL6", "IL6R",
              "STAT3", "MAP3K7", "ABCC2", "ABCB1"]

PARTNER_GENES = EMT_GENES + NFKB_GENES
ALL_GENES     = SOCE_GENES + PARTNER_GENES

GENE_CATEGORY = (
    {g: "EMT" for g in EMT_GENES} |
    {g: "NF-κB / IL-6" for g in NFKB_GENES} |
    {g: "SOCE" for g in SOCE_GENES}
)


def load_expression() -> pd.DataFrame:
    """Load TCGA-LIHC expression from sibling cache or download."""
    for cache in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
    ):
        if os.path.exists(cache):
            df = pd.read_csv(cache, index_col=0)
            avail = [g for g in ALL_GENES if g in df.columns]
            if len(avail) > len(SOCE_GENES):
                print(f"  [cache] {cache}  ({len(avail)} genes available)")
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
    df = pd.DataFrame.from_dict(data, orient="index")
    os.makedirs(OUTDIR, exist_ok=True)
    df.to_csv(os.path.join(OUTDIR, "_tcga_expr_cache.csv"))
    return df


def compute_matrix(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (rho_matrix, long_form_table)."""
    records = []
    rho_matrix = pd.DataFrame(index=SOCE_GENES, columns=PARTNER_GENES, dtype=float)

    for soce in SOCE_GENES:
        if soce not in df.columns:
            continue
        for partner in PARTNER_GENES:
            if partner not in df.columns:
                continue
            common = df[soce].dropna().index.intersection(df[partner].dropna().index)
            if len(common) < 30:
                continue
            rho, p = stats.spearmanr(df.loc[common, soce], df.loc[common, partner])
            rho_matrix.loc[soce, partner] = rho
            records.append({"SOCE": soce, "partner": partner,
                             "rho": round(rho, 4), "p_value": p, "n": len(common)})

    long_df = pd.DataFrame(records)
    if not long_df.empty:
        _, fdr, _, _ = multipletests(long_df["p_value"], method="fdr_bh")
        long_df["FDR"] = fdr.round(4)
        long_df = long_df.sort_values("rho", ascending=False)

    return rho_matrix.astype(float), long_df


def heatmap_clean(rho: pd.DataFrame, long_df: pd.DataFrame, suffix: str = "clean"):
    """Annotated heatmap with category colour bar."""
    # Column annotations by category
    col_cats = [GENE_CATEGORY.get(g, "other") for g in rho.columns]
    cat_colours = {"EMT": "#e67e22", "NF-κB / IL-6": "#8e44ad", "SOCE": "#c0392b", "other": "#95a5a6"}

    # Sort columns: EMT genes first, then NF-κB
    emt_cols  = [g for g in rho.columns if g in EMT_GENES]
    nfkb_cols = [g for g in rho.columns if g in NFKB_GENES]
    sorted_cols = emt_cols + nfkb_cols
    rho_plot = rho[sorted_cols]

    # FDR significance mask
    sig = pd.DataFrame("", index=rho_plot.index, columns=rho_plot.columns)
    if not long_df.empty:
        for _, row in long_df.iterrows():
            if row["SOCE"] in sig.index and row["partner"] in sig.columns:
                stars = ("***" if row["FDR"] < 0.001
                         else ("**" if row["FDR"] < 0.01
                               else ("*" if row["FDR"] < 0.05 else "")))
                if stars:
                    sig.loc[row["SOCE"], row["partner"]] = stars

    fig, axes = plt.subplots(
        2, 1, figsize=(max(14, len(sorted_cols) * 0.65), 7),
        gridspec_kw={"height_ratios": [0.12, 1]}, sharex=True
    )

    # Top annotation bar
    ax_ann = axes[0]
    cat_bar = [cat_colours.get(GENE_CATEGORY.get(g, "other"), "#aaa") for g in sorted_cols]
    ax_ann.bar(range(len(sorted_cols)), [1] * len(sorted_cols),
               color=cat_bar, width=1.0, align="center")
    ax_ann.set_xlim(-0.5, len(sorted_cols) - 0.5)
    ax_ann.set_yticks([])
    ax_ann.set_xticks([])
    ax_ann.set_ylabel("Category", fontsize=8, rotation=0, labelpad=45)
    for patch in [
        mpatches.Patch(color=cat_colours["EMT"], label="EMT"),
        mpatches.Patch(color=cat_colours["NF-κB / IL-6"], label="NF-κB / IL-6"),
    ]:
        ax_ann.add_patch(patch)
    ax_ann.legend(
        handles=[mpatches.Patch(color=v, label=k) for k, v in cat_colours.items()
                 if k != "other"],
        loc="upper right", fontsize=8, framealpha=0.9
    )

    # Heatmap
    ax = axes[1]
    sns.heatmap(
        rho_plot, ax=ax,
        cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
        annot=rho_plot.round(2), fmt=".2f", annot_kws={"size": 8},
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Spearman ρ", "shrink": 0.6},
        xticklabels=True, yticklabels=True,
    )
    # Add significance stars
    for i, soce in enumerate(rho_plot.index):
        for j, partner in enumerate(sorted_cols):
            stars = sig.loc[soce, partner]
            if stars:
                ax.text(j + 0.5, i + 0.15, stars, ha="center", va="bottom",
                        fontsize=7, color="black", fontweight="bold")

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)
    ax.set_title(
        "SOCE Component Co-expression with EMT and NF-κB Pathway Genes\n"
        "(TCGA-LIHC, n ≈ 370, Spearman ρ; * FDR<0.05, ** FDR<0.01, *** FDR<0.001)",
        fontsize=12, fontweight="bold", pad=10
    )

    plt.tight_layout()
    for fname in (
        f"soce_correlation_heatmap_{suffix}",
        "soce_correlation_heatmap_final",
        "soce_correlation_heatmap_v2",
    ):
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(OUTDIR, f"{fname}.{ext}"),
                        dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved soce_correlation_heatmap_clean / _final / _v2  (.png + .svg)")


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 15 — Clean SOCE Correlation Heatmap")
    print("=" * 60)

    print("\n[1/3] Loading expression …")
    df = load_expression()
    print(f"  {df.shape[0]} samples, {df.shape[1]} genes available")

    print("\n[2/3] Computing Spearman correlations …")
    rho_matrix, long_df = compute_matrix(df)
    rho_matrix.to_csv(os.path.join(OUTDIR, "soce_emtf_corr_matrix.csv"))
    long_df.to_csv(os.path.join(OUTDIR, "soce_partner_correlations.csv"), index=False)
    print("  Saved soce_emtf_corr_matrix.csv + soce_partner_correlations.csv")

    # Print key findings
    key_pairs = [("TRPC6", "VIM"), ("TRPC6", "ZEB2"),
                 ("STIM1", "IL6R"), ("STIM1", "RELA")]
    print("\n  Key correlations:")
    for soce, partner in key_pairs:
        row = long_df[(long_df["SOCE"] == soce) & (long_df["partner"] == partner)]
        if not row.empty:
            print(f"    {soce}–{partner}: ρ = {row.iloc[0]['rho']:.3f}, "
                  f"FDR = {row.iloc[0]['FDR']:.4f}")

    print("\n[3/3] Heatmap …")
    heatmap_clean(rho_matrix, long_df)

    print("\n  Done. All outputs saved to:", OUTDIR)
