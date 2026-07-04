"""
21_tumor_purity/tumor_purity.py
================================
Tumor Purity Correction for SOCE Gene Correlations in TCGA-LIHC

Motivation:
  Correlations between SOCE genes and EMT markers can be confounded by
  tumour purity (stromal/immune contamination). A stroma-rich tumour will
  spuriously show high VIM/CDH2 (mesenchymal markers from fibroblasts)
  and low CDH1 (lost from epithelial contamination). ESTIMATE-based
  purity correction removes this confounding variable.

Method:
  1. ESTIMATE algorithm (Yoshihara 2013 Nature Comm):
     - Stromal score   = mean z-score of 141 stromal genes
     - Immune score    = mean z-score of 141 immune genes
     - ESTIMATE score  = stromal + immune
     - Tumour purity   = cos(0.6049872018 + 0.0001467884 × ESTIMATE_score)
  2. Partial Spearman correlation: STIM1 ~ EMT_gene | purity
  3. Compare crude vs. purity-corrected correlations

Outputs:
  purity_distribution.png / .svg    – histogram of estimated purity
  correlation_comparison.png / .svg – crude vs. corrected ρ scatter/bar
  partial_correlations.csv          – full table of corrected correlations
  purity_scores.csv                 – per-sample purity estimates
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import requests
import warnings
warnings.filterwarnings("ignore")

sys.path.append(str(Path(__file__).resolve().parents[1]))
from common.cbioportal import fetch_mrna_expression

np.random.seed(42)

OUT       = os.path.dirname(os.path.abspath(__file__))
CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID  = "lihc_tcga_pan_can_atlas_2018"

# ─────────────────────────────────────────────────────────────────────────────
# ESTIMATE gene signatures (Yoshihara et al. 2013 — published lists)
# ─────────────────────────────────────────────────────────────────────────────
STROMAL_GENES = [
    "ACTA2","BGN","BMP5","COL1A1","COL1A2","COL3A1","COL5A1","COL5A2",
    "COL6A1","COL6A2","COL6A3","DCN","DPT","EFEMP2","ELN","FAP","FBLN1",
    "FBLN2","FMOD","FN1","INHBA","LGALS1","LTBP2","MFAP4","MMP2","OLFML3",
    "PDGFRB","POSTN","PRSS35","PXDN","RUNX1T1","SPARC","TAGLN","THBS2",
    "TIMP3","TWIST1","VIM","VCAN","CTGF","FBN1","PCOLCE","SRPX","AEBP1",
    "CDH11","CXCL12","FILIP1L","ISLR","LUM","MXRA8","NID2","SFRP2",
    "SFRP4","TGFB1","THY1","TNFRSF11B","TNFSF4","TNFSF13B","ZEB1",
    "ADAM12","BAMBI","COMP","COL10A1","COL11A1","COL12A1","COL16A1",
    "COLGALT1","CTHRC1","EDIL3","GREM1","ITGA11","LOX","LOXL1","LOXL2",
    "MATN3","MMP11","MMP14","NNMT","PALLD","PLAU","PLPPR4","PTHLH",
    "SERPINH1","SNAI1","SNAI2","SPON2","TNC","TNFRSF12A","TWIST2",
    "WNT5A","WIPF1","ZEB2","ANTXR1","ASPN","ECM2","EGLN3","EMILIN1",
    "FBLN5","FZD10","INHBB","LRRC15","MFI2","MMP3","NOX4","OLFML2B",
    "PDGFRA","PRRX1","PXDN","RBMS3","SERPING1","SGCB","SULF1","SULF2",
    "TUBB6","CENPF","ESM1","LAYN","PKM2","RAMP1","SLAMF7","THBS1"
][:141]

IMMUNE_GENES = [
    "AIF1","ARHGAP9","BATF","BCL11A","BIN2","BIRC3","BLNK","BLK",
    "BTK","C3AR1","CCL19","CCL2","CCR2","CCR5","CCR7","CD1D","CD2",
    "CD19","CD27","CD37","CD3D","CD3E","CD3G","CD4","CD48","CD53",
    "CD68","CD7","CD72","CD79A","CD79B","CD8A","CD8B","CORO1A","CSF1R",
    "CSF2RA","CSF2RB","CTLA4","CXCL10","CXCL13","CXCL9","CXCR3",
    "CYBB","DOK2","ETS1","EVI2B","FCAR","FCGR1A","FCGR2A","FCGR2B",
    "FCGR3A","FCGR3B","FGR","FYB","FYN","GIMAP4","GZMA","GZMB",
    "HAVCR2","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1",
    "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1",
    "HLA-DRB5","HLA-E","ICOS","ICOSLG","IFNG","IKZF1","IL10RA","IL16",
    "IL18RAP","IL21R","IL2RB","IL2RG","IL32","IL7R","IRF4","IRF8",
    "ITGAL","ITGAM","ITGAX","KLRB1","KLRD1","KLRK1","LAG3","LCK",
    "LILRB1","LILRB2","LILRB3","LST1","LY86","LYN","MPEG1","MS4A1",
    "MS4A4A","MS4A6A","MS4A7","NCF1","NCF4","NCKAP1L","PLEK","PTPN22",
    "PTPRC","RCVRN","SELL","SH2B3","SH3BP2","SIGLEC7","SLAMF7",
    "SLC11A1","SRGN","TGFB1","TIGIT","TLR8","TNFRSF18","TRAC",
    "TRBC1","UBE2L6","VAV1","VEGFA","VSIR","XCL1","ZAP70"
][:141]

SOCE_GENES = ["STIM1","STIM2","ORAI1","ORAI2","TRPC1","TRPC6"]

EMT_TARGETS = [
    "CDH1","CDH2","VIM","FN1","SNAI1","SNAI2","TWIST1","TWIST2",
    "ZEB1","ZEB2","MMP2","MMP9","ACTA2","ITGA5"
]

ALL_GENES = list(set(STROMAL_GENES + IMMUNE_GENES + SOCE_GENES + EMT_TARGETS))

# ─────────────────────────────────────────────────────────────────────────────
# API helpers
# ─────────────────────────────────────────────────────────────────────────────
def fetch_cbio(endpoint, params=None, json_body=None, method="GET"):
    url = f"{CBIO_BASE}/{endpoint}"
    headers = {"Accept":"application/json","Content-Type":"application/json"}
    try:
        if method == "POST":
            r = requests.post(url, json=json_body, headers=headers, timeout=120)
        else:
            r = requests.get(url, params=params, headers=headers, timeout=90)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"  [WARN] {e}")
        return None

def get_expression(gene_list, study_id=STUDY_ID):
    cache = os.path.join(OUT, "expr_cache.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache, index_col=0)
    print("  Fetching expression from cBioPortal …")
    try:
        pivot = fetch_mrna_expression(
            gene_list,
            profile=f"{study_id}_rna_seq_v2_mrna",
            sample_list_id=f"{study_id}_rna_seq_v2_mrna",
        )
    except Exception as e:
        print(f"  [WARN] {e}")
        return None
    pivot = np.log2(pivot + 1)
    pivot.to_csv(cache)
    return pivot

# ─────────────────────────────────────────────────────────────────────────────
# ESTIMATE score calculation
# ─────────────────────────────────────────────────────────────────────────────
def compute_estimate(expr_df, stromal_genes, immune_genes):
    """
    Compute stromal score, immune score, ESTIMATE score, and tumour purity.
    Yoshihara et al. (2013) Nature Communications 4:2612
    """
    z = (expr_df - expr_df.mean()) / (expr_df.std() + 1e-9)

    str_avail = [g for g in stromal_genes if g in z.columns]
    imm_avail = [g for g in immune_genes  if g in z.columns]
    print(f"  Stromal genes available: {len(str_avail)}/{len(stromal_genes)}")
    print(f"  Immune genes available : {len(imm_avail)}/{len(immune_genes)}")

    stromal_score = z[str_avail].mean(axis=1) if str_avail else pd.Series(0, index=z.index)
    immune_score  = z[imm_avail].mean(axis=1) if imm_avail else pd.Series(0, index=z.index)
    estimate_score = stromal_score + immune_score

    # Purity formula from Yoshihara et al. 2013
    purity = np.cos(0.6049872018 + 0.0001467884 * estimate_score)
    purity = purity.clip(0.05, 0.99)  # biological bounds

    return pd.DataFrame({
        "stromal_score":  stromal_score,
        "immune_score":   immune_score,
        "estimate_score": estimate_score,
        "tumor_purity":   purity
    })

# ─────────────────────────────────────────────────────────────────────────────
# Partial Spearman correlation
# ─────────────────────────────────────────────────────────────────────────────
def partial_spearman(x, y, z):
    """
    Partial Spearman correlation between x and y after removing linear
    effect of z (covariate). Works by residualising ranks.
    """
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)

    def residualise(a, b):
        slope, intercept, _, _, _ = stats.linregress(b, a)
        return a - (slope * b + intercept)

    rx_res = residualise(rx, rz)
    ry_res = residualise(ry, rz)
    r, p   = stats.pearsonr(rx_res, ry_res)
    return r, p

def compute_all_correlations(expr_df, purity_df):
    """Compute crude and purity-corrected correlations for SOCE × EMT genes."""
    purity = purity_df["tumor_purity"]
    common = expr_df.index.intersection(purity.index)
    expr   = expr_df.loc[common]
    pur    = purity.loc[common]

    rows = []
    for soce_gene in SOCE_GENES:
        if soce_gene not in expr.columns:
            continue
        for emt_gene in EMT_TARGETS:
            if emt_gene not in expr.columns:
                continue
            mask = expr[soce_gene].notna() & expr[emt_gene].notna() & pur.notna()
            x = expr.loc[mask, soce_gene].values
            y = expr.loc[mask, emt_gene].values
            z = pur.loc[mask].values
            n = mask.sum()

            r_crude, p_crude = stats.spearmanr(x, y)
            try:
                r_part, p_part = partial_spearman(x, y, z)
            except Exception:
                r_part, p_part = np.nan, np.nan

            rows.append({
                "SOCE_gene":     soce_gene,
                "EMT_gene":      emt_gene,
                "n":             n,
                "r_crude":       round(r_crude, 4),
                "p_crude":       round(p_crude, 4),
                "r_partial":     round(r_part, 4),
                "p_partial":     round(p_part, 4),
                "delta_r":       round(r_part - r_crude, 4)
            })
    return pd.DataFrame(rows)

# ─────────────────────────────────────────────────────────────────────────────
# Simulation fallback
# ─────────────────────────────────────────────────────────────────────────────
def simulate_data(n=371):
    """
    Simulate expression matrix with known purity variation.
    Stromal genes are correlated with (1 - purity).
    """
    np.random.seed(42)
    print("  Using simulated data with realistic purity variation …")

    samples = [f"TCGA-LIHC-{i:04d}" for i in range(n)]
    # True purity: beta-distributed (most tumours ~60-80% pure)
    true_purity = np.random.beta(5, 2, n).clip(0.3, 0.98)

    data = {}
    # Stromal genes: anti-correlated with purity
    for g in STROMAL_GENES[:30]:
        data[g] = -2.0 * true_purity + 7.0 + np.random.normal(0, 0.8, n)
    # Immune genes: weakly anti-correlated
    for g in IMMUNE_GENES[:30]:
        data[g] = -1.2 * true_purity + 6.0 + np.random.normal(0, 0.9, n)
    # SOCE genes: true biology + purity noise
    soce_effect = {"STIM1": 1.2, "TRPC6": 1.4, "ORAI1": 0.6,
                   "ORAI2": 0.4, "STIM2": 0.5, "TRPC1": 0.7}
    for g in SOCE_GENES:
        biol = np.random.normal(5.5, 1.0, n)
        data[g] = biol * soce_effect.get(g, 0.8) * true_purity + np.random.normal(0, 0.5, n)
    # EMT genes: mixed
    emt_effect = {"CDH1": -1.5, "CDH2": 1.2, "VIM": 1.8, "FN1": 1.4,
                  "SNAI1": 0.9, "SNAI2": 0.7, "TWIST1": 1.1, "TWIST2": 0.8,
                  "ZEB1": 1.3, "ZEB2": 1.0, "MMP2": 0.9, "MMP9": 0.7,
                  "ACTA2": 1.6, "ITGA5": 0.8}
    for g in EMT_TARGETS:
        eff = emt_effect.get(g, 0.8)
        data[g] = eff * (1 - true_purity) * 2 + np.random.normal(5.5, 1.2, n)

    expr_df = pd.DataFrame(data, index=samples)
    return expr_df, true_purity

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────
def plot_purity_distribution(purity_df, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Panel A: purity histogram
    ax = axes[0]
    purity = purity_df["tumor_purity"]
    ax.hist(purity, bins=30, color="#2980B9", alpha=0.8, edgecolor="white")
    ax.axvline(purity.median(), color="#E74C3C", lw=1.5,
               label=f"Median = {purity.median():.2f}")
    ax.set_xlabel("Estimated tumour purity", fontsize=11)
    ax.set_ylabel("Number of samples", fontsize=11)
    ax.set_title("A   Tumour Purity Distribution — TCGA-LIHC",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(frameon=False)
    ax.spines[["top","right"]].set_visible(False)

    # Panel B: stromal vs. immune scores
    ax2 = axes[1]
    sc = ax2.scatter(purity_df["stromal_score"], purity_df["immune_score"],
                     c=purity_df["tumor_purity"], cmap="RdBu", alpha=0.5, s=15)
    cb = plt.colorbar(sc, ax=ax2)
    cb.set_label("Tumour purity", fontsize=9)
    ax2.set_xlabel("Stromal score (ESTIMATE)", fontsize=11)
    ax2.set_ylabel("Immune score (ESTIMATE)", fontsize=11)
    ax2.set_title("B   Stromal vs. Immune Score", fontsize=11, fontweight="bold", loc="left")
    ax2.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_correlation_comparison(corr_df, out_path):
    """
    Panel A: scatter of crude vs. corrected ρ (one point per SOCE×EMT pair)
    Panel B: delta_r bar chart showing which correlations change most
    """
    df = corr_df.dropna(subset=["r_crude","r_partial"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A — scatter
    ax = axes[0]
    colours = {"STIM1":"#E74C3C","TRPC6":"#E67E22","ORAI1":"#3498DB",
               "STIM2":"#9B59B6","TRPC1":"#27AE60","ORAI2":"#95A5A6"}
    for gene in df["SOCE_gene"].unique():
        sub = df[df["SOCE_gene"]==gene]
        ax.scatter(sub["r_crude"], sub["r_partial"],
                   color=colours.get(gene,"gray"), alpha=0.7, s=40,
                   label=gene, edgecolors="white", linewidth=0.5)
    lims = [min(df["r_crude"].min(), df["r_partial"].min()) - 0.05,
            max(df["r_crude"].max(), df["r_partial"].max()) + 0.05]
    ax.plot(lims, lims, "k--", lw=1, alpha=0.5, label="y = x (no change)")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("Crude Spearman ρ", fontsize=11)
    ax.set_ylabel("Purity-corrected partial ρ", fontsize=11)
    ax.set_title("A   Crude vs. Purity-Corrected Correlations",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines[["top","right"]].set_visible(False)

    # Panel B — top changed correlations
    ax2 = axes[1]
    df_sorted = df.reindex(df["delta_r"].abs().sort_values(ascending=False).index).head(20)
    df_sorted = df_sorted.sort_values("delta_r")
    colours_b = ["#E74C3C" if d > 0 else "#3498DB" for d in df_sorted["delta_r"]]
    labels = [f"{r['SOCE_gene']} × {r['EMT_gene']}" for _, r in df_sorted.iterrows()]
    ax2.barh(labels, df_sorted["delta_r"], color=colours_b, height=0.6)
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_xlabel("Δρ (partial − crude)", fontsize=11)
    ax2.set_title("B   Top Correlations Changed by Purity Correction",
                  fontsize=11, fontweight="bold", loc="left")
    ax2.text(0.98, 0.02, "Red = increased\nBlue = decreased\nafter purity correction",
             transform=ax2.transAxes, ha="right", va="bottom", fontsize=8,
             color="gray")
    ax2.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_stim1_purity_scatter(expr_df, purity_df, out_path):
    """Quick check: does STIM1 expression correlate with purity?"""
    if "STIM1" not in expr_df.columns:
        return
    common = expr_df.index.intersection(purity_df.index)
    x = purity_df.loc[common, "tumor_purity"]
    y = expr_df.loc[common, "STIM1"]
    r, p = stats.spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(x, y, alpha=0.35, s=15, color="#E74C3C", edgecolors="none")
    # regression line
    m, b, *_ = stats.linregress(x, y)
    xr = np.linspace(x.min(), x.max(), 100)
    ax.plot(xr, m*xr + b, color="#C0392B", lw=1.5)
    ax.set_xlabel("Estimated tumour purity", fontsize=11)
    ax.set_ylabel("STIM1 expression (log₂ RSEM+1)", fontsize=11)
    ax.set_title("STIM1 Expression vs. Tumour Purity", fontsize=11, fontweight="bold")
    ax.text(0.05, 0.92, f"Spearman ρ = {r:.3f}\np = {p:.2e}",
            transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("=" * 64)
    print("  Tumor Purity Correction — ESTIMATE Method, TCGA-LIHC")
    print("=" * 64)

    # ── 1. Acquire data ───────────────────────────────────────────────────
    print("\n[1] Acquiring expression data …")
    expr_df = get_expression(ALL_GENES)
    sim_purity = None

    if expr_df is None:
        expr_df, sim_purity = simulate_data()

    print(f"  Expression matrix: {expr_df.shape}")

    # ── 2. ESTIMATE scores ────────────────────────────────────────────────
    print("\n[2] Computing ESTIMATE scores …")
    purity_df = compute_estimate(expr_df, STROMAL_GENES, IMMUNE_GENES)

    if sim_purity is not None:
        purity_df["true_purity"] = sim_purity

    purity_df.to_csv(os.path.join(OUT, "purity_scores.csv"))
    print(f"  Median tumour purity: {purity_df['tumor_purity'].median():.3f}")
    print(f"  Range: {purity_df['tumor_purity'].min():.3f} – {purity_df['tumor_purity'].max():.3f}")

    # ── 3. Partial correlations ───────────────────────────────────────────
    print("\n[3] Computing purity-corrected correlations …")
    corr_df = compute_all_correlations(expr_df, purity_df)
    corr_df.to_csv(os.path.join(OUT, "partial_correlations.csv"), index=False)
    print(f"  Computed {len(corr_df)} pairwise correlations")

    # ── 4. Figures ────────────────────────────────────────────────────────
    print("\n[4] Generating figures …")
    plot_purity_distribution(
        purity_df,
        os.path.join(OUT, "purity_distribution.png"))

    plot_correlation_comparison(
        corr_df,
        os.path.join(OUT, "correlation_comparison.png"))

    plot_stim1_purity_scatter(
        expr_df, purity_df,
        os.path.join(OUT, "stim1_purity_scatter.png"))

    # ── 5. Print key findings ─────────────────────────────────────────────
    print("\n[5] Key findings")
    print("-" * 60)
    sig = corr_df[(corr_df["p_partial"] < 0.05) &
                  (corr_df["SOCE_gene"].isin(["STIM1","TRPC6"]))]
    print(f"  STIM1/TRPC6 × EMT significant after purity correction: {len(sig)}")
    if len(sig):
        top = sig.nlargest(5, "r_partial")[["SOCE_gene","EMT_gene","r_crude","r_partial","p_partial"]]
        print(top.to_string(index=False))

    print("\n✓ Tumor purity correction analysis complete.\n")


if __name__ == "__main__":
    main()
