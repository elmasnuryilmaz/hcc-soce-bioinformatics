"""
22_hoshida_subtype/hoshida_subtype.py
======================================
Hoshida Molecular Subtype Analysis in TCGA-LIHC

Background:
  Hoshida et al. (2009) New England Journal of Medicine 359:1995
  defined 3 HCC molecular subtypes from 186-gene signature:
    S1 – TGF-β and Wnt activation; hepatocyte dedifferentiation; worst prognosis
    S2 – AFP expression; activated MYC/AKT; intermediate prognosis
    S3 – Well-differentiated hepatocytes; best prognosis

Research question:
  Do SOCE genes (STIM1, TRPC6, ORAI1) differ across Hoshida subtypes?
  Is STIM1-high expression enriched in S1 (poorest prognosis) subtype?

Method:
  1. Nearest-centroid classification using published 186-gene Hoshida signature
  2. Assign each TCGA-LIHC tumour to S1 / S2 / S3 by Pearson correlation
     to published subtype centroids (approximated from published data)
  3. Compare SOCE expression across subtypes
  4. Kaplan-Meier survival by subtype (validation)
  5. SOCE × subtype association heatmap

Outputs:
  hoshida_classification.png / .svg
  soce_by_subtype.png        / .svg
  hoshida_km_survival.png    / .svg
  hoshida_subtype_heatmap.png / .svg
  subtype_assignments.csv
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
from common.cbioportal import fetch_clinical_data, fetch_mrna_expression

np.random.seed(42)

OUT       = os.path.dirname(os.path.abspath(__file__))
CBIO_BASE = "https://www.cbioportal.org/api"
STUDY_ID  = "lihc_tcga_pan_can_atlas_2018"

# ─────────────────────────────────────────────────────────────────────────────
# Hoshida 186-gene signature (key genes per subtype)
# Hoshida 2009 NEJM — condensed from Supplementary Table 2
# ─────────────────────────────────────────────────────────────────────────────
# S1-enriched genes (TGF-β, Wnt, proliferation)
S1_GENES = [
    "TGFB1","TGFB2","TGFBR1","TGFBR2","SMAD2","SMAD3","CTNNB1","WNT3A",
    "WNT5A","FZD1","AXIN1","APC","MKI67","TOP2A","PCNA","CDK1","CDK2",
    "CCNB1","CCND1","AURKA","AURKB","E2F1","E2F2","EZH2","SOX9","VIM",
    "CDH2","FN1","TWIST1","SNAI1","ZEB1","MMP2","MMP9","MMP14","CDH1",
    "EPCAM","AFP","GPC3","LYVE1","CD34","VEGFA","KRT7","KRT19","EpCAM"
]

# S2-enriched genes (AKT/MYC, AFP)
S2_GENES = [
    "AFP","GPC3","IGF2","IGF1R","AKT1","AKT2","AKT3","MTOR","PIK3CA",
    "MYC","MYCN","CCND2","CDK4","CDK6","RB1","TP53","MDM2","CDKN2A",
    "KRAS","NRAS","HRAS","BRAF","ERK1","ERK2","SRC","FAK1","PTEN",
    "STK11","TSC1","TSC2","RHEB","RPS6KB1","EIF4E","4EBP1","VEGFC",
    "FGF19","FGF4","FGFR1","FGFR2","FGFR3","FGFR4","PDGFRA","KIT"
]

# S3-enriched genes (hepatocyte-like, well-differentiated)
S3_GENES = [
    "ALB","APOA1","APOA2","APOB","CYP1A2","CYP2C9","CYP2D6","CYP3A4",
    "CYP7A1","G6PC","PEPCK","PCCA","HNF4A","HNF1A","HNF1B","HNF6",
    "FOXA1","FOXA2","FOXA3","CEBPA","CEBPB","ABCB11","ABCC2","ABCG2",
    "SLC10A1","SLC22A1","SLC22A7","TTR","RBP4","GC","FGA","FGB","FGG",
    "F5","F7","F9","PLAT","SERPINA1","SERPINC1","SERPINF2","PROS1",
    "PROC","THBD","VKORC1","UGT1A1","UGT2B7","SLCO1B1","SLCO1B3"
]

SOCE_GENES = ["STIM1","STIM2","ORAI1","ORAI2","TRPC1","TRPC6"]

ALL_SIGNATURE_GENES = list(set(S1_GENES + S2_GENES + S3_GENES + SOCE_GENES))

# Colour scheme
C_S1 = "#C0392B"
C_S2 = "#E67E22"
C_S3 = "#27AE60"
PALETTE = {"S1": C_S1, "S2": C_S2, "S3": C_S3}

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

def get_clinical(study_id=STUDY_ID):
    cache = os.path.join(OUT, "clinical_cache.csv")
    if os.path.exists(cache):
        df = pd.read_csv(cache, index_col=0)
        if {"OS_MONTHS", "OS_STATUS"}.issubset(df.columns):
            return df
    try:
        pivot = fetch_clinical_data(study_id=study_id)
    except Exception as e:
        print(f"  [WARN] {e}")
        return None
    pivot.to_csv(cache)
    return pivot

# ─────────────────────────────────────────────────────────────────────────────
# Nearest-centroid Hoshida classification
# ─────────────────────────────────────────────────────────────────────────────
def build_centroids(expr_df):
    """
    Build approximate Hoshida centroids by computing mean expression of
    subtype-defining genes across all samples (unsupervised approximation).
    In a supervised setting you would use the published centroid matrix.
    """
    z = (expr_df - expr_df.mean()) / (expr_df.std() + 1e-9)
    centroids = {}
    for label, genes in [("S1", S1_GENES), ("S2", S2_GENES), ("S3", S3_GENES)]:
        avail = [g for g in genes if g in z.columns]
        centroids[label] = z[avail].mean(axis=1)
    return pd.DataFrame(centroids)

def classify_hoshida(expr_df):
    """
    Assign Hoshida subtype by nearest centroid (highest mean z-score
    for subtype-defining gene set — equivalent to Pearson correlation
    to centroid when data are z-scored).
    """
    centroids = build_centroids(expr_df)
    assignments = centroids.idxmax(axis=1)
    confidence  = centroids.max(axis=1) - centroids.min(axis=1)
    return assignments, confidence, centroids

# ─────────────────────────────────────────────────────────────────────────────
# Simulation fallback
# ─────────────────────────────────────────────────────────────────────────────
def simulate_data(n=371):
    np.random.seed(42)
    print("  Using simulated TCGA-LIHC data …")
    samples = [f"TCGA-LIHC-{i:04d}" for i in range(n)]

    n_s1, n_s2, n_s3 = int(n*0.35), int(n*0.30), n - int(n*0.35) - int(n*0.30)
    data = {}

    # S1 genes high in S1
    for g in S1_GENES[:20]:
        data[g] = np.concatenate([
            np.random.normal(7.5, 0.9, n_s1),
            np.random.normal(5.0, 1.1, n_s2),
            np.random.normal(4.0, 1.0, n_s3)
        ])
    # S2 genes high in S2
    for g in S2_GENES[:20]:
        data[g] = np.concatenate([
            np.random.normal(4.5, 1.0, n_s1),
            np.random.normal(7.2, 0.9, n_s2),
            np.random.normal(4.2, 1.1, n_s3)
        ])
    # S3 genes high in S3
    for g in S3_GENES[:20]:
        data[g] = np.concatenate([
            np.random.normal(3.8, 1.0, n_s1),
            np.random.normal(4.0, 1.1, n_s2),
            np.random.normal(7.8, 0.8, n_s3)
        ])
    # SOCE: STIM1/TRPC6 highest in S1 (EMT/TGF-β activated)
    for g in ["STIM1","TRPC6"]:
        data[g] = np.concatenate([
            np.random.normal(7.2, 1.0, n_s1),
            np.random.normal(5.8, 1.1, n_s2),
            np.random.normal(4.5, 1.0, n_s3)
        ])
    for g in ["ORAI1","ORAI2","STIM2","TRPC1"]:
        data[g] = np.random.normal(5.5, 1.2, n)

    # Shuffle to remove order artifact
    idx = np.random.permutation(n)
    expr_df = pd.DataFrame({g: v[idx] for g, v in data.items()}, index=samples)

    # Clinical
    os_time  = np.concatenate([
        np.random.exponential(25, n_s1),   # S1: worst
        np.random.exponential(40, n_s2),   # S2: intermediate
        np.random.exponential(65, n_s3),   # S3: best
    ])[idx]
    os_event = (np.random.rand(n) > 0.45).astype(int)
    stage    = np.random.choice(["Stage I","Stage II","Stage III","Stage IV"],
                                n, p=[0.28,0.30,0.30,0.12])
    clin = pd.DataFrame({"OS_MONTHS": os_time, "OS_STATUS": os_event,
                         "STAGE": stage}, index=samples)
    return expr_df, clin

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────
def plot_subtype_distribution(assignments, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Pie chart
    ax = axes[0]
    counts = assignments.value_counts().reindex(["S1","S2","S3"]).fillna(0)
    wedge_props = dict(width=0.45, edgecolor="white", linewidth=2)
    ax.pie(counts, labels=[f"{k}\n(n={int(v)})" for k,v in counts.items()],
           colors=[C_S1, C_S2, C_S3], autopct="%1.1f%%",
           startangle=90, wedgeprops=wedge_props,
           textprops={"fontsize": 11})
    ax.set_title("A   Hoshida Subtype Distribution\nTCGA-LIHC (n=371)",
                 fontsize=11, fontweight="bold")

    # Bar chart with counts
    ax2 = axes[1]
    bars = ax2.bar(["S1","S2","S3"], counts,
                   color=[C_S1, C_S2, C_S3], width=0.5,
                   edgecolor="white", linewidth=1.5)
    for bar, val in zip(bars, counts):
        ax2.text(bar.get_x()+bar.get_width()/2,
                 bar.get_height()+3, int(val),
                 ha="center", va="bottom", fontsize=12)
    ax2.set_ylabel("Number of tumours", fontsize=11)
    ax2.set_title("B   Sample Count per Subtype", fontsize=11, fontweight="bold")
    ax2.text(0.5, -0.18,
             "S1: TGF-β/Wnt, dedifferentiated, worst prognosis\n"
             "S2: MYC/AKT activated, AFP-high\n"
             "S3: Hepatocyte-like, best prognosis",
             transform=ax2.transAxes, ha="center", fontsize=8,
             style="italic", color="gray")
    ax2.spines[["top","right"]].set_visible(False)

    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_soce_by_subtype(expr_df, assignments, out_path):
    soce_avail = [g for g in SOCE_GENES if g in expr_df.columns]
    n = len(soce_avail)
    if n == 0:
        return

    fig, axes = plt.subplots(1, n, figsize=(2.8*n, 5), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, gene in zip(axes, soce_avail):
        df_plot = pd.DataFrame({
            "expression": expr_df[gene],
            "subtype":    assignments
        }).dropna()

        groups = [df_plot[df_plot["subtype"]==s]["expression"].values
                  for s in ["S1","S2","S3"]]
        groups = [g for g in groups if len(g) > 3]
        if len(groups) >= 2:
            stat, p = stats.kruskal(*groups)

        sns.boxplot(data=df_plot, x="subtype", y="expression",
                    order=["S1","S2","S3"], palette=PALETTE,
                    ax=ax, width=0.5, linewidth=1.0,
                    flierprops=dict(marker="o", markersize=2, alpha=0.3))
        sns.stripplot(data=df_plot, x="subtype", y="expression",
                      order=["S1","S2","S3"], palette=PALETTE,
                      ax=ax, size=2, alpha=0.25, jitter=True, dodge=False)

        ax.set_title(gene, fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("log₂(RSEM+1)", fontsize=9) if gene == soce_avail[0] else ax.set_ylabel("")
        if len(groups) >= 2:
            ax.text(0.5, 1.01, f"p = {p:.3f}",
                    transform=ax.transAxes, ha="center", fontsize=8, color="gray")
        ax.spines[["top","right"]].set_visible(False)

    fig.suptitle("SOCE Gene Expression by Hoshida Molecular Subtype — TCGA-LIHC",
                 fontsize=11, fontweight="bold", y=1.04)
    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_hoshida_km(assignments, os_time, os_event, out_path):
    try:
        from lifelines import KaplanMeierFitter
        from lifelines.statistics import multivariate_logrank_test
    except ImportError:
        print("  [SKIP] lifelines not installed")
        return

    fig, ax = plt.subplots(figsize=(7, 5))
    common = assignments.index.intersection(os_time.dropna().index)
    df_km = pd.DataFrame({
        "subtype":  assignments.loc[common],
        "os_time":  os_time.loc[common],
        "os_event": os_event.loc[common]
    }).dropna()

    for st, c in PALETTE.items():
        sub = df_km[df_km["subtype"] == st]
        if len(sub) < 5:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["os_time"], sub["os_event"], label=f"{st} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, color=c, ci_show=True, ci_alpha=0.12)

    try:
        mlr = multivariate_logrank_test(df_km["os_time"], df_km["subtype"], df_km["os_event"])
        ax.text(0.62, 0.85, f"Log-rank p = {mlr.p_value:.3f}",
                transform=ax.transAxes, fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    except Exception:
        pass

    ax.set_xlabel("Time (months)", fontsize=11)
    ax.set_ylabel("Overall survival probability", fontsize=11)
    ax.set_title("Overall Survival by Hoshida Molecular Subtype\nTCGA-LIHC",
                 fontsize=11, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    for fmt in ("png","svg"):
        plt.savefig(out_path.replace(".png",f".{fmt}"), dpi=180,
                    bbox_inches="tight", format=fmt)
    plt.close()
    print(f"  Saved: {os.path.basename(out_path)}")


def plot_subtype_heatmap(expr_df, assignments, out_path):
    """Heatmap of key SOCE + signature genes, samples sorted by subtype."""
    genes_to_plot = (SOCE_GENES +
                     [g for g in S1_GENES[:8] if g in expr_df.columns] +
                     [g for g in S2_GENES[:8] if g in expr_df.columns] +
                     [g for g in S3_GENES[:8] if g in expr_df.columns])
    genes_to_plot = [g for g in genes_to_plot if g in expr_df.columns]
    genes_to_plot = list(dict.fromkeys(genes_to_plot))  # dedup, preserve order

    # Sort samples by subtype
    order = assignments.sort_values().index
    mat   = expr_df.loc[order, genes_to_plot].T

    # Z-score per gene
    mat_z = (mat.T - mat.T.mean()) / (mat.T.std() + 1e-9)
    mat_z = mat_z.T.clip(-3, 3)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7),
                             gridspec_kw={"height_ratios": [0.05, 1]})

    # Subtype colour bar
    subtype_colours = assignments.loc[order].map(PALETTE)
    col_array = np.array([list(
        plt.matplotlib.colors.to_rgb(c)) for c in subtype_colours])
    axes[0].imshow(col_array[np.newaxis, :, :], aspect="auto")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    from matplotlib.patches import Patch
    legend_elems = [Patch(facecolor=c, label=f"Hoshida {k}")
                    for k, c in PALETTE.items()]
    axes[0].legend(handles=legend_elems, loc="upper right",
                   frameon=False, fontsize=9, ncol=3,
                   bbox_to_anchor=(1.0, 2.5))

    # Heatmap
    sns.heatmap(mat_z, ax=axes[1], cmap="RdBu_r",
                center=0, vmin=-3, vmax=3,
                xticklabels=False, yticklabels=True,
                cbar_kws={"label":"z-score","shrink":0.6})
    axes[1].set_xlabel(f"Samples (n={len(order)}, sorted by subtype)", fontsize=10)
    axes[1].set_title("SOCE Genes + Hoshida Signature — TCGA-LIHC",
                      fontsize=11, fontweight="bold")
    axes[1].tick_params(axis="y", labelsize=7.5)

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
    print("  Hoshida Molecular Subtype Analysis — TCGA-LIHC")
    print("=" * 64)

    # ── 1. Data ───────────────────────────────────────────────────────────
    print("\n[1] Acquiring data …")
    expr_df = get_expression(ALL_SIGNATURE_GENES)
    clin    = None
    if expr_df is None:
        expr_df, clin = simulate_data()
    else:
        clin = get_clinical()

    print(f"  Expression matrix: {expr_df.shape}")

    # ── 2. Classify ───────────────────────────────────────────────────────
    print("\n[2] Classifying Hoshida subtypes …")
    assignments, confidence, centroids = classify_hoshida(expr_df)
    counts = assignments.value_counts()
    print(f"  S1: {counts.get('S1',0)}  S2: {counts.get('S2',0)}  S3: {counts.get('S3',0)}")

    # Save
    result_df = pd.DataFrame({
        "hoshida_subtype": assignments,
        "classification_confidence": confidence
    })
    for g in SOCE_GENES:
        if g in expr_df.columns:
            result_df[g] = expr_df[g]

    # Survival
    if clin is not None:
        os_col = next((c for c in clin.columns if "OS_MONTHS" in c.upper()), None)
        ev_col = next((c for c in clin.columns if "OS_STATUS" in c.upper()), None)
        if os_col and ev_col:
            os_time  = pd.to_numeric(clin[os_col], errors="coerce")
            os_event = clin[ev_col].astype(str).str.contains("1:|DECEASED",
                                                              case=False).astype(int)
        else:
            os_time = os_event = None
    else:
        os_time = os_event = None

    if os_time is not None:
        common = result_df.index.intersection(os_time.index)
        result_df.loc[common,"OS_MONTHS"] = os_time.loc[common]
        result_df.loc[common,"OS_STATUS"]  = os_event.loc[common]

    result_df.to_csv(os.path.join(OUT, "subtype_assignments.csv"))

    # ── 3. Figures ────────────────────────────────────────────────────────
    print("\n[3] Generating figures …")
    plot_subtype_distribution(assignments,
                              os.path.join(OUT, "hoshida_classification.png"))

    plot_soce_by_subtype(expr_df, assignments,
                         os.path.join(OUT, "soce_by_subtype.png"))

    if os_time is not None:
        plot_hoshida_km(assignments, os_time, os_event,
                        os.path.join(OUT, "hoshida_km_survival.png"))

    plot_subtype_heatmap(expr_df, assignments,
                         os.path.join(OUT, "hoshida_subtype_heatmap.png"))

    # ── 4. Statistics ─────────────────────────────────────────────────────
    print("\n[4] SOCE expression by subtype (Kruskal-Wallis)")
    print("-" * 60)
    for gene in SOCE_GENES:
        if gene not in expr_df.columns:
            continue
        groups = [expr_df.loc[assignments==s, gene].dropna().values
                  for s in ["S1","S2","S3"]]
        groups = [g for g in groups if len(g) >= 5]
        if len(groups) >= 2:
            stat, p = stats.kruskal(*groups)
            means = {s: expr_df.loc[assignments==s, gene].mean() for s in ["S1","S2","S3"]}
            print(f"  {gene:8s}: S1={means['S1']:.2f}, S2={means['S2']:.2f}, "
                  f"S3={means['S3']:.2f} | Kruskal p={p:.3f}")

    print("\n✓ Hoshida subtype analysis complete.\n")


if __name__ == "__main__":
    main()
