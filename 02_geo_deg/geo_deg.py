#!/usr/bin/env python3
"""
02_geo_deg.py  –  Differential Expression Analysis (GSE140202)
==============================================================
Identifies differentially expressed genes between sorafenib-resistant
and sorafenib-sensitive Huh7 hepatocellular carcinoma cells.

Dataset
-------
  GSE140202 (Roessler et al. / GEO)
  Platform : Affymetrix HG-U133A 2.0 (GPL571)
  Samples  : sorafenib-resistant Huh7  vs.  parental Huh7
  Genes    : ~12,998 unique gene symbols after filtering

Method
------
  1. Download series matrix from NCBI GEO via GEOparse
  2. Log2-transform expression values (if not already)
  3. Compute fold-change (resistant / sensitive) and t-test p-value
  4. Apply BH correction for FDR
  5. Volcano plot (|log2FC| > 1 and FDR < 0.05 highlighted)

Outputs
-------
  deg_full_results.csv       – All ~12,998 genes
  deg_significant.csv        – FDR < 0.05 subset
  volcano_plot.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn GEOparse
  python geo_deg.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import GEOparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# ── Config ──────────────────────────────────────────────────────────────────
OUTDIR   = os.path.dirname(os.path.abspath(__file__))
GEO_ID   = "GSE140202"
FC_THRESH = 1.0     # log2 fold-change threshold for highlighting
FDR_THRESH = 0.05   # FDR threshold

# Manually define sample group labels based on GSE140202 metadata.
# Group labels: 'resistant' vs 'sensitive'
# Update RESISTANT_PATTERN / SENSITIVE_PATTERN to match GSM titles if needed.
RESISTANT_PATTERN = "resistant"
SENSITIVE_PATTERN = "parental"

# ── 1. Download GEO dataset ──────────────────────────────────────────────────
def load_geo(geo_id: str) -> tuple[pd.DataFrame, dict]:
    """
    Downloads GSE matrix via GEOparse (cached in geo_cache/).
    Returns (expression_df, sample_metadata_dict).
    expression_df: rows = probes, columns = GSM IDs
    """
    cache_dir = os.path.join(OUTDIR, "geo_cache")
    os.makedirs(cache_dir, exist_ok=True)

    print(f"  Downloading {geo_id} from NCBI GEO …")
    gse = GEOparse.get_GEO(geo_id, destdir=cache_dir, silent=True)

    # Build expression matrix
    pivot_tables = []
    for gsm_id, gsm in gse.gsms.items():
        if gsm.table is not None and not gsm.table.empty:
            col = gsm.table.set_index("ID_REF")["VALUE"].rename(gsm_id)
            pivot_tables.append(col)

    if not pivot_tables:
        raise RuntimeError(
            f"{geo_id} GEO series matrix contains no sample expression tables. "
            "Use the checked-in DEG table or provide the original gene-level expression matrix."
        )

    expr = pd.concat(pivot_tables, axis=1)
    expr = expr.apply(pd.to_numeric, errors="coerce")

    # Sample metadata
    meta = {gsm_id: gsm.metadata for gsm_id, gsm in gse.gsms.items()}
    return expr, meta


# ── 2. Assign groups ─────────────────────────────────────────────────────────
def assign_groups(meta: dict) -> dict[str, str]:
    """
    Returns {gsm_id: 'resistant' | 'sensitive'} based on sample title.
    """
    groups = {}
    for gsm_id, m in meta.items():
        title = " ".join(m.get("title", [""])).lower()
        if RESISTANT_PATTERN in title:
            groups[gsm_id] = "resistant"
        elif SENSITIVE_PATTERN in title:
            groups[gsm_id] = "sensitive"
    return groups


# ── 3. Probe → gene symbol mapping ──────────────────────────────────────────
def map_probes_to_genes(gse) -> pd.Series:
    """Returns Series: probe_id → gene_symbol from the platform table."""
    for _, gpl in gse.gpls.items():
        tbl = gpl.table
        # Common column names for gene symbol in Affy platforms
        sym_col = next(
            (c for c in tbl.columns if "gene_symbol" in c.lower() or c.upper() == "GENE_SYMBOL"),
            None
        )
        if sym_col is None:
            sym_col = next(
                (c for c in tbl.columns if "symbol" in c.lower()), None
            )
        if sym_col and "ID" in tbl.columns:
            return tbl.set_index("ID")[sym_col].squeeze()
    return pd.Series(dtype=str)


# ── 4. Differential expression ───────────────────────────────────────────────
def run_deg(expr: pd.DataFrame, groups: dict) -> pd.DataFrame:
    """Welch t-test + BH correction on probe-level expression."""
    res_samples = [s for s, g in groups.items() if g == "resistant" and s in expr.columns]
    sen_samples = [s for s, g in groups.items() if g == "sensitive" and s in expr.columns]

    print(f"  Resistant samples: {len(res_samples)}")
    print(f"  Sensitive samples: {len(sen_samples)}")

    res_expr = expr[res_samples]
    sen_expr = expr[sen_samples]

    # Log2-transform if values look like linear counts (max > 100)
    if res_expr.values.max() > 100:
        res_expr = np.log2(res_expr + 1)
        sen_expr = np.log2(sen_expr + 1)

    mean_res = res_expr.mean(axis=1)
    mean_sen = sen_expr.mean(axis=1)
    log2fc   = mean_res - mean_sen

    # Welch t-test (unequal variance)
    t_stats, pvals = stats.ttest_ind(
        res_expr.values, sen_expr.values, axis=1, equal_var=False
    )

    # FDR correction (Benjamini-Hochberg)
    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")

    df_deg = pd.DataFrame({
        "log2FC":  log2fc.values,
        "mean_resistant": mean_res.values,
        "mean_sensitive": mean_sen.values,
        "t_stat":  t_stats,
        "p_value": pvals,
        "FDR":     fdr,
    }, index=expr.index)

    return df_deg


# ── 5. Aggregate to gene level ────────────────────────────────────────────────
def aggregate_to_gene(df_deg: pd.DataFrame, probe_map: pd.Series) -> pd.DataFrame:
    """Keep the probe with the largest |log2FC| per gene symbol."""
    df = df_deg.copy()
    df["gene"] = probe_map.reindex(df.index)
    df = df.dropna(subset=["gene"])
    df = df[df["gene"].str.strip() != ""]
    # Keep one probe per gene (highest |log2FC|)
    df["abs_fc"] = df["log2FC"].abs()
    df = df.sort_values("abs_fc", ascending=False).groupby("gene").first()
    df = df.drop(columns=["abs_fc"])
    df = df.sort_values("FDR")
    return df


def normalise_deg_schema(df: pd.DataFrame) -> pd.DataFrame:
    """Accept either this module's schema or the checked-in limma-style DEG table."""
    df = df.copy()
    if "gene" not in df.columns:
        df["gene"] = df.index.astype(str)
    if "log2FC" not in df.columns and "logFC" in df.columns:
        df["log2FC"] = df["logFC"]
    if "p_value" not in df.columns and "P.Value" in df.columns:
        df["p_value"] = df["P.Value"]
    if "FDR" not in df.columns and "adj.P.Val" in df.columns:
        df["FDR"] = df["adj.P.Val"]
    if "t_stat" not in df.columns and "t" in df.columns:
        df["t_stat"] = df["t"]
    if "gene" in df.columns:
        df["gene"] = df["gene"].astype(str)
        df = df.set_index("gene", drop=False)
    required = {"log2FC", "p_value", "FDR"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DEG table is missing required columns: {sorted(missing)}")
    return df


def load_existing_deg_table() -> pd.DataFrame:
    """Load the repository's DEG table when GEO no longer exposes matrix values."""
    deg_file = os.path.join(OUTDIR, "deg_full_results.csv")
    if not os.path.exists(deg_file):
        raise FileNotFoundError(f"Fallback DEG table not found: {deg_file}")
    print(f"  [fallback] Loading existing DEG table: {deg_file}")
    return normalise_deg_schema(pd.read_csv(deg_file))


# ── 6. Volcano plot ───────────────────────────────────────────────────────────
def volcano_plot(df: pd.DataFrame):
    """Publication-quality volcano plot."""
    df = normalise_deg_schema(df)
    df = df.dropna(subset=["p_value", "log2FC"])
    neg_log10p = -np.log10(df["p_value"].clip(lower=1e-300))

    colours = np.where(
        (df["FDR"] < FDR_THRESH) & (df["log2FC"] >  FC_THRESH), "#c0392b",
        np.where(
            (df["FDR"] < FDR_THRESH) & (df["log2FC"] < -FC_THRESH), "#2980b9",
            "#bdc3c7"
        )
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(df["log2FC"], neg_log10p, c=colours, s=10, alpha=0.7, linewidths=0)

    # Threshold lines
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8)
    ax.axvline( FC_THRESH,  color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(-FC_THRESH,  color="grey", linestyle="--", linewidth=0.8)

    # Label top DE genes
    soce_genes = ["STIM1", "TRPC6", "TRPC1", "ORAI1",
                  "VIM", "ZEB1", "ZEB2", "SNAI1", "TWIST1",
                  "RELA", "IL6R", "TNF", "ABCC2"]
    for gene in soce_genes:
        if gene in df.index:
            x = df.loc[gene, "log2FC"]
            y = -np.log10(df.loc[gene, "p_value"])
            ax.annotate(gene, (x, y), fontsize=7, ha="center",
                        xytext=(5, 3), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", lw=0.5))

    n_up   = ((df["FDR"] < FDR_THRESH) & (df["log2FC"] >  FC_THRESH)).sum()
    n_down = ((df["FDR"] < FDR_THRESH) & (df["log2FC"] < -FC_THRESH)).sum()
    ax.text(0.02, 0.98, f"↑ {n_up} up", transform=ax.transAxes,
            color="#c0392b", fontsize=10, va="top")
    ax.text(0.98, 0.98, f"↓ {n_down} down", transform=ax.transAxes,
            color="#2980b9", fontsize=10, va="top", ha="right")

    ax.set_xlabel("log₂ Fold-Change (Resistant / Sensitive)", fontsize=12)
    ax.set_ylabel("−log₁₀(p-value)", fontsize=12)
    ax.set_title(
        "Differential Expression: Sorafenib-Resistant vs. Sensitive Huh7\n"
        f"(GSE140202, FDR < {FDR_THRESH}, |log₂FC| > {FC_THRESH})",
        fontsize=12, fontweight="bold"
    )
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"volcano_plot.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved volcano_plot.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 02 — GEO Differential Expression (GSE140202)")
    print("=" * 60)

    try:
        expr, meta = load_geo(GEO_ID)
    except Exception as exc:
        print(f"  ⚠  GEO expression matrix unavailable: {exc}")
        df_gene = load_existing_deg_table()
        df_sig = df_gene[(df_gene["FDR"] < FDR_THRESH) & (df_gene["log2FC"].abs() > FC_THRESH)]
        df_sig.to_csv(os.path.join(OUTDIR, "deg_significant.csv"), index=False)
        print(f"  Significant (FDR<0.05, |log2FC|>1) : {len(df_sig)}")
        print("\n[3/3] Volcano plot …")
        volcano_plot(df_gene)
        print("\n  Done using existing DEG table. All outputs saved to:", OUTDIR)
        import sys; sys.exit(0)

    print(f"  Expression matrix: {expr.shape[0]} probes × {expr.shape[1]} samples")

    groups = assign_groups(meta)
    if len(groups) == 0:
        print("  ⚠  No samples matched resistant/sensitive patterns.")
        df_gene = load_existing_deg_table()
        df_sig = df_gene[(df_gene["FDR"] < FDR_THRESH) & (df_gene["log2FC"].abs() > FC_THRESH)]
        df_sig.to_csv(os.path.join(OUTDIR, "deg_significant.csv"), index=False)
        print(f"  Significant (FDR<0.05, |log2FC|>1) : {len(df_sig)}")
        print("\n[3/3] Volcano plot …")
        volcano_plot(df_gene)
        print("\n  Done using existing DEG table. All outputs saved to:", OUTDIR)
        import sys; sys.exit(0)

    print("\n[1/3] Differential expression analysis …")
    df_deg_probe = run_deg(expr, groups)

    # Reload GSE for platform mapping
    cache_dir = os.path.join(OUTDIR, "geo_cache")
    gse = GEOparse.get_GEO(GEO_ID, destdir=cache_dir, silent=True)
    probe_map = map_probes_to_genes(gse)

    print("\n[2/3] Aggregating probes → genes …")
    df_gene = aggregate_to_gene(df_deg_probe, probe_map)
    df_gene.to_csv(os.path.join(OUTDIR, "deg_full_results.csv"))

    df_sig = df_gene[(df_gene["FDR"] < FDR_THRESH) & (df_gene["log2FC"].abs() > FC_THRESH)]
    df_sig.to_csv(os.path.join(OUTDIR, "deg_significant.csv"))
    print(f"  Total genes tested : {len(df_gene)}")
    print(f"  Significant (FDR<0.05, |log2FC|>1) : {len(df_sig)}")

    print("\n[3/3] Volcano plot …")
    volcano_plot(df_gene)

    print("\n  Done. All outputs saved to:", OUTDIR)
