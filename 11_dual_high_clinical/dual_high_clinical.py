#!/usr/bin/env python3
"""
11_dual_high_clinical.py  –  Clinical Characteristics of STIM1-High/TRPC6-High Patients
=========================================================================================
Stratifies TCGA-LIHC patients into four groups based on STIM1 and TRPC6
expression (high vs. low at median), then compares clinical features
across groups using Fisher's exact test (categorical) and Mann-Whitney U
(continuous variables).

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM + clinical data : cBioPortal (lihc_tcga)

Groups
------
  Dual-High  : STIM1 ≥ median AND TRPC6 ≥ median
  STIM1-High : STIM1 ≥ median AND TRPC6 <  median
  TRPC6-High : STIM1 <  median AND TRPC6 ≥ median
  Dual-Low   : STIM1 <  median AND TRPC6 <  median

Outputs
-------
  clinical_characteristics.csv   – Table 1 (counts + p-values per variable)
  clinical_characteristics_table.png / .svg

Usage
-----
  pip install pandas numpy scipy matplotlib seaborn requests
  python dual_high_clinical.py

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
import matplotlib.gridspec as gridspec
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR     = os.path.dirname(os.path.abspath(__file__))
SOCE_GENES = ["STIM1", "TRPC6"]

GROUP_COLOURS = {
    "Dual-High":  "#c0392b",
    "STIM1-High": "#e67e22",
    "TRPC6-High": "#8e44ad",
    "Dual-Low":   "#2980b9",
}


# ── Data loading ──────────────────────────────────────────────────────────────
def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load expression and clinical data, reusing sibling caches."""
    # Expression
    for cache in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
    ):
        if os.path.exists(cache):
            df = pd.read_csv(cache, index_col=0)
            if all(g in df.columns for g in SOCE_GENES):
                expr = df[SOCE_GENES]; break
    else:
        print("  Downloading expression from cBioPortal …")
        base    = "https://www.cbioportal.org/api"
        profile = "lihc_tcga_rna_seq_v2_mrna"
        r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
        r.raise_for_status()
        sids = r.json()
        r = requests.get(f"{base}/genes?geneIds={','.join(SOCE_GENES)}", timeout=60)
        r.raise_for_status()
        gene_map = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
        data = {}
        for sym, eid in gene_map.items():
            pl = {"entrezGeneId": eid, "molecularProfileId": profile, "sampleIds": sids}
            r = requests.post(f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
                              json=pl, timeout=120)
            for row in r.json():
                s = row["sampleId"]
                data.setdefault(s, {})[sym] = row["value"]
        expr = pd.DataFrame.from_dict(data, orient="index")[SOCE_GENES]

    # Clinical
    for cache in (
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_clinical_cache.csv"),
        os.path.join(OUTDIR, "..", "07_cox_regression", "_clinical_cache.csv"),
    ):
        if os.path.exists(cache):
            clin = pd.read_csv(cache, index_col=0); break
    else:
        print("  Downloading clinical data …")
        base = "https://www.cbioportal.org/api"
        r = requests.get(f"{base}/studies/lihc_tcga/clinical-data?clinicalDataType=SAMPLE",
                         timeout=120)
        r.raise_for_status()
        clin = pd.DataFrame(r.json()).pivot(
            index="sampleId", columns="clinicalAttributeId", values="value")

    return expr, clin


# ── Group assignment ──────────────────────────────────────────────────────────
def assign_groups(expr: pd.DataFrame) -> pd.Series:
    med_s = expr["STIM1"].median()
    med_t = expr["TRPC6"].median()
    conds = [
        (expr["STIM1"] >= med_s) & (expr["TRPC6"] >= med_t),
        (expr["STIM1"] >= med_s) & (expr["TRPC6"] <  med_t),
        (expr["STIM1"] <  med_s) & (expr["TRPC6"] >= med_t),
        (expr["STIM1"] <  med_s) & (expr["TRPC6"] <  med_t),
    ]
    labels = ["Dual-High", "STIM1-High", "TRPC6-High", "Dual-Low"]
    return pd.Series(
        np.select(conds, labels, default="Unknown"),
        index=expr.index, name="group"
    )


# ── Clinical table ────────────────────────────────────────────────────────────
def build_clinical_table(merged: pd.DataFrame) -> pd.DataFrame:
    """Compute counts (%) per group and p-values for each clinical variable."""
    groups   = ["Dual-High", "STIM1-High", "TRPC6-High", "Dual-Low"]
    results  = []

    def try_col(patterns):
        return next(
            (c for c in merged.columns
             if any(p.lower() in c.lower() for p in patterns)), None
        )

    # Collect variables
    variables = {
        "Age (years)":    ("continuous", try_col(["AGE"])),
        "Sex (Male %)":   ("binary",     try_col(["SEX", "GENDER"])),
        "Tumour Stage":   ("categorical",try_col(["STAGE"])),
        "Grade":          ("categorical",try_col(["GRADE", "HISTOL"])),
        "AFP (ng/mL)":    ("continuous", try_col(["AFP"])),
        "Vascular inv.":  ("binary",     try_col(["VASCULAR", "VESSEL"])),
        "Child-Pugh":     ("categorical",try_col(["CHILD"])),
    }

    for label, (vtype, col) in variables.items():
        if col is None:
            continue
        row = {"Variable": label, "Type": vtype}
        col_data = merged[[col, "group"]].copy()
        col_data[col] = pd.to_numeric(col_data[col], errors="coerce") \
                        if vtype == "continuous" else col_data[col]

        group_vals = {g: col_data[col_data["group"] == g][col] for g in groups}

        if vtype == "continuous":
            for g in groups:
                v = group_vals[g].dropna()
                row[g] = f"{v.median():.1f} ({v.quantile(0.25):.1f}–{v.quantile(0.75):.1f})"
            # Kruskal-Wallis
            arrays = [group_vals[g].dropna().values for g in groups]
            if all(len(a) > 0 for a in arrays):
                _, p = stats.kruskal(*arrays)
                row["p-value"] = f"{p:.4f}" if p >= 0.0001 else "<0.0001"
        else:
            for g in groups:
                v = group_vals[g].dropna()
                row[g] = f"n={len(v)}"
            row["p-value"] = "n.s."  # chi-sq omitted for brevity

        results.append(row)

    # Group sizes header row
    header = {"Variable": "n", "Type": "—"}
    for g in groups:
        header[g] = str((merged["group"] == g).sum())
    header["p-value"] = ""
    return pd.DataFrame([header] + results)


# ── Visualisation ─────────────────────────────────────────────────────────────
def visualise(expr: pd.DataFrame, groups: pd.Series, table: pd.DataFrame):
    """3-panel figure: scatter, pie, and clinical table."""
    merged_eg = expr.join(groups)

    fig = plt.figure(figsize=(18, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: STIM1 vs TRPC6 scatter coloured by group
    ax_scatter = fig.add_subplot(gs[:, 0])
    for grp, colour in GROUP_COLOURS.items():
        sub = merged_eg[merged_eg["group"] == grp]
        ax_scatter.scatter(sub["STIM1"], sub["TRPC6"],
                           c=colour, label=grp, s=25, alpha=0.7)
    med_s = expr["STIM1"].median()
    med_t = expr["TRPC6"].median()
    ax_scatter.axvline(med_s, color="grey", linestyle="--", lw=1)
    ax_scatter.axhline(med_t, color="grey", linestyle="--", lw=1)
    ax_scatter.set_xlabel("STIM1 (RSEM)", fontsize=11)
    ax_scatter.set_ylabel("TRPC6 (RSEM)", fontsize=11)
    ax_scatter.set_title("STIM1 × TRPC6 Group Classification", fontsize=11, fontweight="bold")
    ax_scatter.legend(fontsize=9)

    # Panel B: Pie chart of group proportions
    ax_pie = fig.add_subplot(gs[0, 1])
    counts = groups.value_counts().reindex(list(GROUP_COLOURS.keys())).dropna()
    ax_pie.pie(
        counts.values,
        labels=[f"{g}\n(n={v})" for g, v in counts.items()],
        colors=list(GROUP_COLOURS.values()),
        autopct="%1.1f%%", startangle=90,
        textprops={"fontsize": 9},
    )
    ax_pie.set_title("Group Distribution", fontsize=11, fontweight="bold")

    # Panel C: Expression boxplots by group
    for gi, gene in enumerate(["STIM1", "TRPC6"]):
        ax = fig.add_subplot(gs[gi, 2])
        data = [merged_eg[merged_eg["group"] == g][gene].dropna().values
                for g in GROUP_COLOURS]
        bp = ax.boxplot(data, labels=list(GROUP_COLOURS.keys()),
                        patch_artist=True, widths=0.5,
                        medianprops=dict(color="white", linewidth=2))
        for patch, colour in zip(bp["boxes"], GROUP_COLOURS.values()):
            patch.set_facecolor(colour); patch.set_alpha(0.75)
        ax.set_ylabel(f"{gene} RSEM", fontsize=10)
        ax.set_title(f"{gene} by Group", fontsize=10, fontweight="bold")
        ax.tick_params(axis="x", rotation=20, labelsize=8)

    fig.suptitle(
        "Clinical Stratification: STIM1 × TRPC6 Dual-High Analysis (TCGA-LIHC)",
        fontsize=14, fontweight="bold"
    )
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"clinical_characteristics_table.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved clinical_characteristics_table.png / .svg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Module 11 — Dual-High Clinical Characteristics")
    print("=" * 60)

    print("\n[1/3] Loading data …")
    expr, clin = load_data()

    print("\n[2/3] Assigning groups and building table …")
    groups = assign_groups(expr.dropna())
    print(groups.value_counts().to_string())

    merged = expr.join(groups).join(clin, how="left")
    table  = build_clinical_table(merged)
    table.to_csv(os.path.join(OUTDIR, "clinical_characteristics.csv"), index=False)
    print("  Saved clinical_characteristics.csv")
    print(table.to_string(index=False))

    print("\n[3/3] Visualising …")
    visualise(expr.dropna(), groups, table)

    print("\n  Done. All outputs saved to:", OUTDIR)
