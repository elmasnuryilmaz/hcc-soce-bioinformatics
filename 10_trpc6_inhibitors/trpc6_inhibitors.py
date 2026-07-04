#!/usr/bin/env python3
"""
10_trpc6_inhibitors.py  –  TRPC6 Inhibitor Landscape (ChEMBL)
=============================================================
  Queries the ChEMBL REST API for compounds that inhibit TRPC6 (CHEMBL2417347),
filters for drug-like molecules, and visualises the potency landscape.

Data source
-----------
  ChEMBL v33 REST API : https://www.ebi.ac.uk/chembl/api/data/

Outputs
-------
  chembl_trpc6_inhibitors.csv          – All TRPC6 bioactivities
  chembl_trpc6_specific_inhibitors.csv – Filtered drug-like subset
  trpc6_inhibitor_summary.csv          – Summary by assay type
  trpc6_inhibitor_landscape.png / .svg

Usage
-----
  pip install pandas numpy matplotlib seaborn requests
  python trpc6_inhibitors.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings, math
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

OUTDIR          = os.path.dirname(os.path.abspath(__file__))
TRPC6_CHEMBL_ID = "CHEMBL2417347"      # human TRPC6 ChEMBL target ID
CHEMBL_API      = "https://www.ebi.ac.uk/chembl/api/data"
MAX_PAGES       = 20
PAGE_SIZE       = 1000


def fetch_bioactivities(target_id: str) -> pd.DataFrame:
    """Fetches all bioactivities for the given ChEMBL target."""
    cache = os.path.join(OUTDIR, "_chembl_raw_cache.csv")
    if os.path.exists(cache):
        print(f"  [cache] {cache}")
        return pd.read_csv(cache)

    records = []
    url = f"{CHEMBL_API}/activity.json"
    params = {
        "target_chembl_id": target_id,
        "limit":            PAGE_SIZE,
        "offset":           0,
        "format":           "json",
    }

    print(f"  Fetching bioactivities for {target_id} from ChEMBL …")
    for page in range(MAX_PAGES):
        params["offset"] = page * PAGE_SIZE
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        activities = data.get("activities", [])
        records.extend(activities)
        print(f"    Page {page+1}: {len(activities)} activities (total: {len(records)})")
        if len(activities) < PAGE_SIZE:
            break  # last page

    df = pd.DataFrame(records)
    df.to_csv(cache, index=False)
    print(f"  Total bioactivities fetched: {len(df)}")
    return df


def filter_inhibitors(df: pd.DataFrame) -> pd.DataFrame:
    """Keep IC50 / Ki measurements with defined pChEMBL ≥ 5."""
    if df.empty or "standard_type" not in df.columns:
        return pd.DataFrame()
    keep_types = {"IC50", "Ki", "Kd", "EC50", "Potency"}
    df = df[df["standard_type"].isin(keep_types)].copy()
    df["pchembl_value"] = pd.to_numeric(df.get("pchembl_value", np.nan), errors="coerce")
    df = df[df["pchembl_value"] >= 5]   # pIC50 ≥ 5 → IC50 ≤ 10 µM
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["pchembl_value"])
    df = df.sort_values("pchembl_value", ascending=False)
    return df


def summarise(df_filt: pd.DataFrame) -> pd.DataFrame:
    """Summary by standard_type."""
    if df_filt.empty or "standard_type" not in df_filt.columns:
        return pd.DataFrame(columns=["n_compounds", "median_pchembl", "best_pchembl"])
    summary = df_filt.groupby("standard_type").agg(
        n_compounds=("molecule_chembl_id", "nunique"),
        median_pchembl=("pchembl_value", "median"),
        best_pchembl=("pchembl_value", "max"),
    ).round(2)
    return summary


def landscape_plot(df_filt: pd.DataFrame):
    """Scatter: pChEMBL value per compound, coloured by assay type."""
    if df_filt.empty:
        print("  No filtered data to plot.")
        return

    # Rank compounds by potency
    df_plot = df_filt.copy().reset_index(drop=True)
    df_plot = df_plot.sort_values("pchembl_value", ascending=False).head(100)
    df_plot["rank"] = range(1, len(df_plot) + 1)

    colours = sns.color_palette("Set2", n_colors=df_plot["standard_type"].nunique())
    type_colour = dict(zip(df_plot["standard_type"].unique(), colours))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Left: scatter of top compounds
    ax = axes[0]
    for st, grp in df_plot.groupby("standard_type"):
        ax.scatter(grp["rank"], grp["pchembl_value"], label=st,
                   color=type_colour[st], s=50, alpha=0.8)
    ax.axhline(6, color="grey", linestyle="--", linewidth=0.8, label="pIC50 = 6 (100 nM)")
    ax.axhline(7, color="orange", linestyle="--", linewidth=0.8, label="pIC50 = 7 (10 nM)")
    ax.set_xlabel("Compound rank (by potency)", fontsize=11)
    ax.set_ylabel("pChEMBL value (−log₁₀ IC50 / Ki / EC50 in M)", fontsize=11)
    ax.set_title("TRPC6 Inhibitor Potency Landscape\n(ChEMBL, top 100 compounds)", fontsize=11)
    ax.legend(fontsize=9)

    # Right: histogram of pChEMBL distribution
    ax2 = axes[1]
    for st, grp in df_filt.groupby("standard_type"):
        ax2.hist(grp["pchembl_value"], bins=20, alpha=0.6, label=st,
                 color=type_colour.get(st, "grey"), edgecolor="white")
    ax2.set_xlabel("pChEMBL value", fontsize=11)
    ax2.set_ylabel("Count", fontsize=11)
    ax2.set_title(
        f"pChEMBL Distribution\n(n = {df_filt['molecule_chembl_id'].nunique()} unique compounds)",
        fontsize=11
    )
    ax2.legend(fontsize=9)

    fig.suptitle(
        "TRPC6 (CHEMBL2364817) — Inhibitor Landscape from ChEMBL v33",
        fontsize=13, fontweight="bold"
    )
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"trpc6_inhibitor_landscape.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved trpc6_inhibitor_landscape.png / .svg")


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 10 — TRPC6 Inhibitor Landscape (ChEMBL)")
    print("=" * 60)

    print(f"\n[1/3] Fetching TRPC6 bioactivities (target: {TRPC6_CHEMBL_ID}) …")
    df_raw = fetch_bioactivities(TRPC6_CHEMBL_ID)
    df_raw.to_csv(os.path.join(OUTDIR, "chembl_trpc6_inhibitors.csv"), index=False)
    print(f"  Saved chembl_trpc6_inhibitors.csv  ({len(df_raw)} records)")

    print("\n[2/3] Filtering drug-like inhibitors …")
    df_filt = filter_inhibitors(df_raw)
    df_filt.to_csv(os.path.join(OUTDIR, "chembl_trpc6_specific_inhibitors.csv"), index=False)
    n_compounds = df_filt["molecule_chembl_id"].nunique() if "molecule_chembl_id" in df_filt.columns else 0
    print(f"  Filtered to {len(df_filt)} measurements "
          f"({n_compounds} unique compounds, pIC50 ≥ 5)")

    summary = summarise(df_filt)
    summary.to_csv(os.path.join(OUTDIR, "trpc6_inhibitor_summary.csv"))
    print("  Summary:\n", summary.to_string())

    print("\n[3/3] Landscape plot …")
    landscape_plot(df_filt)

    print("\n  Done. All outputs saved to:", OUTDIR)
