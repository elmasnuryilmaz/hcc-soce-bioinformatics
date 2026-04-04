#!/usr/bin/env python3
"""
13_stage_stratified_km.py  –  Stage-Stratified Kaplan-Meier Survival
=====================================================================
Performs Kaplan-Meier survival analysis for STIM1 and TRPC6 within
each AJCC tumour stage (I, II, III, IV) to assess whether the prognostic
associations are stage-dependent.

Data source
-----------
  TCGA-LIHC RNA-seq v2 RSEM + clinical data : cBioPortal (lihc_tcga)

Outputs
-------
  stage_km_results.csv               – Log-rank p-values per gene × stage
  stage_stratified_km.png / .svg     – 2×4 grid of KM curves

Usage
-----
  pip install pandas numpy scipy matplotlib lifelines requests
  python stage_stratified_km.py

Author : Elmasnur Yilmaz (elmasnrylmz@gmail.com)
"""

import os, warnings
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

OUTDIR = os.path.dirname(os.path.abspath(__file__))
GENES  = ["STIM1", "TRPC6"]
STAGES = ["I", "II", "III", "IV"]

COLOURS = {"High": "#c0392b", "Low": "#2980b9"}


def load_data():
    """Load expression and clinical data from sibling caches or cBioPortal."""
    # Expression
    for cache in (
        os.path.join(OUTDIR, "..", "05_coexpression", "_tcga_expr_cache.csv"),
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_expr_cache.csv"),
    ):
        if os.path.exists(cache):
            df = pd.read_csv(cache, index_col=0)
            if all(g in df.columns for g in GENES):
                expr = df[GENES]
                print(f"  [expr cache] {cache}")
                break
    else:
        print("  Downloading expression …")
        base    = "https://www.cbioportal.org/api"
        profile = "lihc_tcga_rna_seq_v2_mrna"
        r = requests.get(f"{base}/sample-lists/lihc_tcga_all/sample-ids", timeout=60)
        r.raise_for_status()
        sids = r.json()
        r = requests.get(f"{base}/genes?geneIds={','.join(GENES)}", timeout=60)
        r.raise_for_status()
        gmap = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in r.json()}
        data = {}
        for sym, eid in gmap.items():
            pl = {"entrezGeneId": eid, "molecularProfileId": profile, "sampleIds": sids}
            r = requests.post(
                f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
                json=pl, timeout=120)
            for row in r.json():
                s = row["sampleId"]
                data.setdefault(s, {})[sym] = row["value"]
        expr = pd.DataFrame.from_dict(data, orient="index")[GENES]

    # Clinical
    for cache in (
        os.path.join(OUTDIR, "..", "01_tcga_survival", "_tcga_lihc_clinical_cache.csv"),
        os.path.join(OUTDIR, "..", "07_cox_regression", "_clinical_cache.csv"),
    ):
        if os.path.exists(cache):
            clin = pd.read_csv(cache, index_col=0)
            print(f"  [clin cache] {cache}")
            break
    else:
        print("  Downloading clinical data …")
        base = "https://www.cbioportal.org/api"
        r = requests.get(
            f"{base}/studies/lihc_tcga/clinical-data?clinicalDataType=SAMPLE",
            timeout=120)
        r.raise_for_status()
        clin = pd.DataFrame(r.json()).pivot(
            index="sampleId", columns="clinicalAttributeId", values="value")

    return expr, clin


def prepare_merged(expr, clin):
    """Merge and add stage + OS columns."""
    merged = expr.join(clin, how="inner")

    os_t = next((c for c in merged.columns if "OS_MONTHS" in c.upper()), None)
    os_e = next((c for c in merged.columns if "OS_STATUS" in c.upper()), None)
    if os_t is None:
        raise RuntimeError("OS columns not found.")

    merged["os_time"]  = pd.to_numeric(merged[os_t], errors="coerce")
    merged["os_event"] = (
        merged[os_e].astype(str).str.extract(r"(\d+)")[0]
        .astype(float).fillna(0).astype(int)
    )

    stage_col = next((c for c in merged.columns if "STAGE" in c.upper()), None)
    if stage_col:
        merged["stage"] = (
            merged[stage_col].astype(str).str.upper()
            .str.extract(r"(I{1,3}V?|IV)$")[0]
        )
    else:
        merged["stage"] = np.nan

    return merged.dropna(subset=["os_time", "os_event"])


def run_stratified_km(merged):
    """Returns results DataFrame and fig."""
    results = []
    fig, axes = plt.subplots(len(GENES), len(STAGES),
                             figsize=(5 * len(STAGES), 5 * len(GENES)))
    fig.suptitle(
        "Stage-Stratified Kaplan-Meier Survival (TCGA-LIHC)\n"
        "STIM1 and TRPC6 — High vs. Low (split at stage-specific median)",
        fontsize=14, fontweight="bold"
    )

    for gi, gene in enumerate(GENES):
        med_global = merged[gene].median()
        for si, stage in enumerate(STAGES):
            ax = axes[gi][si]
            sub = merged[merged["stage"] == stage].copy()

            if len(sub) < 15:
                ax.set_title(f"Stage {stage}\n(n={len(sub)}, insufficient)", fontsize=9)
                ax.axis("off")
                results.append({"gene": gene, "stage": stage, "n": len(sub),
                                 "logrank_p": np.nan, "n_high": 0, "n_low": 0})
                continue

            med = sub[gene].median()
            grp_h = sub[sub[gene] >= med]
            grp_l = sub[sub[gene] <  med]

            kmf_h = KaplanMeierFitter()
            kmf_l = KaplanMeierFitter()
            kmf_h.fit(grp_h["os_time"], grp_h["os_event"],
                      label=f"High (n={len(grp_h)})")
            kmf_l.fit(grp_l["os_time"], grp_l["os_event"],
                      label=f"Low  (n={len(grp_l)})")

            lr = logrank_test(grp_h["os_time"], grp_l["os_time"],
                              grp_h["os_event"], grp_l["os_event"])
            p = lr.p_value

            kmf_h.plot_survival_function(ax=ax, ci_show=True,
                                         color=COLOURS["High"])
            kmf_l.plot_survival_function(ax=ax, ci_show=True,
                                         color=COLOURS["Low"])

            p_str = f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.2e}"
            ax.text(0.97, 0.97, p_str, transform=ax.transAxes,
                    ha="right", va="top", fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.8))
            ax.set_title(f"{gene} — Stage {stage}", fontsize=10, fontweight="bold")
            ax.set_xlabel("Months", fontsize=9)
            ax.set_ylabel("Survival" if si == 0 else "", fontsize=9)

            results.append({"gene": gene, "stage": stage,
                             "n": len(sub), "n_high": len(grp_h), "n_low": len(grp_l),
                             "logrank_p": round(p, 4)})

    plt.tight_layout()
    return pd.DataFrame(results), fig


if __name__ == "__main__":
    print("=" * 60)
    print("  Module 13 — Stage-Stratified Kaplan-Meier Survival")
    print("=" * 60)

    print("\n[1/3] Loading data …")
    expr, clin = load_data()

    print("\n[2/3] Running stage-stratified KM …")
    merged  = prepare_merged(expr, clin)
    results, fig = run_stratified_km(merged)

    results.to_csv(os.path.join(OUTDIR, "stage_km_results.csv"), index=False)
    print("  Saved stage_km_results.csv")
    print(results.to_string(index=False))

    print("\n[3/3] Saving figures …")
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUTDIR, f"stage_stratified_km.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved stage_stratified_km.png / .svg")

    print("\n  Done. All outputs saved to:", OUTDIR)
