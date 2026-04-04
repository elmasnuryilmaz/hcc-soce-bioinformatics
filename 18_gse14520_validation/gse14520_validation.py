"""
GSE14520 Independent Validation of SOCE-EMT Findings
=====================================================
Dataset: GSE14520 (Roessler et al., HCC microarray, n=247 tumor + 40 normal)
Platform: Affymetrix HG-U133A 2.0 (GPL571)

Analyses:
1. SOCE gene expression: tumor vs paired-normal (Mann-Whitney U)
2. Kaplan-Meier survival (STIM1, TRPC6 high vs low, OS)
3. Spearman co-expression heatmap (STIM1, TRPC6 vs EMT panel)
4. Compare with TCGA-LIHC results
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
import GEOparse
import warnings
warnings.filterwarnings('ignore')

OUT = "mnt/Beste_Bioinformatics/18_gse14520_validation"
os.makedirs(OUT, exist_ok=True)

# ─── COLOUR PALETTE ──────────────────────────────────────────────────────────
TUMOR_COL  = "#C0392B"
NORMAL_COL = "#2980B9"
HIGH_COL   = "#E74C3C"
LOW_COL    = "#3498DB"

# ─── GENE PANELS ─────────────────────────────────────────────────────────────
SOCE_GENES = ["STIM1", "TRPC6", "TRPC1", "ORAI1"]
EMT_GENES  = ["CDH1", "CDH2", "VIM", "FN1", "SNAI1", "ZEB1", "ZEB2",
               "TWIST1", "MMP2", "NFKB1", "RELA", "IL6R", "ABCC2",
               "STAT3", "AKT1", "HIF1A"]
ALL_GENES  = list(set(SOCE_GENES + EMT_GENES))

# ─────────────────────────────────────────────────────────────────────────────
print("=" * 65)
print("  GSE14520 Independent Validation — SOCE & EMT in HCC")
print("=" * 65)

# ─── 1. DOWNLOAD GSE14520 ────────────────────────────────────────────────────
print("\n[1] Downloading GSE14520 from GEO …")
cache = os.path.join(OUT, "geo_cache")
os.makedirs(cache, exist_ok=True)

gse = GEOparse.get_GEO("GSE14520", destdir=cache, silent=True)
print(f"    GSMs: {len(gse.gsms)}")

# ─── 2. EXTRACT EXPRESSION MATRIX ────────────────────────────────────────────
print("\n[2] Extracting expression data …")

# Pivot: rows = probes, columns = samples
expr_list = []
for gsm_name, gsm in gse.gsms.items():
    df_tmp = gsm.table[["ID_REF", "VALUE"]].copy()
    df_tmp.columns = ["probe", gsm_name]
    df_tmp = df_tmp.set_index("probe")
    expr_list.append(df_tmp)

expr_raw = pd.concat(expr_list, axis=1).astype(float)
print(f"    Shape: {expr_raw.shape[0]} probes × {expr_raw.shape[1]} samples")

# ─── 3. MAP PROBES → GENE SYMBOLS ────────────────────────────────────────────
print("\n[3] Mapping probes to gene symbols …")

# Use GPL571 annotation (platform table embedded in GEO object)
platform_key = list(gse.gpls.keys())[0]
gpl = gse.gpls[platform_key]
annot = gpl.table[["ID", "Gene Symbol"]].dropna()
annot.columns = ["probe", "symbol"]
annot["symbol"] = annot["symbol"].str.strip()
# Keep only probes matching our genes of interest
annot_filt = annot[annot["symbol"].isin(ALL_GENES)]

# For multi-mapped probes, keep highest-variance probe per gene
expr_annot = expr_raw.loc[expr_raw.index.isin(annot_filt["probe"])]
probe2gene = annot_filt.set_index("probe")["symbol"].to_dict()
expr_annot.index = expr_annot.index.map(probe2gene)
expr_annot = expr_annot[~expr_annot.index.isna()]

# Select highest-variance probe per gene
var = expr_annot.var(axis=1)
expr_annot = expr_annot.loc[var.groupby(expr_annot.index).idxmax()]
print(f"    Genes mapped: {expr_annot.shape[0]} / {len(ALL_GENES)} requested")
print(f"    Available SOCE: {[g for g in SOCE_GENES if g in expr_annot.index]}")
print(f"    Available EMT:  {[g for g in EMT_GENES if g in expr_annot.index]}")

expr_annot = expr_annot.T   # samples × genes

# ─── 4. PARSE CLINICAL / SAMPLE METADATA ─────────────────────────────────────
print("\n[4] Parsing clinical metadata …")

meta_rows = []
for gsm_name, gsm in gse.gsms.items():
    ch = gsm.metadata
    title      = ch.get("title", [""])[0]
    tissue     = ch.get("source_name_ch1", [""])[0].lower()
    chars      = {k: v[0] if v else "" for k, v in ch.items()}
    # characteristics_ch1 is a list of "key: value" strings
    char_dict = {}
    for item in ch.get("characteristics_ch1", []):
        if ":" in item:
            k, v = item.split(":", 1)
            char_dict[k.strip().lower()] = v.strip()
    meta_rows.append({
        "gsm": gsm_name,
        "title": title,
        "tissue": tissue,
        **char_dict
    })

meta = pd.DataFrame(meta_rows).set_index("gsm")
print(f"    Total samples: {len(meta)}")
print(f"    Columns: {list(meta.columns[:10])}")
print(meta["tissue"].value_counts().head(10))

# Determine tumor vs normal
meta["sample_type"] = meta["tissue"].apply(
    lambda x: "Normal" if "non-tumor" in x or "normal" in x else "Tumor"
)
print(f"\n    Tumor:  {(meta['sample_type']=='Tumor').sum()}")
print(f"    Normal: {(meta['sample_type']=='Normal').sum()}")

# Survival columns — try multiple field names
for col in meta.columns:
    if "surv" in col.lower() or "time" in col.lower() or "os" in col.lower() or "month" in col.lower():
        print(f"    Survival-related col: {col} | sample: {meta[col].iloc[0]}")

# ─── 5. TUMOR vs NORMAL EXPRESSION ───────────────────────────────────────────
print("\n[5] Tumor vs Normal expression (Mann-Whitney U) …")

common = expr_annot.index.intersection(meta.index)
expr = expr_annot.loc[common]
meta = meta.loc[common]

tumor_idx  = meta[meta["sample_type"] == "Tumor"].index
normal_idx = meta[meta["sample_type"] == "Normal"].index

results_tvn = []
for gene in SOCE_GENES:
    if gene not in expr.columns:
        continue
    t_vals = expr.loc[tumor_idx, gene].dropna()
    n_vals = expr.loc[normal_idx, gene].dropna()
    stat, pval = stats.mannwhitneyu(t_vals, n_vals, alternative="two-sided")
    log2fc = np.log2(t_vals.median() + 1) - np.log2(n_vals.median() + 1)
    results_tvn.append({
        "gene": gene,
        "n_tumor": len(t_vals),
        "n_normal": len(n_vals),
        "median_tumor": t_vals.median(),
        "median_normal": n_vals.median(),
        "log2FC": log2fc,
        "p_value": pval,
        "significant": pval < 0.001
    })
    stars = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
    print(f"    {gene:6s}: log2FC={log2fc:+.2f}, p={pval:.2e} {stars}")

df_tvn = pd.DataFrame(results_tvn)
df_tvn.to_csv(f"{OUT}/tumor_vs_normal_results.csv", index=False)

# ─── 6. PLOT: TUMOR vs NORMAL BOXPLOT ────────────────────────────────────────
print("\n[6] Plotting tumor vs normal boxplot …")
avail_soce = [g for g in SOCE_GENES if g in expr.columns]
fig, axes = plt.subplots(1, len(avail_soce), figsize=(3.5 * len(avail_soce), 4.5))
if len(avail_soce) == 1:
    axes = [axes]

for ax, gene in zip(axes, avail_soce):
    t_vals = expr.loc[tumor_idx, gene].dropna()
    n_vals = expr.loc[normal_idx, gene].dropna()
    data = pd.DataFrame({
        "Expression": pd.concat([n_vals, t_vals]),
        "Group": ["Normal"] * len(n_vals) + ["Tumor"] * len(t_vals)
    })
    sns.boxplot(data=data, x="Group", y="Expression",
                palette={"Normal": NORMAL_COL, "Tumor": TUMOR_COL},
                width=0.5, fliersize=2, ax=ax)
    stat, pval = stats.mannwhitneyu(t_vals, n_vals, alternative="two-sided")
    stars = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
    y_max = data["Expression"].max()
    ax.plot([0, 1], [y_max * 1.05, y_max * 1.05], color="black", lw=1.2)
    ax.text(0.5, y_max * 1.08, stars, ha="center", fontsize=13, fontweight="bold")
    ax.set_title(gene, fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Expression (log2)" if gene == avail_soce[0] else "")
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("GSE14520: SOCE Expression in HCC vs Normal Liver",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()
fig.savefig(f"{OUT}/tumor_vs_normal_boxplot.png", dpi=180, bbox_inches="tight")
fig.savefig(f"{OUT}/tumor_vs_normal_boxplot.svg", bbox_inches="tight")
plt.close()
print("    Saved tumor_vs_normal_boxplot.png")

# ─── 7. KAPLAN-MEIER SURVIVAL ─────────────────────────────────────────────────
print("\n[7] Kaplan-Meier survival analysis …")

tumor_meta = meta[meta["sample_type"] == "Tumor"].copy()

# Try to identify OS time and event columns
print("    Tumor meta columns:", list(tumor_meta.columns))

# GSE14520 survival columns: "overall survival time (month)" and "survival status"
os_time_col   = None
os_event_col  = None
for col in tumor_meta.columns:
    cl = col.lower()
    if "survival time" in cl or ("time" in cl and "month" in cl) or "os_month" in cl:
        os_time_col = col
    if "survival status" in cl or "status" in cl or "event" in cl or "censor" in cl:
        os_event_col = col

print(f"    OS time col  : {os_time_col}")
print(f"    OS event col : {os_event_col}")

if os_time_col and os_event_col:
    tumor_meta["os_time"] = pd.to_numeric(tumor_meta[os_time_col], errors="coerce")
    tumor_meta["os_event"] = tumor_meta[os_event_col].apply(
        lambda x: 1 if str(x).lower() in ["1", "dead", "death", "deceased", "died", "yes"] else 0
    )
    surv_df = tumor_meta[["os_time", "os_event"]].dropna()
    surv_df = surv_df[surv_df["os_time"] > 0]
    surv_df = surv_df.join(expr[avail_soce])
    print(f"    Survival samples: {len(surv_df)}")
    print(f"    Events: {surv_df['os_event'].sum()} / {len(surv_df)}")

    km_results = []
    fig, axes = plt.subplots(1, min(2, len(avail_soce)),
                             figsize=(5 * min(2, len(avail_soce)), 4.5))
    if len(avail_soce) < 2:
        axes = [axes]

    for ax, gene in zip(axes, [g for g in ["STIM1", "TRPC6"] if g in avail_soce]):
        med = surv_df[gene].median()
        high = surv_df[gene] >= med
        low  = surv_df[gene] < med
        kmf = KaplanMeierFitter()
        kmf.fit(surv_df.loc[high, "os_time"], surv_df.loc[high, "os_event"],
                label=f"{gene}-high (n={high.sum()})")
        kmf.plot_survival_function(ax=ax, color=HIGH_COL, ci_show=True)
        kmf.fit(surv_df.loc[low, "os_time"], surv_df.loc[low, "os_event"],
                label=f"{gene}-low (n={low.sum()})")
        kmf.plot_survival_function(ax=ax, color=LOW_COL, ci_show=True)
        lr = logrank_test(surv_df.loc[high, "os_time"], surv_df.loc[low, "os_time"],
                          event_observed_A=surv_df.loc[high, "os_event"],
                          event_observed_B=surv_df.loc[low, "os_event"])
        # HR from univariate Cox
        cph = CoxPHFitter()
        cox_df = surv_df[["os_time", "os_event", gene]].dropna()
        cph.fit(cox_df, duration_col="os_time", event_col="os_event")
        hr  = np.exp(cph.params_[gene])
        hr_lo = np.exp(cph.confidence_intervals_.iloc[0, 0])
        hr_hi = np.exp(cph.confidence_intervals_.iloc[0, 1])
        p = lr.p_value
        stars = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        ax.set_title(f"{gene} — GSE14520\nHR={hr:.2f} ({hr_lo:.2f}–{hr_hi:.2f}), p={p:.3f} {stars}",
                     fontsize=10)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Overall Survival Probability")
        ax.spines[["top", "right"]].set_visible(False)
        km_results.append({"gene": gene, "p_logrank": p, "HR": hr, "HR_lo": hr_lo, "HR_hi": hr_hi})
        print(f"    {gene}: HR={hr:.2f} ({hr_lo:.2f}–{hr_hi:.2f}), p={p:.4f} {stars}")

    plt.tight_layout()
    fig.savefig(f"{OUT}/km_survival_gse14520.png", dpi=180, bbox_inches="tight")
    fig.savefig(f"{OUT}/km_survival_gse14520.svg", bbox_inches="tight")
    plt.close()
    pd.DataFrame(km_results).to_csv(f"{OUT}/km_survival_results.csv", index=False)
    print("    Saved km_survival_gse14520.png")
else:
    print("    ⚠ Could not identify OS columns — skipping KM")
    print("    Available cols:", list(tumor_meta.columns))
    km_results = []

# ─── 8. SPEARMAN CO-EXPRESSION ────────────────────────────────────────────────
print("\n[8] Spearman co-expression (STIM1, TRPC6 vs EMT panel) …")

tumor_expr = expr.loc[tumor_idx].copy()
avail_emt  = [g for g in EMT_GENES if g in tumor_expr.columns]
rows = []
for soce in avail_soce:
    for emt in avail_emt:
        rho, pval = stats.spearmanr(tumor_expr[soce], tumor_expr[emt])
        rows.append({"SOCE": soce, "partner": emt, "rho": rho, "p": pval,
                     "sig": pval < 0.01 and abs(rho) >= 0.2})
corr_df = pd.DataFrame(rows)
corr_df.to_csv(f"{OUT}/spearman_correlations.csv", index=False)

# Print key correlations
for soce in avail_soce:
    top = corr_df[corr_df["SOCE"] == soce].sort_values("rho", ascending=False).head(5)
    print(f"\n    {soce} top partners:")
    for _, r in top.iterrows():
        print(f"      {r['partner']:8s}: ρ={r['rho']:.3f}, p={r['p']:.2e}")

# Pivot for heatmap
heatmap_soce = [g for g in ["STIM1", "TRPC6"] if g in avail_soce]
piv = corr_df[corr_df["SOCE"].isin(heatmap_soce)].pivot(index="SOCE", columns="partner", values="rho")
piv_p = corr_df[corr_df["SOCE"].isin(heatmap_soce)].pivot(index="SOCE", columns="partner", values="p")

# Sort columns by STIM1 rho
if "STIM1" in piv.index:
    order = piv.loc["STIM1"].sort_values(ascending=False).index
else:
    order = piv.columns

piv = piv[order]
piv_p = piv_p[order]

fig, ax = plt.subplots(figsize=(max(10, len(order) * 0.7), 3.5))
sns.heatmap(piv, cmap="RdBu_r", center=0, vmin=-0.6, vmax=0.6,
            annot=True, fmt=".2f", linewidths=0.5,
            annot_kws={"size": 8}, ax=ax)

# Add significance asterisks
for i, soce in enumerate(piv.index):
    for j, gene in enumerate(piv.columns):
        p = piv_p.loc[soce, gene]
        if p < 0.01:
            ax.text(j + 0.75, i + 0.25, "*", ha="center", fontsize=9,
                    color="black", fontweight="bold")

ax.set_title("GSE14520: Spearman Co-expression (SOCE vs EMT Panel)\n"
             "* p < 0.01", fontsize=11, fontweight="bold")
ax.set_xlabel("")
ax.set_ylabel("")
plt.tight_layout()
fig.savefig(f"{OUT}/spearman_heatmap.png", dpi=180, bbox_inches="tight")
fig.savefig(f"{OUT}/spearman_heatmap.svg", bbox_inches="tight")
plt.close()
print("\n    Saved spearman_heatmap.png")

# ─── 9. COMPARISON TABLE: TCGA vs GSE14520 ───────────────────────────────────
print("\n[9] Building TCGA vs GSE14520 comparison …")

tcga_corr = {
    ("TRPC6", "VIM"):   0.48, ("TRPC6", "ZEB2"):  0.45,
    ("TRPC6", "ZEB1"):  0.39, ("TRPC6", "TWIST1"): 0.39,
    ("TRPC6", "SNAI1"): 0.31, ("STIM1", "IL6R"):  0.41,
    ("STIM1", "ABCC2"): 0.37, ("STIM1", "RELA"):  0.25,
    ("STIM1", "NFKB1"): 0.16,
}

comp_rows = []
for (soce, partner), tcga_rho in tcga_corr.items():
    gse_row = corr_df[(corr_df["SOCE"] == soce) & (corr_df["partner"] == partner)]
    gse_rho = gse_row["rho"].values[0] if len(gse_row) else np.nan
    gse_p   = gse_row["p"].values[0] if len(gse_row) else np.nan
    direction_match = (tcga_rho > 0) == (gse_rho > 0) if not np.isnan(gse_rho) else None
    comp_rows.append({
        "SOCE": soce, "Partner": partner,
        "TCGA_rho": tcga_rho,
        "GSE14520_rho": round(gse_rho, 3) if not np.isnan(gse_rho) else "N/A",
        "GSE14520_p": f"{gse_p:.3e}" if not np.isnan(gse_p) else "N/A",
        "Direction_match": direction_match,
        "Validated": direction_match and gse_p < 0.05 if direction_match is not None else False
    })

comp_df = pd.DataFrame(comp_rows)
comp_df.to_csv(f"{OUT}/tcga_vs_gse14520_comparison.csv", index=False)

print("\n    TCGA vs GSE14520 key correlations:")
print(comp_df[["SOCE", "Partner", "TCGA_rho", "GSE14520_rho", "Direction_match", "Validated"]].to_string(index=False))

# ─── 10. SAVE SUMMARY ────────────────────────────────────────────────────────
validated = comp_df["Validated"].sum()
total = len(comp_df)
print(f"\n{'='*65}")
print(f"  VALIDATION SUMMARY")
print(f"  Correlations validated: {validated}/{total} ({100*validated/total:.0f}%)")
print(f"{'='*65}")

with open(f"{OUT}/validation_summary.txt", "w") as f:
    f.write("GSE14520 INDEPENDENT VALIDATION SUMMARY\n")
    f.write("=" * 50 + "\n\n")
    f.write("Dataset: GSE14520 (Roessler et al.)\n")
    f.write(f"Samples: {(meta['sample_type']=='Tumor').sum()} tumor, "
            f"{(meta['sample_type']=='Normal').sum()} normal\n\n")
    f.write("TUMOR vs NORMAL:\n")
    for _, r in df_tvn.iterrows():
        stars = "***" if r['p_value'] < 0.001 else ("**" if r['p_value'] < 0.01 else "*")
        f.write(f"  {r['gene']}: log2FC={r['log2FC']:+.2f}, p={r['p_value']:.2e} {stars}\n")
    f.write("\nCO-EXPRESSION VALIDATION:\n")
    f.write(comp_df[["SOCE","Partner","TCGA_rho","GSE14520_rho","Validated"]].to_string(index=False))
    f.write(f"\n\nOverall: {validated}/{total} key correlations validated\n")

print("\nAll outputs saved to:", OUT)
print("Done.")
