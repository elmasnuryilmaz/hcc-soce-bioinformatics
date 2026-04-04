# Bioinformatic Analysis of SOCE Components in Hepatocellular Carcinoma

> **Companion code and data for:**
> Yurdacan Yasar B, Erac Y. *Store-Operated Calcium Entry Regulates Early Adaptive Responses to Sorafenib in Hepatocellular Carcinoma.* (2026, under review)

---

## Overview

This repository contains all bioinformatic analyses supporting the manuscript. The analyses characterise the clinical relevance of Store-Operated Calcium Entry (SOCE) components — **STIM1**, **TRPC6**, **TRPC1**, and **ORAI1** — in hepatocellular carcinoma (HCC) using publicly available multi-omics datasets.

The analyses are organised into **18 independent modules**, each self-contained with its own input data, analysis script(s), and output files. All results are fully reproducible from public data sources listed below.

### Core Datasets Used

| Dataset | Source | n | Access |
|---------|--------|---|--------|
| TCGA-LIHC RNA-seq v2 RSEM | cBioPortal (`lihc_tcga`) | 370 tumours | https://www.cbioportal.org |
| GTEx Liver (v8) | GTEx Portal | 226 normal | https://gtexportal.org |
| GSE140202 (Huh7 sorafenib resistance) | NCBI GEO | 12,998 genes | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140202 |
| GSE14520 (Independent HCC cohort) | NCBI GEO | 247 tumours + 40 normal | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520 |
| DepMap (22Q4) | Cancer Dependency Map | 21 HCC lines | https://depmap.org |
| ChEMBL | EMBL-EBI | — | https://www.ebi.ac.uk/chembl |
| STRING (v12) | EMBL | — | https://string-db.org |

---

## Summary of Key Findings

### Tumour vs. Normal Expression (TCGA-LIHC vs. GTEx Liver)

All four SOCE components are significantly overexpressed in HCC tumours compared to normal liver tissue (Mann-Whitney U-test, all *p* < 0.001):

| Gene | Median (Normal) | Median (Tumour) | log₂FC | *p*-value |
|------|----------------|-----------------|--------|-----------|
| STIM1 | 4.20 | 11.24 | **+7.04** | 1.16 × 10⁻⁹³ |
| ORAI1 | 3.24 | 8.99 | **+5.75** | 1.16 × 10⁻⁹³ |
| TRPC1 | 0.45 | 4.81 | **+4.36** | 1.81 × 10⁻⁹² |
| TRPC6 | 0.06 | 4.11 | **+4.05** | 9.13 × 10⁻⁹² |

### Survival Analysis (Kaplan-Meier, n = 370)

| Gene | High-group median OS | Low-group median OS | Log-rank *p* | HR (95% CI) |
|------|----------------------|----------------------|--------------|-------------|
| STIM1 | 83.5 months | 41.8 months | **0.014** | 0.68 (0.50–0.94) |
| TRPC6 | 80.7 months | 45.1 months | **0.033** | 0.73 |
| TRPC1 | 53.3 months | 58.8 months | 0.190 | ns |
| ORAI1 | 69.5 months | 53.3 months | 0.843 | ns |

> ⚠️ **Important:** Multivariate Cox regression adjusting for tumour stage, grade, age, sex, and AFP showed that STIM1 and TRPC6 are **not** independent prognostic factors (*p* > 0.05 in both models). Tumour stage is the dominant predictor. See `07_cox_regression/`.

### Spearman Co-expression (TCGA-LIHC, n = 370)

Two distinct co-expression axes were identified (all *p* < 0.01):

**TRPC6 — EMT Transcription Factor Axis:**

| Partner | ρ | Biological Relevance |
|---------|---|----------------------|
| VIM | **0.48** | Core mesenchymal marker |
| ZEB2 | **0.45** | EMT master regulator |
| ZEB1 | 0.39 | EMT master regulator |
| TWIST1 | 0.39 | EMT transcription factor |
| SNAI1 | 0.31 | EMT transcription factor |

**STIM1 — NF-κB / IL-6 / Drug Resistance Axis:**

| Partner | ρ | Biological Relevance |
|---------|---|----------------------|
| IL6R | **0.41** | IL-6 receptor / JAK-STAT |
| ABCC2 | 0.37 | ABC drug transporter (MDR) |
| RELA | 0.25 | NF-κB p65 subunit |

### Gene Set Enrichment Analysis (GSE140202 — Huh7 Sorafenib Resistance)

Genes ranked by limma *t*-statistic (sorafenib-resistant vs. sensitive Huh7 cells, n = 12,998 genes):

| Pathway | NES | FDR *q*-value | Direction |
|---------|-----|--------------|-----------|
| TNFα / NF-κB Signalling | **+1.85** | **< 0.001** | Enriched in resistant |
| TGF-β Signalling | −1.52 | 0.097 | Enriched in sensitive |
| IL-6 / JAK / STAT3 | +1.18 | 0.527 | ns |

### CRISPR Essentiality (DepMap, 21 HCC Cell Lines)

TRPC6 is **non-essential** in all 21 tested HCC cell lines (mean gene effect = −0.028; all lines above the −0.5 essentiality threshold), supporting it as a pharmacologically tractable target with a potentially favourable therapeutic window.

---

## Repository Structure

```
.
├── 01_tcga_survival/           # Tumour vs. normal expression + Kaplan-Meier survival
├── 02_geo_deg/                 # Differential expression analysis (GSE140202)
├── 03_ppi_network/             # Protein-protein interaction network (preliminary)
├── 04_gsea/                    # GSEA on sorafenib-resistance DEGs (preliminary)
├── 05_coexpression/            # Spearman co-expression (preliminary version)
├── 06_stim1_trpc6_correlation/ # Direct STIM1–TRPC6 correlation in TCGA
├── 07_cox_regression/          # Multivariate Cox proportional hazards
├── 08_depmap_drug_sensitivity/ # DepMap drug sensitivity correlations
├── 09_soce_risk_score/         # Composite SOCE risk score + Kaplan-Meier
├── 10_trpc6_inhibitors/        # ChEMBL TRPC6 inhibitor landscape
├── 11_dual_high_clinical/      # TRPC6-high + STIM1-high clinical characterisation
├── 12_crispr_essentiality/     # CRISPR gene essentiality (DepMap)
├── 13_stage_stratified_km/     # Stage-stratified Kaplan-Meier analysis
├── 14_tcga_gsea_stim1/         # GSEA (final, GSE140202 + TCGA pathway correlations)
├── 15_clean_correlation/       # Spearman heatmap (final, publication-ready)
├── 16_emtf_correlation/        # EMT transcription factor bubble/dot plot
├── 17_ppi_improved/            # PPI network (final, with Ca²⁺-TF bridge nodes)
├── 18_gse14520_validation/     # Independent cohort validation (GSE14520)
├── FIGURE_comprehensive_6panel.png    # Manuscript Figure 7 (panels A–E)
├── FIGURE_supplementary_4panel.png    # Supplementary figure
├── BIOINFORMATICS_KEY_FINDINGS.txt    # Numerical results summary
└── README.md
```

---

## Module Descriptions

### `01_tcga_survival/` — Tumour vs. Normal Expression & Survival

**Purpose:** Compare SOCE gene expression between TCGA-LIHC tumours and GTEx normal liver; perform Kaplan-Meier overall survival analysis stratified by SOCE expression level.

**Methods:**
- Tumour vs. normal: Mann-Whitney U-test (two-sided); effect size reported as log₂(median tumour / median normal + 1)
- Survival: Kaplan-Meier estimator; high/low split at median expression; log-rank test for group comparison; hazard ratios from univariate Cox proportional hazards model (`lifelines` v0.30.0)

**Outputs:**

| File | Description |
|------|-------------|
| `tumor_vs_normal_stats.csv` | log₂FC and Mann-Whitney statistics per SOCE gene |
| `tumor_vs_normal_boxplot.png/svg` | Publication-ready boxplot (Fig. 7A) |
| `km_survival_results.csv` | Log-rank *p*-values, median OS, significance flags |
| `km_survival_curves.png/svg` | Kaplan-Meier survival curves (Fig. 7B) |

---

### `02_geo_deg/` — Differential Expression (GSE140202)

**Purpose:** Identify genes differentially expressed between sorafenib-resistant and sorafenib-sensitive Huh7 cells to generate a pre-ranked gene list for GSEA.

**Methods:**
- Dataset: GSE140202 (NCBI GEO); RNA-seq of sorafenib-resistant versus parental Huh7 cells
- Differential expression via limma-style moderated *t*-statistic
- Volcano plot with thresholds: |log₂FC| > 1 and adjusted *p* < 0.05

**Outputs:**

| File | Description |
|------|-------------|
| `deg_full_results.csv` | Full differential expression table (12,998 genes) |
| `deg_significant.csv` | Significant DEGs only (|log₂FC| > 1, adj. *p* < 0.05) |
| `volcano_plot.png/svg` | Volcano plot with labelled top hits |

---

### `03_ppi_network/` — Protein-Protein Interaction Network (Preliminary)

**Purpose:** Initial PPI network construction for SOCE and EMT proteins from the STRING database.

> **Note:** See `17_ppi_improved/` for the final, publication-ready network used in Figure 7E.

**Outputs:** `ppi_edges.csv`, `ppi_node_stats.csv`, `ppi_network.png/svg`

---

### `04_gsea/` — Gene Set Enrichment Analysis (Preliminary)

**Purpose:** Initial GSEA run on GSE140202 DEGs against Hallmark gene sets.

> **Note:** See `14_tcga_gsea_stim1/` for the final version with extended pathway panel used in Figure 7D.

**Outputs:** `gsea_barplot.png/svg`

---

### `05_coexpression/` — Spearman Co-expression (Preliminary)

**Purpose:** Initial Spearman correlation analysis between SOCE components and EMT-related genes in TCGA-LIHC.

> **Note:** See `15_clean_correlation/` for the final heatmap used in Figure 7C.

**Outputs:** `coexpression_results.csv`, `correlation_heatmap.png/svg`, `scatter_plots.png/svg`

---

### `06_stim1_trpc6_correlation/` — Direct STIM1–TRPC6 Correlation

**Purpose:** Characterise the direct co-expression relationship between STIM1 and TRPC6 across the TCGA-LIHC cohort.

**Outputs:** `stim1_trpc6_analysis.png/svg` — scatter plot with Spearman ρ and regression line

---

### `07_cox_regression/` — Multivariate Cox Proportional Hazards

**Purpose:** Test whether STIM1 and TRPC6 are independent prognostic factors after adjustment for established clinical confounders.

**Models fitted:**

| Model | Covariates | n |
|-------|-----------|---|
| Model 1 | SOCE gene + Stage + Grade + Age + Sex | 344 |
| Model 2 | Model 1 + AFP (log₂) | 264 |

**Key finding:** Neither STIM1 (*p* = 0.711) nor TRPC6 (*p* = 0.314) reached significance in multivariate analysis. **Tumour stage** is the dominant independent predictor (HR = 1.68, *p* < 0.001). This finding is reported in the manuscript Discussion as an important limitation of univariate survival analyses.

**Outputs:**

| File | Description |
|------|-------------|
| `cox_results.csv` | Full coefficient table (HR, 95% CI, *p*-value) for both models |
| `cox_forest_plot.png/svg` | Forest plot of HR and confidence intervals |

---

### `08_depmap_drug_sensitivity/` — DepMap Drug Sensitivity

**Purpose:** Investigate whether SOCE gene expression correlates with drug sensitivity scores in HCC cell lines using DepMap PRISM compound profiling data.

**Outputs:** `depmap_correlation_results.csv`, `depmap_sensitivity_scatter.png/svg`

---

### `09_soce_risk_score/` — Composite SOCE Risk Score

**Purpose:** Construct a composite SOCE expression score and evaluate its prognostic value in TCGA-LIHC.

**Method:** Principal component analysis (PCA) on z-score normalised expression of STIM1, TRPC6, TRPC1, and ORAI1. PC1 used as composite risk score. Patients split at median into risk-high and risk-low groups for KM analysis.

**Key finding:** Dual TRPC6-high + STIM1-high group showed a trend towards improved OS (median: 21.6 vs. 18.3 months; log-rank *p* = 0.057, borderline significant).

**Outputs:**

| File | Description |
|------|-------------|
| `soce_risk_score_patients.csv` | Per-patient risk scores and survival data |
| `soce_risk_score_summary.csv` | Group-level summary statistics |
| `soce_risk_score_km.png/svg` | KM curves: risk-high vs. risk-low |
| `soce_risk_score_km_v2.png/svg` | Alternative stratification |

---

### `10_trpc6_inhibitors/` — TRPC6 Inhibitor Landscape

**Purpose:** Catalogue known small-molecule TRPC6 inhibitors from ChEMBL to contextualise the translational relevance of pharmacological TRPC6 targeting in HCC.

**Method:** ChEMBL REST API query for TRPC6 (UniProt Q9Y210) bioactive compounds; filtered by activity type (IC₅₀, Ki) and assay confidence score ≥ 7.

**Outputs:**

| File | Description |
|------|-------------|
| `chembl_trpc6_inhibitors.csv` | All bioactive ChEMBL hits |
| `chembl_trpc6_specific_inhibitors.csv` | High-confidence TRPC6-selective compounds |
| `trpc6_inhibitor_summary.csv` | Summary by pharmacological class |
| `trpc6_inhibitor_landscape.png/svg` | Inhibitor activity landscape |

---

### `11_dual_high_clinical/` — Dual-High Patient Clinical Characteristics

**Purpose:** Characterise clinical and pathological features of TCGA-LIHC patients with simultaneous high STIM1 and high TRPC6 expression.

**Outputs:** `clinical_characteristics.csv`, `clinical_characteristics_table.png/svg`

---

### `12_crispr_essentiality/` — CRISPR Gene Essentiality (DepMap)

**Purpose:** Evaluate whether SOCE genes — particularly TRPC6 — are essential for HCC cell viability using CRISPR-Cas9 loss-of-function screening data from the Cancer Dependency Map.

**Method:** DepMap gene effect scores (Chronos algorithm). Threshold: gene effect < −0.5 indicates essentiality (benchmarked against common essential genes such as *RPL14*, *PCNA*). Non-essential genes typically score near 0.

**Key finding:** TRPC6 mean gene effect across 21 HCC cell lines = **−0.028** (range: −0.204 to +0.082). All 21 lines are above the −0.5 threshold → **non-essential** in all tested HCC models. This supports TRPC6 as a candidate therapeutic target with an inherently favourable safety margin.

**Outputs:**

| File | Description |
|------|-------------|
| `trpc6_gene_effect_hcc.csv` | Per-cell-line DepMap gene effect scores with cell line names |
| `crispr_essentiality.png/svg` | Gene effect distribution across HCC lines |

---

### `13_stage_stratified_km/` — Stage-Stratified Kaplan-Meier

**Purpose:** Assess whether the STIM1/TRPC6 univariate survival signal persists after stratification by pathological tumour stage (Stage I–II vs. Stage III–IV).

**Key finding:** Univariate KM associations weaken substantially within stage-stratified subgroups, consistent with stage acting as the primary confounder of SOCE expression–OS associations.

**Outputs:** `stage_km_results.csv`, `stage_stratified_km.png/svg`

---

### `14_tcga_gsea_stim1/` — GSEA: Sorafenib Resistance Pathways (Final)

**Purpose:** Identify oncogenic signalling pathways enriched in sorafenib-resistant Huh7 cells using preranked GSEA. This is the **final version** used in manuscript Figure 7D.

**Method:**
- Gene list: 12,998 genes ranked by limma *t*-statistic (resistant − sensitive; GSE140202)
- Gene sets: MSigDB Hallmark + custom KEGG/literature-curated gene sets (9 pathways total)
- Tool: `gseapy` v1.1.13; 1,000 permutations; preranked mode
- Significance threshold: FDR *q* < 0.25

**Key finding:** TNFα/NF-κB signalling is the top enriched pathway in sorafenib-resistant cells (NES = 1.85, FDR < 0.001), independently corroborating the STIM1–RELA co-expression axis identified in TCGA-LIHC.

**Outputs:**

| File | Description |
|------|-------------|
| `gsea_geo140202_results.csv` | Full GSEA results (NES, FDR, leading-edge genes) |
| `gsea_geo140202_barplot.png/svg` | Ranked NES barplot (Fig. 7D) |
| `tcga_stim1_pathway_correlation.csv` | Mean Spearman ρ between STIM1 and pathway gene members in TCGA |
| `tcga_stim1_pathway_enrichment.png/svg` | Supplementary pathway correlation figure |
| `figure_gsea_combined.png/svg` | Combined GSEA + TCGA pathway panel |

---

### `15_clean_correlation/` — Spearman Co-expression Heatmap (Final)

**Purpose:** Final, publication-ready clustered Spearman co-expression heatmap of STIM1 and TRPC6 versus EMT markers, transcription factors, and drug resistance genes in TCGA-LIHC. Used in manuscript **Figure 7C**.

**Method:**
- Spearman rank correlation for each SOCE–partner gene pair in tumour samples (n = 370)
- Significance threshold: ρ ≥ 0.2 and *p* < 0.01 (not adjusted for multiple testing; reported as exploratory)
- Hierarchical clustering (Ward linkage) on correlation matrix
- Asterisks (*) mark *p* < 0.01 pairs

**Outputs:**

| File | Description |
|------|-------------|
| `soce_emtf_corr_matrix.csv` | Full ρ matrix (wide format) |
| `soce_partner_correlations.csv` | Long-format table: SOCE gene, partner, ρ, *p*-value |
| `soce_correlation_heatmap_final.png/svg` | Final heatmap (Fig. 7C) |
| `soce_correlation_heatmap_clean.png/svg` | Clean version (alternative layout) |
| `soce_correlation_heatmap_v2.png/svg` | Version with adjusted colour scale |

---

### `16_emtf_correlation/` — EMT Transcription Factor Bubble Plot

**Purpose:** Visualise the magnitude and significance of correlations between SOCE components and EMT transcription factors as a bubble/dot plot for supplementary figures.

**Outputs:** `emtf_correlation_data.csv`, `emtf_bubble_dotplot.png/svg` (dot size = −log₁₀ *p*; colour = ρ)

---

### `17_ppi_improved/` — PPI Network with Ca²⁺-TF Bridge Nodes (Final)

**Purpose:** Final protein-protein interaction network integrating SOCE components, EMT transcription factors, and Ca²⁺-dependent transcription factor bridge nodes (NFATC1, NFKB1, RELA). Used in manuscript **Figure 7E**.

**Method:**
- STRING v12 database; combined interaction score ≥ 0.7
- Ca²⁺-dependent bridge nodes (NFATC1, NFKB1, RELA) added from literature curation (orange edges)
- Visualised with `networkx` v3.4.2; layout: spring algorithm
- Node colour encodes functional group; edge thickness encodes STRING confidence score

**Outputs:**

| File | Description |
|------|-------------|
| `ppi_node_stats_improved.csv` | Degree centrality, betweenness centrality, functional annotation |
| `ppi_network_improved.png/svg` | Final network figure (Fig. 7E) |
| `ppi_network_v2.png/svg` | Alternative layout |

---

### `18_gse14520_validation/` — Independent Cohort Validation (GSE14520)

**Purpose:** Validate the key co-expression and overexpression findings from TCGA-LIHC in an independent HCC dataset.

**Dataset:** Roessler et al. (Zhongshan Hospital, Shanghai); n = 247 HCC tumours + 40 paired non-tumour liver samples; Affymetrix HG-U133A 2.0 Array (GPL571).

> ⚠️ **This script requires internet access** to automatically download raw GEO data (~200 MB SOFT file). Run time: ~10–15 minutes.

**Analyses:**
1. SOCE gene overexpression: tumour vs. non-tumour (Mann-Whitney U, two-sided)
2. Kaplan-Meier overall survival for STIM1 and TRPC6 (split at median expression)
3. Spearman co-expression: STIM1 and TRPC6 vs. EMT gene panel
4. Cross-cohort comparison table: TCGA-LIHC ρ vs. GSE14520 ρ for 9 key gene pairs

**Running the validation:**

```bash
python 18_gse14520_validation/gse14520_validation.py
```

**Expected outputs (generated upon running):**

| File | Description |
|------|-------------|
| `tumor_vs_normal_results.csv` | Mann-Whitney statistics for SOCE genes |
| `tumor_vs_normal_boxplot.png/svg` | Tumour vs. non-tumour boxplots |
| `km_survival_gse14520.png/svg` | Kaplan-Meier curves in GSE14520 |
| `km_survival_results.csv` | Log-rank *p* and HR for STIM1, TRPC6 |
| `spearman_heatmap.png/svg` | Co-expression heatmap in GSE14520 |
| `spearman_correlations.csv` | Full Spearman ρ table |
| `tcga_vs_gse14520_comparison.csv` | Cross-cohort validation table |
| `validation_summary.txt` | Plain-text summary of validation results |

---

## Reproducing the Analyses

### 1. Environment Setup

```bash
# Clone this repository
git clone https://github.com/<USERNAME>/hcc-soce-bioinformatics.git
cd hcc-soce-bioinformatics

# Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate      # Linux / macOS
# venv\Scripts\activate       # Windows

# Install all dependencies
pip install -r requirements.txt
```

### 2. Download Primary Data

Most analysis scripts load pre-computed results from the CSV files in this repository. To run analyses from scratch, download the following:

**TCGA-LIHC** (cBioPortal):
```
https://www.cbioportal.org → lihc_tcga → Download RNA-seq v2 RSEM + clinical data
```

**GTEx Liver** (GTEx Portal):
```
https://gtexportal.org/home/datasets
→ GTEx_Analysis_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
→ Filter for tissue = "Liver"
```

**GSE140202** (NCBI GEO) — downloaded automatically by the analysis scripts via `GEOparse`.

**GSE14520** (NCBI GEO) — downloaded automatically by `18_gse14520_validation/gse14520_validation.py`.

### 3. Running Modules

```bash
# Core analyses (run in this order for reproducibility)
python 01_tcga_survival/tcga_survival_analysis.py
python 02_geo_deg/geo_deg_analysis.py
python 07_cox_regression/cox_regression.py
python 12_crispr_essentiality/crispr_essentiality.py
python 14_tcga_gsea_stim1/gsea_analysis.py
python 15_clean_correlation/clean_correlation.py
python 17_ppi_improved/ppi_network.py

# Independent validation (requires internet access)
python 18_gse14520_validation/gse14520_validation.py
```

---

## Software Requirements

All analyses were performed in **Python 3.10** on Ubuntu 22.04 LTS.

```
pandas >= 2.0
numpy >= 1.24
scipy >= 1.11
matplotlib >= 3.7
seaborn >= 0.12
lifelines == 0.30.0
gseapy == 1.1.13
networkx == 3.4.2
GEOparse >= 2.0.4
scikit-learn >= 1.3
```

Install via pip:

```bash
pip install pandas numpy scipy matplotlib seaborn "lifelines==0.30.0" "gseapy==1.1.13" "networkx==3.4.2" GEOparse scikit-learn
```

---

## Statistical Methods Summary

| Analysis | Test / Method | Library | Threshold |
|----------|--------------|---------|-----------|
| Tumour vs. normal | Mann-Whitney U (two-sided) | `scipy.stats` | *p* < 0.001 |
| Overall survival (univariate) | Kaplan-Meier + log-rank | `lifelines` v0.30.0 | *p* < 0.05 |
| Hazard ratio (univariate) | Cox proportional hazards | `lifelines` v0.30.0 | *p* < 0.05 |
| Hazard ratio (multivariate) | Multivariate Cox PH | `lifelines` v0.30.0 | *p* < 0.05 |
| Gene co-expression | Spearman rank correlation | `scipy.stats` | ρ ≥ 0.2, *p* < 0.01 |
| Pathway enrichment | Preranked GSEA | `gseapy` v1.1.13 | FDR *q* < 0.25 |
| PPI network | STRING + networkx | `networkx` v3.4.2 | Score ≥ 0.7 |
| Gene essentiality | Chronos (DepMap) | — | Gene effect < −0.5 |

---

## Manuscript Figures

| Figure | Module(s) | Description |
|--------|-----------|-------------|
| **Fig. 7A** | `01_tcga_survival/` | SOCE overexpression in HCC vs. normal (boxplots) |
| **Fig. 7B** | `01_tcga_survival/` | Kaplan-Meier overall survival (STIM1, TRPC6) |
| **Fig. 7C** | `15_clean_correlation/` | Spearman co-expression heatmap |
| **Fig. 7D** | `14_tcga_gsea_stim1/` | GSEA barplot (sorafenib resistance pathways) |
| **Fig. 7E** | `17_ppi_improved/` | PPI network with Ca²⁺-TF bridge nodes |
| **Supplementary** | `07_cox_regression/`, `12_crispr_essentiality/` | Cox forest plot, CRISPR essentiality |

The combined 6-panel figure used in the manuscript is provided as `FIGURE_comprehensive_6panel.png`.

---

## Authors

**Beste Yurdacan Yasar**
Department of Pharmacology, Faculty of Pharmacy, Istinye University, Istanbul, Turkey
ORCID: [0000-0003-0937-4394](https://orcid.org/0000-0003-0937-4394)

**Yasemin Erac** *(Corresponding Author)*
Department of Pharmacology, Faculty of Pharmacy, Ege University, Izmir, Turkey
ORCID: [0000-0002-3522-7921](https://orcid.org/0000-0002-3522-7921)
✉ yasemin.erac@ege.edu.tr

---

## Funding

This study was supported by the Ege University Scientific Research Projects Coordination Unit (Project No: TDK-2021-23075) and by The Scientific and Technological Research Council of Türkiye (TÜBİTAK) through the 1002 Short-Term R&D Funding Program (Project No: 223S649).

---

## Data Availability

All raw data used in this study are publicly available:

| Resource | URL |
|----------|-----|
| TCGA-LIHC | https://www.cbioportal.org (study: `lihc_tcga`) |
| GTEx Portal | https://gtexportal.org/home/datasets |
| GEO GSE140202 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140202 |
| GEO GSE14520 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520 |
| DepMap Portal | https://depmap.org/portal/ |
| STRING | https://string-db.org |
| ChEMBL | https://www.ebi.ac.uk/chembl |

---

## Citation

If you use this code or data in your own research, please cite:

```
Yurdacan Yasar B, Erac Y. Store-Operated Calcium Entry Regulates Early Adaptive
Responses to Sorafenib in Hepatocellular Carcinoma. 2026 (under review).
```

---

## License

This repository is released under the [MIT License](LICENSE). The underlying biological data are subject to the terms of use of the respective public repositories (TCGA data use policy, GTEx data use agreement, NCBI GEO terms of service, DepMap terms of use).
