# Data Sources, Tools, and Methods
## SOCE Signaling in HCC — Comprehensive Bioinformatics Analysis

**Date:** February 2026


---

## 1. Databases and Data Sources

### 1.1 Clinical Genomics

| Database | Version / Release | Data Used | Access URL | Citation |
|----------|------------------|-----------|------------|---------|
| **cBioPortal** | TCGA-LIHC (2024 freeze) | RNA-seq RSEM expression (n=371 tumors), overall survival (OS) months/status, clinical variables (stage, grade, AFP, age, sex, vascular invasion) | https://www.cbioportal.org | Cerami et al., Cancer Discov 2012; Gao et al., Sci Signal 2013 |
| **GTEx** | v8 (dbGaP phs000424.v8) | Normal liver RNA-seq TPM (n=226 donors) | https://gtexportal.org | GTEx Consortium, Science 2020 |
| **NCBI GEO** | GSE140202 | Sorafenib-resistant vs. sensitive HuH-7 RNA-seq (n=3 replicates/group, raw counts) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140202 | Hagiwara et al., 2021 |

### 1.2 Protein Interaction & Network

| Database | Version | Data Used | Access URL |
|----------|---------|-----------|------------|
| **STRING** | v11.5 | Protein-protein interactions (confidence score ≥ 0.7), network edges and scores for 61 seed + partner genes | https://string-db.org |

### 1.3 Functional Genomics (CRISPR)

| Database | Release | Data Used | Access URL | DOI / release record |
|----------|---------|-----------|------------|-------------|
| **DepMap** | 24Q4 (December 2024) | CRISPRGeneEffect.csv — Chronos gene effect scores (1,178 cell lines × 17,916 genes) | https://depmap.org/portal | doi:10.25452/figshare.plus.27993248.v1 |
| **DepMap** | 25Q3 (2025) | Model and omics/drug-response mapping files used for HCC cell-line drug sensitivity integration | https://depmap.org/portal/download/all/ | https://forum.depmap.org/t/announcing-the-25q3-release/4476 |

### 1.4 Drug Sensitivity

| Database | Version | Data Used | Access URL |
|----------|---------|-----------|------------|
| **GDSC (Genomics of Drug Sensitivity in Cancer)** | GDSC2 (Sanger Institute) | Drug viability scores for sorafenib and thapsigargin across HCC cell lines; COSMIC cell line IDs for mapping | https://www.cancerrxgene.org |
| **DepMap** | 25Q3 | TPM expression data for HCC cell lines (OmicsExpressionProteinCodingGenesTPMLogp1.csv) and cell-line metadata | https://depmap.org/portal/download/all/ |

### 1.5 Drug / Chemical

| Database | Version | Data Used | Access URL |
|----------|---------|-----------|------------|
| **ChEMBL** | v33 (REST API) | TRPC6 inhibitor activity records: 278 unique compounds, 844 activity records across targets CHEMBL2417347, CHEMBL1795081, CHEMBL4523582 | https://www.ebi.ac.uk/chembl/ |

### 1.6 Clinical Trials

| Database | Access Date | Query Terms | URL |
|----------|------------|-------------|-----|
| **ClinicalTrials.gov** | February 2026 | TRPC6, STIM1, ORAI1, SOCE, store-operated calcium, HCC, hepatocellular carcinoma, liver cancer | https://clinicaltrials.gov |

---

## 2. Software Tools and Versions

### 2.1 Programming Languages

| Language | Version |
|----------|---------|
| Python | 3.11.8 |
| R | 4.4.0 (for DESeq2) |

### 2.2 Python Libraries

| Library | Version | Used For |
|---------|---------|---------|
| pandas | 2.1.4 | Data manipulation, merging, filtering |
| numpy | 1.26.3 | Numerical operations, array handling |
| scipy | 1.12.0 | Statistical tests (Mann-Whitney U, Spearman correlation, chi-squared, Fisher's exact) |
| matplotlib | 3.8.2 | All figure generation and visualization |
| seaborn | 0.13.2 | Correlation heatmaps |
| lifelines | 0.27.8 | Kaplan-Meier survival analysis, log-rank test, Cox proportional hazards regression |
| gseapy | 1.1.3 | Gene Set Enrichment Analysis (GSEA prerank) |
| networkx | 3.2.1 | PPI network construction and centrality analysis |
| GEOparse | 2.0.4 | GEO dataset download and parsing |
| requests | 2.31.0 | REST API queries (cBioPortal, STRING, ChEMBL, GTEx) |
| Pillow (PIL) | 10.2.0 | Figure caption overlay processing |

### 2.3 R / Bioconductor Packages

| Package | Version | Used For |
|---------|---------|---------|
| DESeq2 | 1.42.0 | Differential expression analysis (GSE140202) |
| ggplot2 | 3.5.0 | Supplementary visualization |

### 2.4 External Tools / APIs

| Tool | Version / Endpoint | Used For |
|------|-------------------|---------|
| cBioPortal REST API | v3 (2024) | Fetching TCGA-LIHC expression, clinical, and survival data |
| GTEx REST API | v8 | Fetching normal liver expression data |
| STRING REST API | v11.5 | Fetching PPI network data |
| ChEMBL REST API | v33 | Fetching TRPC6 inhibitor activity data |
| DepMap portal / Figshare+ | 24Q4 / 25Q3 | Downloading public DepMap release files |
| ClinicalTrials.gov API | v2 | Searching clinical trials by gene/condition |
| NCBI GEO FTP | — | Downloading GSE140202 raw count matrix |

---

## 3. Analysis Methods

### 3.1 Tumor vs. Normal Expression (Section 1)
- **Input:** TCGA-LIHC RNA-seq RSEM (log₂+1 transformed) + GTEx v8 liver TPM (log₂+1 transformed)
- **Test:** Two-sided Mann-Whitney U test (scipy.stats.mannwhitneyu)
- **Genes:** TRPC1, TRPC6, STIM1, ORAI1
- **Visualization:** Box plots with jittered individual points (matplotlib)

### 3.2 Kaplan-Meier Survival Analysis (Sections 2, 9, 13)
- **Input:** TCGA-LIHC OS months + OS status (DECEASED/LIVING)
- **Stratification:** Median expression split (high ≥ median, low < median)
- **Test:** Log-rank test (lifelines.statistics.logrank_test)
- **Estimator:** Kaplan-Meier (lifelines.KaplanMeierFitter)
- **Dual-High definition:** TRPC6 ≥ median AND STIM1 ≥ median
- **Stage stratification:** AJCC pathologic tumor stage grouped as Stage I-II vs. Stage III-IV

### 3.3 Differential Expression Analysis (Section 3)
- **Input:** GEO GSE140202 raw count matrix (sorafenib-resistant HuH7-SR vs. sensitive HuH7, n=3/group)
- **Method:** DESeq2 (R 4.4.0, Bioconductor 3.18)
- **Significance threshold:** FDR (Benjamini-Hochberg) ≤ 0.05 AND |log₂FC| ≥ 0.5
- **Total DEGs identified:** 4,936 (2,288 up, 2,648 down)

### 3.4 PPI Network Analysis (Section 4)
- **Input:** 8 seed genes (TRPC1, TRPC6, STIM1, ORAI1, VIM, CDH1, CDH2, ABCB1)
- **Database:** STRING v11.5, confidence score ≥ 0.7 (high confidence)
- **Network metrics:** Degree centrality, betweenness centrality (networkx 3.2)
- **Final network:** 61 nodes, 75 edges

### 3.5 Gene Set Enrichment Analysis (Section 5)
- **Input:** Ranked gene list from GSE140202 (ranking metric: log₂FC × −log₁₀(p-value))
- **Gene sets:** MSigDB Hallmark (v2023.2) + KEGG (via gseapy built-in)
- **Method:** GSEApy prerank, 1,000 permutations
- **Significance threshold:** FDR ≤ 0.25 (standard GSEA threshold)

### 3.6 Co-expression Analysis (Section 6)
- **Input:** TCGA-LIHC RNA-seq RSEM (log₂+1), n=371 tumor samples
- **Method:** Spearman rank correlation (scipy.stats.spearmanr)
- **Gene panel:** 40 genes across SOCE, EMT, MDR, signaling, and apoptosis categories

### 3.7 Multivariate Cox Regression (Section 8)
- **Input:** TCGA-LIHC expression + clinical covariates (stage, grade, age, sex, AFP)
- **Method:** Cox proportional hazards (lifelines.CoxPHFitter, Breslow baseline)
- **Model 1:** n=344 (no AFP); **Model 2:** n=264 (with AFP)
- **Evaluation:** Concordance index (C-index)

### 3.8 Drug Sensitivity Correlation (Section 10)
- **Input:** DepMap 24Q4 TPM expression + GDSC2 drug viability scores
- **Cell lines:** 21 HCC lines (OncotreePrimaryDisease = "Hepatocellular Carcinoma")
- **Drugs:** Sorafenib (n=13 lines), thapsigargin (n=15 lines)
- **Method:** Spearman rank correlation (scipy.stats.spearmanr)

### 3.9 CRISPR Essentiality Analysis (Section 14)
- **Input:** DepMap 24Q4 CRISPRGeneEffect.csv (Chronos algorithm)
- **Cell lines:** 21 HCC lines with CRISPR data
- **Thresholds:** Selective essential < −0.5; common essential < −1.0 (DepMap standard)
- **Comparison:** HCC lines vs. pan-cancer (n=1,178) using Mann-Whitney U test

### 3.10 Clinical Characteristics Table (Section 12)
- **Input:** TCGA-LIHC clinical data via cBioPortal
- **Groups:** Dual-High (TRPC6 ≥ median AND STIM1 ≥ median, n=99) vs. Other (n=272)
- **Tests:** Mann-Whitney U (continuous), χ²/Fisher's exact (categorical)

---

## 4. Figure Index

| Figure | File | Analysis | Key Finding |
|--------|------|---------|-------------|
| Fig 1 | Fig01_tumor_vs_normal.png | Tumor vs. normal expression | All SOCE genes overexpressed (p<10⁻⁹²) |
| Fig 2 | Fig02_km_survival.png | Kaplan-Meier OS | STIM1 p=0.014, TRPC6 p=0.033 |
| Fig 3 | Fig03_volcano_deg.png | DESeq2 DEG | TRPC6↑ in resistant cells (padj=0.0001) |
| Fig 4 | Fig04_ppi_network.png | STRING PPI network | TRPC6 highest betweenness centrality |
| Fig 5 | Fig05_gsea.png | GSEA | IL-17 NES=2.03, Wnt suppressed |
| Fig 6A | Fig06a_coexpression_heatmap.png | Spearman correlation | TRPC6-VIM r=0.484, STIM1-TRPC6 r=0.043 |
| Fig 6B | Fig06b_coexpression_scatter.png | Co-expression scatter | TRPC6-ZEB1/ZEB2 EMT axis |
| Fig 7 | Fig07_stim1_trpc6_relationship.png | STIM1-TRPC6 relationship | Independent regulation, channel-switching model |
| Fig 8 | Fig08_cox_forest_plot.png | Multivariate Cox regression | Stage only independent factor (HR=1.678) |
| Fig 9 | Fig09_depmap_drug_sensitivity.png | DepMap drug sensitivity | TRPC6-sorafenib trend r=−0.454 |
| Fig 10 | Fig10_soce_risk_score_km.png | SOCE risk score KM | Dual-High p=0.025, Δ55.9 months |
| Fig 11 | Fig11_trpc6_inhibitor_landscape.png | ChEMBL inhibitor landscape | SAR7334 IC₅₀=7.9 nM, no clinical trials |
| Table 1 | Fig12_dual_high_clinical_table.png | Clinical characteristics | No confounding by stage/grade/AFP |
| Fig 13 | Fig13_crispr_essentiality.png | DepMap CRISPR essentiality | TRPC6 non-essential in all 21 HCC lines |
| Fig 14 | Fig14_stage_stratified_km.png | Stage-stratified KM | Prognostic only in Stage I-II (p<0.05) |

---

## 5. Reproducibility Notes

- All analyses were performed in a containerized Python 3.11 environment
- Random seeds set to 42 where applicable (jitter plots, permutation tests)
- All API queries used public, freely accessible endpoints (no authentication required)
- Raw data files are available from the original databases listed in Section 1
- All intermediate CSV files are saved in the corresponding numbered result folders


---

## 6. Key Literature References

| # | Reference | Relevance |
|---|-----------|-----------|
| 1 | Cerami et al., Cancer Discov 2012 | cBioPortal platform |
| 2 | GTEx Consortium, Science 2020 | GTEx v8 normal tissue data |
| 3 | Qiao et al., Future Med Chem 2023 | BP3112 — first HCC-specific TRPC6 inhibitor |
| 4 | Maier et al., Br J Pharmacol 2015 | SAR7334 — IC₅₀=7.9 nM |
| 5 | Lin et al., PNAS 2019 | BI-749327 — IC₅₀=13 nM, in vivo efficacy |
| 6 | Xu et al., Cancer Res 2018 | TRPC6/NCX1 complex in HCC invasion |
| 7 | Wen et al., Sci Rep 2016 | TRPC6 in HCC multi-drug resistance reversal |
| 8 | Bai et al., eLife 2020 | TRPC6 cryo-EM structure |
| 9 | Mukhopadhyay et al., Cell Rep 2023 | TRPC6 in breast cancer chemoresistance |
| 10 | Dimitrov et al., Cell Rep 2025 | TRPC6 in ferroptosis and metastasis |
| 11 | He et al., Biomolecules 2024 | ORAI Ca²⁺ channels in cancer review |
| 12 | Moccia et al., Curr Med Chem 2024 | Targeting STIM/ORAI in anticancer therapy |
| 13 | Tsherniak et al., Cell 2017 | DepMap CRISPR essentiality framework |
| 14 | Hagiwara et al., 2021 | GSE140202 — sorafenib resistance transcriptomics |
