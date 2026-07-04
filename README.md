# TRPC6-STIM1 SOCE Bioinformatics in Hepatocellular Carcinoma

Companion analysis code for the manuscript:

**Prognostic Significance of TRPC6 and STIM1 in Hepatocellular Carcinoma: An Integrative Bioinformatic Analysis**

Authors: Beste Yurdacan Yasar, Elmasnur Yilmaz, Yasemin Erac

Repository: <https://github.com/elmasnuryilmaz/hcc-soce-bioinformatics>

## Overview

This repository contains the bioinformatics code, analysis outputs, and reproducibility notes supporting an integrative analysis of store-operated calcium entry (SOCE) genes in hepatocellular carcinoma (HCC). The main genes evaluated are `TRPC6`, `STIM1`, `TRPC1`, `ORAI1`, and related SOCE, EMT, MDR, signaling, apoptosis, and clinical covariate panels.

The analyses use public datasets from TCGA-LIHC, GTEx, cBioPortal, DepMap, GDSC, GEO, STRING, ChEMBL, and related resources. No protected patient-level data are included in this repository.

## Code Availability

All analysis scripts and reproducibility documentation are maintained in this GitHub repository:

<https://github.com/elmasnuryilmaz/hcc-soce-bioinformatics>

The repository is organized as numbered analysis modules. Each module contains the script, generated summary tables, and figure outputs required to trace the corresponding manuscript result. See also:

- [CODE_AVAILABILITY.md](CODE_AVAILABILITY.md)
- [REPRODUCIBILITY_RUNBOOK.md](REPRODUCIBILITY_RUNBOOK.md)
- [ANALYSIS_ROADMAP.md](ANALYSIS_ROADMAP.md)
- [data_sources_and_methods.md](data_sources_and_methods.md)

## Main Public Data Sources

| Dataset | Release / version | Use in analysis | Access |
|---|---:|---|---|
| TCGA-LIHC | cBioPortal / GDC | tumor expression, clinical covariates, survival | <https://www.cbioportal.org>, <https://portal.gdc.cancer.gov> |
| GTEx liver | v8 | normal liver expression comparison | <https://gtexportal.org> |
| DepMap | 24Q4 Public | CRISPR Chronos gene effect scores | <https://doi.org/10.25452/figshare.plus.27993248.v1> |
| DepMap | 25Q3 Public | HCC cell-line expression / drug-response integration | <https://depmap.org/portal/download/all/> |
| GDSC | GDSC2 | sorafenib and thapsigargin sensitivity context | <https://www.cancerrxgene.org> |
| GEO | GSE140202 | sorafenib-resistant Huh7 transcriptomics | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140202> |
| GEO | GSE14520 | independent HCC validation cohort | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520> |
| STRING | v12 / API | protein interaction network | <https://string-db.org> |
| ChEMBL | REST API | TRPC6 inhibitor landscape | <https://www.ebi.ac.uk/chembl/> |

DepMap 25Q3 is cited through the DepMap portal and official release notes because recent DepMap releases are distributed through the portal rather than as Figshare records.

## Repository Structure

```text
.
├── 01_tcga_survival/              # TCGA/GTEx tumor-normal expression and Kaplan-Meier analysis
├── 02_geo_deg/                    # GSE140202 differential expression
├── 03_ppi_network/                # preliminary STRING PPI network
├── 04_gsea/                       # preliminary GSEA
├── 05_coexpression/               # preliminary SOCE co-expression analysis
├── 06_stim1_trpc6_correlation/    # STIM1-TRPC6 direct expression relationship
├── 07_cox_regression/             # multivariable Cox proportional hazards models
├── 08_depmap_drug_sensitivity/    # DepMap/GDSC drug sensitivity correlations
├── 09_soce_risk_score/            # combined SOCE risk score analysis
├── 10_trpc6_inhibitors/           # ChEMBL TRPC6 inhibitor landscape
├── 11_dual_high_clinical/         # TRPC6/STIM1 dual-high clinical characteristics
├── 12_crispr_essentiality/        # DepMap CRISPR essentiality
├── 13_stage_stratified_km/        # stage-stratified Kaplan-Meier analysis
├── 14_tcga_gsea_stim1/            # final pathway enrichment panels
├── 15_clean_correlation/          # final SOCE correlation heatmap
├── 16_emtf_correlation/           # EMT transcription factor correlation plot
├── 17_ppi_improved/               # final SOCE/EMT/Ca2+ bridge PPI network
├── 18_gse14520_validation/        # independent GSE14520 validation workflow
├── 19_emp_score/                  # epithelial-mesenchymal plasticity scoring
├── 20_mutation_cna/               # mutation and copy-number context
├── 21_tumor_purity/               # tumor purity adjustment
├── 22_hoshida_subtype/            # Hoshida subtype analysis
├── requirements.txt
├── environment.yml
├── REPRODUCIBILITY_RUNBOOK.md
├── ANALYSIS_ROADMAP.md
├── CODE_AVAILABILITY.md
└── data_sources_and_methods.md
```

## Environment

The manuscript methods specify Python 3.11. The environment below was used for the final manuscript-facing analyses.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

An equivalent conda-style environment is provided in [environment.yml](environment.yml).

## Quick Reproduction

Run the core manuscript analyses from the repository root:

```bash
python 01_tcga_survival/tcga_survival.py
python 07_cox_regression/cox_regression.py
python 15_clean_correlation/clean_correlation.py
python 16_emtf_correlation/emtf_correlation.py
python 09_soce_risk_score/soce_risk_score.py
python 13_stage_stratified_km/stage_stratified_km.py
python 12_crispr_essentiality/crispr_essentiality.py
python 08_depmap_drug_sensitivity/depmap_drug_sensitivity.py
```

Run the extended validation and mechanistic analyses:

```bash
python 02_geo_deg/geo_deg.py
python 14_tcga_gsea_stim1/tcga_gsea_stim1.py
python 18_gse14520_validation/gse14520_validation.py
python 19_emp_score/emp_score.py
python 20_mutation_cna/mutation_cna.py
python 21_tumor_purity/tumor_purity.py
python 22_hoshida_subtype/hoshida_subtype.py
python 17_ppi_improved/ppi_improved.py
python 10_trpc6_inhibitors/trpc6_inhibitors.py
```

The GSE14520 workflow downloads GEO data and may take several minutes depending on network speed.

Reproducibility notes:

- `02_geo_deg/geo_deg.py` first attempts to read the GSE140202 GEO series matrix. If GEO does not expose sample-level expression tables through the SOFT record, it regenerates the volcano plot from the checked-in `02_geo_deg/deg_full_results.csv` table.
- `08_depmap_drug_sensitivity/depmap_drug_sensitivity.py` and `12_crispr_essentiality/crispr_essentiality.py` can use local DepMap raw files when supplied. Without those large files, they use small deterministic manuscript-summary fallbacks and state this in the console output.
- `18_gse14520_validation/gse14520_validation.py` downloads the GSE14520 SOFT file over HTTPS. In the current SOFT metadata, overall-survival fields are not present, so the script generates tumor-vs-normal and SOCE/EMT validation outputs and skips the KM panel.

## Module-To-Manuscript Map

| Manuscript result | Main module(s) | Primary output |
|---|---|---|
| SOCE expression in HCC vs normal liver | `01_tcga_survival/` | `tumor_vs_normal_boxplot.png`, `tumor_vs_normal_stats.csv` |
| Clinical characteristics by TRPC6/STIM1 dual-high status | `11_dual_high_clinical/` | `clinical_characteristics.csv`, `clinical_characteristics_table.png` |
| Co-expression landscape | `15_clean_correlation/`, `16_emtf_correlation/` | `soce_correlation_heatmap_final.png`, `emtf_bubble_dotplot.png` |
| STIM1-TRPC6 direct relationship | `06_stim1_trpc6_correlation/` | `stim1_trpc6_analysis.png` |
| SOCE risk score | `09_soce_risk_score/` | `soce_risk_score_km_v2.png`, `soce_risk_score_summary.csv` |
| Stage-stratified survival | `13_stage_stratified_km/` | `stage_stratified_km.png`, `stage_km_results.csv` |
| Multivariable Cox models | `07_cox_regression/` | `cox_forest_plot.png`, `cox_results.csv` |
| CRISPR essentiality | `12_crispr_essentiality/` | `crispr_essentiality.png`, `trpc6_gene_effect_hcc.csv` |
| Drug sensitivity | `08_depmap_drug_sensitivity/` | `depmap_sensitivity_scatter.png`, `depmap_correlation_results.csv` |
| Independent validation | `18_gse14520_validation/` | validation tables and figures generated by the script |

## Statistical Methods

| Analysis | Method | Python package |
|---|---|---|
| tumor vs normal expression | two-sided Mann-Whitney U test | `scipy` |
| survival curves | Kaplan-Meier and log-rank test | `lifelines` |
| multivariable survival models | Cox proportional hazards | `lifelines` |
| gene co-expression | Spearman rank correlation | `scipy` |
| pathway enrichment | preranked GSEA | `gseapy` |
| PPI network | STRING API and graph metrics | `networkx` |
| visualization | static publication figures | `matplotlib`, `seaborn` |

## Data And Privacy Notes

This project uses de-identified public data only. The repository should not include raw protected health information, institutional credentials, API tokens, local `.venv` directories, local Word manuscript drafts, or temporary render folders.

Large raw downloads can be regenerated from the public sources listed above. The GEO cache for GSE14520 is intentionally excluded through `.gitignore`.

## Citation

If you use this repository before a formal journal citation is available, cite the GitHub repository URL and the manuscript title above. A versioned archival DOI can be minted later through Zenodo if required by the target journal.
