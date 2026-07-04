# Reproducibility Runbook

This runbook records where the code lives, in what order the analyses should be run, and which outputs each step produces. It is the operational companion to [ANALYSIS_ROADMAP.md](ANALYSIS_ROADMAP.md).

## Environment

### Python

- Version target: `Python 3.11`
- Install dependencies:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### R

- Required for the GEO differential expression workflow if the script uses `DESeq2`
- Suggested packages:

```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Code Inventory

| Step | Module | Script | Main output folder |
|------|--------|--------|--------------------|
| 1 | Tumour vs normal + KM | `01_tcga_survival/tcga_survival.py` | `01_tcga_survival/` |
| 2 | GEO differential expression | `02_geo_deg/geo_deg.py` | `02_geo_deg/` |
| 3 | Preliminary PPI | `03_ppi_network/ppi_network.py` | `03_ppi_network/` |
| 4 | Preliminary GSEA | `04_gsea/gsea.py` | `04_gsea/` |
| 5 | Preliminary coexpression | `05_coexpression/coexpression.py` | `05_coexpression/` |
| 6 | STIM1-TRPC6 relationship | `06_stim1_trpc6_correlation/stim1_trpc6_correlation.py` | `06_stim1_trpc6_correlation/` |
| 7 | Multivariable Cox | `07_cox_regression/cox_regression.py` | `07_cox_regression/` |
| 8 | Drug sensitivity | `08_depmap_drug_sensitivity/depmap_drug_sensitivity.py` | `08_depmap_drug_sensitivity/` |
| 9 | SOCE risk score | `09_soce_risk_score/soce_risk_score.py` | `09_soce_risk_score/` |
| 10 | TRPC6 inhibitor landscape | `10_trpc6_inhibitors/trpc6_inhibitors.py` | `10_trpc6_inhibitors/` |
| 11 | Dual-high clinical profile | `11_dual_high_clinical/dual_high_clinical.py` | `11_dual_high_clinical/` |
| 12 | CRISPR essentiality | `12_crispr_essentiality/crispr_essentiality.py` | `12_crispr_essentiality/` |
| 13 | Stage-stratified KM | `13_stage_stratified_km/stage_stratified_km.py` | `13_stage_stratified_km/` |
| 14 | Final pathway enrichment | `14_tcga_gsea_stim1/tcga_gsea_stim1.py` | `14_tcga_gsea_stim1/` |
| 15 | Final SOCE-EMT correlation | `15_clean_correlation/clean_correlation.py` | `15_clean_correlation/` |
| 16 | EMT transcription factor panel | `16_emtf_correlation/emtf_correlation.py` | `16_emtf_correlation/` |
| 17 | Final PPI | `17_ppi_improved/ppi_improved.py` | `17_ppi_improved/` |
| 18 | External validation | `18_gse14520_validation/gse14520_validation.py` | `18_gse14520_validation/` |
| 19 | EMP / hybrid EMT scoring | `19_emp_score/emp_score.py` | `19_emp_score/` |
| 20 | Mutation / CNA | `20_mutation_cna/mutation_cna.py` | `20_mutation_cna/` |
| 21 | Tumour purity adjustment | `21_tumor_purity/tumor_purity.py` | `21_tumor_purity/` |
| 22 | Hoshida subtype | `22_hoshida_subtype/hoshida_subtype.py` | `22_hoshida_subtype/` |

## Recommended Execution Order

Run the analyses in this order when rebuilding the manuscript evidence from scratch:

### Core manuscript build

```bash
python 01_tcga_survival/tcga_survival.py
python 07_cox_regression/cox_regression.py
python 15_clean_correlation/clean_correlation.py
python 16_emtf_correlation/emtf_correlation.py
python 19_emp_score/emp_score.py
python 02_geo_deg/geo_deg.py
python 14_tcga_gsea_stim1/tcga_gsea_stim1.py
python 21_tumor_purity/tumor_purity.py
python 18_gse14520_validation/gse14520_validation.py
```

### Extended clinical and translational build

```bash
python 13_stage_stratified_km/stage_stratified_km.py
python 09_soce_risk_score/soce_risk_score.py
python 11_dual_high_clinical/dual_high_clinical.py
python 20_mutation_cna/mutation_cna.py
python 22_hoshida_subtype/hoshida_subtype.py
python 08_depmap_drug_sensitivity/depmap_drug_sensitivity.py
python 12_crispr_essentiality/crispr_essentiality.py
python 17_ppi_improved/ppi_improved.py
python 10_trpc6_inhibitors/trpc6_inhibitors.py
```

## Which Results Should Be Used in the Manuscript

Use these modules as the main evidence base:

1. `01_tcga_survival`
2. `07_cox_regression`
3. `15_clean_correlation`
4. `16_emtf_correlation`
5. `19_emp_score`
6. `02_geo_deg`
7. `14_tcga_gsea_stim1`
8. `21_tumor_purity`
9. `18_gse14520_validation`

Treat these as supporting analyses:

1. `09_soce_risk_score`
2. `11_dual_high_clinical`
3. `13_stage_stratified_km`
4. `20_mutation_cna`
5. `22_hoshida_subtype`
6. `08_depmap_drug_sensitivity`
7. `12_crispr_essentiality`
8. `17_ppi_improved`
9. `10_trpc6_inhibitors`

## File Preservation Rules

- Keep each module self-contained: script, tables, figures, and any intermediate files should stay in the same numbered folder.
- Do not overwrite summary files without regenerating the associated figures and tables.
- Commit both code and output files together for each module so that each result remains traceable.

## Git Notes

- The current repository contains untracked analysis outputs; review before staging.
- The current [push_to_github.sh](push_to_github.sh) is a safe status/helper script. It does not reset history, force-push, or stage the whole tree.
- Do not stage local manuscript drafts, `tmp/`, `.venv/`, or raw GEO cache files.
- When we are ready to publish, use a safer staged workflow: review, stage, commit, push, then open the GitHub PR.

## Suggested Commit Strategy

Use small commits grouped by analysis family:

1. `docs: add manuscript analysis roadmap and reproducibility runbook`
2. `analysis: finalize core SOCE EMT modules`
3. `analysis: add validation and purity-controlled results`
4. `analysis: add extended clinical and translational modules`

## Final Check Before GitHub Push

Before pushing:

1. Confirm that each primary figure is backed by a script in the same module.
2. Confirm that the manuscript only claims what the code directly demonstrates.
3. Confirm that exploratory analyses are labeled exploratory.
4. Confirm that no temporary files from `tmp/` are accidentally staged unless intentionally needed.

## Current External-Source Caveats

- `02_geo_deg` uses the checked-in `deg_full_results.csv` as a fallback when the GSE140202 SOFT record does not expose sample-level expression tables.
- `08_depmap_drug_sensitivity` and `12_crispr_essentiality` can run from local DepMap raw files; otherwise they generate deterministic manuscript-summary fallback outputs and say so in the console.
- `18_gse14520_validation` downloads GSE14520 through HTTPS because GEO FTP size checks can fail in GEOparse. The current SOFT metadata lacks OS columns, so KM survival is skipped while tumor-normal and SOCE/EMT validation outputs are generated.
