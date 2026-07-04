# HCC SOCE Manuscript Analysis Roadmap

This document defines the bioinformatics analyses required to complete the manuscript story as a professional, publication-oriented workflow. It also maps each priority analysis to the existing reproducible code in this repository.

## Main Story We Need to Support

The manuscript argues that SOCE contributes to early adaptive responses to sorafenib in HCC, promotes EMT-associated plasticity, and may be therapeutically targetable. To support that story convincingly, the bioinformatics layer must answer five questions:

1. Are SOCE genes dysregulated in human HCC?
2. Are SOCE genes linked to EMT, partial/hybrid EMT, and plasticity rather than only single-marker changes?
3. Is the signal clinically relevant across prognosis, subtype, and tumour context?
4. Is the SOCE program linked to sorafenib adaptation/resistance and downstream pathway activity?
5. Are the findings robust after confounder control and independent validation?

## Priority Order

### Tier 1: Required Core Analyses

These analyses are the minimum set needed to make the manuscript story coherent and defensible.

1. **Tumour vs. normal SOCE expression**
   - Biological question: Are `TRPC1`, `TRPC6`, `STIM1`, and `ORAI1` upregulated in HCC?
   - Why it matters: Establishes clinical relevance beyond the Huh7 model.
   - Existing code: [01_tcga_survival/tcga_survival.py](01_tcga_survival/tcga_survival.py)
   - Outputs: `01_tcga_survival/tumor_vs_normal_stats.csv`, `tumor_vs_normal_boxplot.png/svg`

2. **Survival analysis of SOCE genes**
   - Biological question: Do SOCE genes associate with patient outcome?
   - Why it matters: Shows clinical consequence of SOCE dysregulation.
   - Existing code: [01_tcga_survival/tcga_survival.py](01_tcga_survival/tcga_survival.py)
   - Outputs: `01_tcga_survival/km_survival_results.csv`, `km_survival_curves.png/svg`
   - Interpretation note: Use as univariate evidence only.

3. **Multivariable Cox regression**
   - Biological question: Are SOCE associations independent of stage and other clinical factors?
   - Why it matters: Prevents overclaiming prognostic value.
   - Existing code: [07_cox_regression/cox_regression.py](07_cox_regression/cox_regression.py)
   - Outputs: `07_cox_regression/cox_results.csv`, `cox_forest_plot.png/svg`

4. **Expanded SOCE-EMT correlation analysis**
   - Biological question: Is SOCE linked to EMT and EMT-transcription-factor programs, not just `CDH1/CDH2/VIM`?
   - Why it matters: The current manuscript under-captures EMT plasticity.
   - Existing code:
   - [15_clean_correlation/clean_correlation.py](15_clean_correlation/clean_correlation.py)
   - [16_emtf_correlation/emtf_correlation.py](16_emtf_correlation/emtf_correlation.py)
   - Outputs:
   - `15_clean_correlation/soce_partner_correlations.csv`
   - `15_clean_correlation/soce_emtf_corr_matrix.csv`
   - `16_emtf_correlation/emtf_correlation_data.csv`
   - Recommendation: This should replace the current 7-gene TCGA panel in the manuscript main text.

5. **Hybrid / partial EMT or EMP scoring**
   - Biological question: Does the SOCE-high state track with epithelial-mesenchymal plasticity rather than a simple binary EMT model?
   - Why it matters: The manuscript currently claims hybrid E/M states, but the existing text is too weak for that claim.
   - Existing code: [19_emp_score/emp_score.py](19_emp_score/emp_score.py)
   - Outputs:
   - `19_emp_score/emp_scores.csv`
   - `19_emp_score/emp_score_distribution.png/svg`
   - `19_emp_score/emp_soce_boxplot.png/svg`
   - `19_emp_score/emp_km_survival.png/svg`

6. **Sorafenib-resistance transcriptomics**
   - Biological question: What transcriptomic program characterizes sorafenib adaptation/resistance in Huh7?
   - Why it matters: Connects cell-line pharmacology to a broader resistance program.
   - Existing code: [02_geo_deg/geo_deg.py](02_geo_deg/geo_deg.py)
   - Outputs: `02_geo_deg/deg_full_results.csv`, `deg_significant.csv`, `volcano_plot.png/svg`

7. **Pathway enrichment on resistance data**
   - Biological question: Are EMT, NF-kB, IL6/JAK/STAT, TGF-beta, calcium, apoptosis, and adaptation pathways enriched?
   - Why it matters: Supports the mechanism-level story instead of relying on single genes.
   - Existing code:
   - [04_gsea/gsea.py](04_gsea/gsea.py)
   - [14_tcga_gsea_stim1/tcga_gsea_stim1.py](14_tcga_gsea_stim1/tcga_gsea_stim1.py)
   - Outputs:
   - `04_gsea/gsea_barplot.png/svg`
   - `14_tcga_gsea_stim1/tcga_stim1_pathway_correlation.csv`
   - `14_tcga_gsea_stim1/tcga_stim1_pathway_enrichment.png/svg`
   - Recommendation: Use the extended pathway version as the primary figure.

8. **Independent cohort validation**
   - Biological question: Do major SOCE findings reproduce outside TCGA?
   - Why it matters: Validation is essential for a publishable human-data layer.
   - Existing code: [18_gse14520_validation/gse14520_validation.py](18_gse14520_validation/gse14520_validation.py)
   - Output folder: `18_gse14520_validation/`

### Tier 2: Strongly Recommended Analyses

These analyses materially strengthen the manuscript and address likely reviewer questions.

9. **Tumour purity / stromal confounding control**
   - Biological question: Are SOCE-EMT associations still present after purity adjustment?
   - Why it matters: Bulk RNA-seq hybrid EMT signals are vulnerable to stromal admixture.
   - Existing code: [21_tumor_purity/tumor_purity.py](21_tumor_purity/tumor_purity.py)
   - Outputs:
   - `21_tumor_purity/purity_scores.csv`
   - `21_tumor_purity/partial_correlations.csv`
   - `21_tumor_purity/correlation_comparison.png/svg`

10. **Stage-stratified survival**
   - Biological question: Is prognostic value stage-specific?
   - Why it matters: Helps reconcile univariate and multivariate findings.
   - Existing code: [13_stage_stratified_km/stage_stratified_km.py](13_stage_stratified_km/stage_stratified_km.py)
   - Outputs: `13_stage_stratified_km/stage_km_results.csv`, `stage_stratified_km.png/svg`

11. **SOCE composite subgrouping / dual-high classification**
   - Biological question: Does a combined SOCE state identify a clinically meaningful subgroup?
   - Why it matters: Better reflects pathway state than a single-gene readout.
   - Existing code:
   - [09_soce_risk_score/soce_risk_score.py](09_soce_risk_score/soce_risk_score.py)
   - [11_dual_high_clinical/dual_high_clinical.py](11_dual_high_clinical/dual_high_clinical.py)
   - Outputs:
   - `09_soce_risk_score/soce_risk_score_summary.csv`
   - `09_soce_risk_score/soce_risk_score_patients.csv`
   - `11_dual_high_clinical/clinical_characteristics.csv`

12. **Mutation / copy-number context**
   - Biological question: Are SOCE alterations transcription-driven only, or linked to genomic alterations?
   - Why it matters: Adds mechanistic context and translational depth.
   - Existing code: [20_mutation_cna/mutation_cna.py](20_mutation_cna/mutation_cna.py)
   - Outputs: `20_mutation_cna/alteration_summary.csv`, `mutation_summary.csv`, `cna_summary.csv`

13. **Molecular subtype analysis**
   - Biological question: Do SOCE-high tumours enrich in known HCC molecular subtypes?
   - Why it matters: Places the findings within established HCC biology.
   - Existing code: [22_hoshida_subtype/hoshida_subtype.py](22_hoshida_subtype/hoshida_subtype.py)
   - Outputs: `22_hoshida_subtype/subtype_assignments.csv`, `soce_by_subtype.png/svg`, `hoshida_km_survival.png/svg`

### Tier 3: Translational / Supportive Analyses

These are not strictly required to complete the mechanistic manuscript, but they strengthen translational framing.

14. **Drug sensitivity associations in HCC cell lines**
   - Existing code: [08_depmap_drug_sensitivity/depmap_drug_sensitivity.py](08_depmap_drug_sensitivity/depmap_drug_sensitivity.py)

15. **CRISPR essentiality**
   - Existing code: [12_crispr_essentiality/crispr_essentiality.py](12_crispr_essentiality/crispr_essentiality.py)

16. **PPI / network reconstruction**
   - Existing code:
   - [17_ppi_improved/ppi_improved.py](17_ppi_improved/ppi_improved.py)
   - [03_ppi_network/ppi_network.py](03_ppi_network/ppi_network.py)

17. **TRPC6 inhibitor landscape**
   - Existing code: [10_trpc6_inhibitors/trpc6_inhibitors.py](10_trpc6_inhibitors/trpc6_inhibitors.py)

## What Is Still Missing or Needs Tightening

Even with the existing scripts, these points need careful handling in the final manuscript:

1. **Do not overstate sorafenib-treated TCGA survival.**
   - The old manuscript text uses a very small `n=29` sorafenib-treated subgroup. Keep this exploratory at most.

2. **Do not claim MDR1 / NFAT / NF-kB / CREB as demonstrated unless directly shown.**
   - If those appear in the graphical abstract, support them with explicit analyses from GSE140202/TCGA or tone the figure down.

3. **Use EMP/hybrid EMT scoring instead of inferring hybrid state from `CDH1`-`VIM` correlation alone.**

4. **Use purity-controlled results whenever discussing hybrid EMT from bulk tumours.**

5. **Lead with pathway- and program-level findings, not isolated gene-gene correlations.**

## Recommended Figure Logic for the Manuscript

If we rebuild the bioinformatics figure set, the cleanest narrative order is:

1. Tumour vs normal SOCE expression
2. Expanded SOCE-EMT / EMT-TF heatmap
3. EMP score association with SOCE
4. Sorafenib-resistance DEG + pathway enrichment
5. Purity-adjusted confirmation
6. Independent cohort validation
7. Survival / subtype / dual-high clinical context

## Publication-Safe Bottom Line

If we want a conservative, reviewer-resistant manuscript, the essential evidence stack should be:

1. `01_tcga_survival`
2. `07_cox_regression`
3. `15_clean_correlation`
4. `16_emtf_correlation`
5. `19_emp_score`
6. `02_geo_deg`
7. `14_tcga_gsea_stim1`
8. `21_tumor_purity`
9. `18_gse14520_validation`

That set is the strongest core package for the SOCE-EMT-sorafenib adaptation story.
