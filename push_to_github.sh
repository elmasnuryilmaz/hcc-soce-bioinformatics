#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

echo "=== HCC SOCE GitHub Publish Helper ==="
echo "Repository: $(git remote get-url origin 2>/dev/null || echo 'no origin configured')"
echo ""

if ! command -v git >/dev/null 2>&1; then
  echo "git is required but was not found."
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "This folder is not a git repository."
  exit 1
fi

echo "Current branch: $(git branch --show-current)"
echo ""
git status -sb
echo ""
echo "This helper does not reset history, force-push, or stage the whole tree."
echo "Review the status above, then stage intentional files explicitly, for example:"
echo ""
echo "  git add README.md CODE_AVAILABILITY.md REPRODUCIBILITY_RUNBOOK.md ANALYSIS_ROADMAP.md"
echo "  git add requirements.txt environment.yml CITATION.cff .gitignore"
echo "  git add 01_tcga_survival 02_geo_deg 07_cox_regression 08_depmap_drug_sensitivity"
echo "  git add 09_soce_risk_score 11_dual_high_clinical 12_crispr_essentiality"
echo "  git add 13_stage_stratified_km 15_clean_correlation 16_emtf_correlation"
echo "  git add 18_gse14520_validation 19_emp_score 20_mutation_cna 21_tumor_purity 22_hoshida_subtype"
echo ""
echo "Then commit and push:"
echo ""
echo "  git commit -m \"docs: prepare manuscript code availability package\""
echo "  git push origin $(git branch --show-current)"
