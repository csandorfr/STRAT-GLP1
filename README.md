# Power Analysis for Metabolically Stratified Parkinson’s Disease Studies

This repository contains a reproducible power analysis framework supporting the MJFF proposal:
**“Metabolically stratified Parkinson’s disease: genetics, biomarkers, and GLP-1RA response.”**

The analyses quantify statistical power across three linked aims:
1. Metabolically stratified PD GWAS
2. Drug-target Mendelian randomisation (GLP1R → PD subtypes)
3. Real‑world GLP‑1RA effectiveness analyses in CPRD

---

## Contents

- `power_report.Rmd`  
  Fully reproducible R Markdown file generating all power calculations, tables, and figures.

- `README.md`  
  This document.

---

## Cohort Sizes Used

| Cohort | PD cases | Controls / at-risk | Notes |
|------|---------|-------------------|------|
| UK Biobank | ~4,000 | ~450,000 | Genotype + Olink |
| FinnGen | ~5,000 | ~300,000 | Registry-linked |
| All of Us | ~3,000 | ~200,000 | Ancestry-diverse |
| Tracking / OPDC | ~2,500 | – | Deep phenotyping |
| CPRD Aurum | ~500,000 | >10M | Drug exposure |
| CPRD GOLD | ~250,000 | >5M | Drug exposure |

Metabolic stratification assumes **30–50%** of PD cases fall into insulin‑resistant / metabolic subtypes.

---

## Power Scenarios Covered

### Aim 1 – Metabolically Stratified PD GWAS
- Case–control GWAS
- Genome‑wide significance: α = 5×10⁻⁸
- Detectable odds ratios:
  - OR ≈ 1.15–1.20 for common variants (MAF ≥ 0.2)
  - OR ≈ 1.25–1.30 for moderate-frequency variants (MAF ≥ 0.1)

### Aim 1 (Extension) – Drug‑target MR (GLP1R)
- cis‑eQTL and cis‑pQTL instruments
- Power to detect OR ≤ 0.85 in metabolically enriched PD subtypes
- Heterogeneity testing between responder vs non‑responder subtypes

### Aim 3 – CPRD GLP‑1RA Effectiveness
- Active‑comparator, new‑user design
- ≥80% power to detect:
  - HR ≤ 0.85 for incident PD
  - HR ≤ 0.90 for LEDD‑based progression
- Powered stratified analyses by HbA1c, TG:HDL, BMI, hypertension

---

## Key Conclusions

- Available cohort sizes provide **adequate power for all proposed analyses**
- Metabolic stratification remains well powered even under conservative assumptions
- CPRD enables uniquely powered real‑world validation not previously possible
- Results directly support biologically targeted trial enrichment strategies

---

## Software Requirements

- R ≥ 4.2
- Packages:
  - `genpwr`
  - `pwr`
  - `TwoSampleMR`
  - `ggplot2`
  - `dplyr`

---

## Intended Use

This repository supports:
- MJFF proposal review
- Power justification sections
- GitHub‑hosted reproducibility
- Extension to future metabolic or pharmacogenetic analyses

---

## Contact

Cynthia Sandor, DVM PhD  
Imperial College London  
UK Dementia Research Institute
