# MJFF GLP-1RA × PD: Power calculations (Aims 1–3)

Reproducible power calculations supporting:
- **Aim 1**: metabolically stratified PD GWAS (case-control; MDE OR at GWAS alpha)
- **Aim 1 add-on**: GLP1R drug-target MR placeholder (requires instrument R²)
- **Aim 2**: proteomics DE power (two-sample approximation; Bonferroni/FDR planning)
- **Aim 3**: CPRD RWE power for incident PD (Schoenfeld/Cox) and progression proxy (LEDD slope)

## Quick start
```r
# from repo root
install.packages("renv")
renv::init()
renv::restore()

# run power calculations
source("R/run_power.R")
