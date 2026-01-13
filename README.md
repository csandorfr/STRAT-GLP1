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

rmarkdown::render("reports/power_report.Rmd")
The report will be saved to reports/power_report.html.

Inputs

Edit config/params.yml to update cohort sizes, subtype fractions, MAF grid, alpha thresholds,
expected event counts, and LEDD slope assumptions.

Notes / assumptions
	•	GWAS uses standard effective sample size + allelic test approximation.
	•	Cox power uses Schoenfeld event-based approximation.
	•	LEDD slope power is a planning approximation; replace with empirical slope SD if available.
	•	MR power requires instrument strength (R²) from cis-eQTL/pQTL instruments.

Contact

Cynthia Sandor lb / Webber group collaboration.
