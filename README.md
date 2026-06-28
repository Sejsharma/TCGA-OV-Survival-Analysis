# TCGA Ovarian Cancer — Survival Analysis

## Overview

This project investigates the relationship between **TP53 gene expression** and **overall survival** in ovarian cancer patients using data from The Cancer Genome Atlas (TCGA-OV).

TP53 is mutated in approximately 96% of high-grade serous ovarian cancers, making it the most clinically significant molecular marker in this disease. This analysis explores whether variation in TP53 expression level — even within a predominantly mutated cohort — correlates with patient outcomes.

---

## Methods

### Data Source
- **Dataset:** TCGA-OV (The Cancer Genome Atlas — Ovarian Serous Cystadenocarcinoma)
- **Access:** Downloaded via the `TCGAbiolinks` R/Bioconductor package
- **Data types:** RNA-seq gene expression (STAR counts, TPM) + clinical data

### Analysis Pipeline

1. **Data retrieval** — Clinical and RNA-seq data downloaded from GDC portal via TCGAbiolinks
2. **Expression extraction** — TP53 TPM values extracted and log2-transformed
3. **Patient stratification** — Patients split into High/Low TP53 expression groups at the median
4. **Survival analysis** — Kaplan-Meier curves with log-rank test for group comparison
5. **Multivariate modelling** — Cox proportional hazards model adjusting for age at diagnosis

### Tools & Packages

| Tool | Purpose |
|------|---------|
| R (v4.x) | Primary analysis language |
| TCGAbiolinks | TCGA data download and preparation |
| survival | Survival object creation and modelling |
| survminer | Kaplan-Meier visualisation |
| ggplot2 | Publication-quality plotting |
| dplyr | Data wrangling |

---

## Key Outputs

### 1. Kaplan-Meier Survival Curve
![KM Plot](KM_TP53_OvarianCancer.png)

Kaplan-Meier curves comparing overall survival between patients with high vs low TP53 expression. The log-rank p-value indicates whether the difference between groups is statistically significant.

### 2. Cox Proportional Hazards Forest Plot
![Cox Plot](CoxModel_TP53_OvarianCancer.png)

Multivariate Cox model adjusting for age at diagnosis, showing hazard ratios and 95% confidence intervals for each variable.

---

## How to Run

### Prerequisites
- R version 4.0 or above
- Internet connection (for TCGA data download)

### Installation & Execution

```r
# Clone this repository
# Open TCGA_OV_Survival_Analysis.R in RStudio

# The script will automatically install required packages on first run
# Then download TCGA-OV data (~500MB) and run the full analysis

# Expected runtime: 10-20 minutes (mostly data download)
```

### Output Files Generated
- `KM_TP53_OvarianCancer.png` — Kaplan-Meier survival plot
- `CoxModel_TP53_OvarianCancer.png` — Cox model forest plot
- `TCGA_OV_survival_data.csv` — Merged clinical + expression dataset

---

## Biological Context

Ovarian cancer has the highest mortality rate of all gynaecological cancers. High-grade serous ovarian carcinoma (HGSOC), the most common subtype, is characterised by near-universal TP53 mutation. Understanding how TP53 expression variation relates to clinical outcomes may inform prognostic stratification and treatment decisions.

This analysis forms part of an independent bioinformatics portfolio developed alongside an MSc in Biotechnology (University of Greenwich, 2025–2026).

---

## Author

**Sejal Sharma**
MSc Biotechnology, University of Greenwich
Member, British Association for Cancer Research

---

## License

This project uses publicly available TCGA data. All data is subject to the [TCGA Data Use Certification](https://gdc.cancer.gov/access-data/data-use-certification-agreement).
