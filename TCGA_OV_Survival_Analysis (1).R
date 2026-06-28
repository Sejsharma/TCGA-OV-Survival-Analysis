# =============================================================================
# TCGA Ovarian Cancer (OV) - Survival Analysis
# Author: Sejal Sharma
# Dataset: TCGA-OV (The Cancer Genome Atlas - Ovarian Cancer)
# Packages: TCGAbiolinks, survival, survminer, ggplot2, dplyr
# Description: Kaplan-Meier survival analysis correlating TP53 gene expression
#              with overall survival in ovarian cancer patients
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 1: Install and load required packages
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE)
install.packages(c("survival", "survminer", "ggplot2", "dplyr", "tidyr"), 
                 repos = "https://cran.r-project.org")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# STEP 2: Download TCGA-OV clinical data
# -----------------------------------------------------------------------------

cat("Downloading TCGA-OV clinical data...\n")

clinical <- GDCquery_clinic(project = "TCGA-OV", type = "clinical")

# Preview
cat("Clinical data dimensions:", dim(clinical), "\n")
head(clinical[, c("submitter_id", "vital_status", "days_to_death", 
                  "days_to_last_follow_up", "age_at_index")])

# -----------------------------------------------------------------------------
# STEP 3: Download TCGA-OV gene expression data (RNA-seq)
# -----------------------------------------------------------------------------

cat("Querying TCGA-OV gene expression data...\n")

query <- GDCquery(
  project      = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type    = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query)

cat("Expression data dimensions:", dim(assay(data)), "\n")

# -----------------------------------------------------------------------------
# STEP 4: Extract TP53 expression
# (TP53 is the most commonly mutated gene in ovarian cancer ~96% of cases)
# -----------------------------------------------------------------------------

# Get expression matrix (log2 TPM or counts — we'll use log-normalised counts)
expr_matrix <- assay(data, "tpm_unstrand")

# Log2 transform
expr_log2 <- log2(expr_matrix + 1)

# Find TP53 row (gene symbol in rowData)
gene_info <- rowData(data)
tp53_id   <- rownames(gene_info)[gene_info$gene_name == "TP53"]

if (length(tp53_id) == 0) stop("TP53 not found in dataset")

tp53_expr <- as.numeric(expr_log2[tp53_id, ])

# Patient barcodes (truncate to 12 characters to match clinical)
patient_ids <- substr(colnames(expr_log2), 1, 12)

expr_df <- data.frame(
  submitter_id = patient_ids,
  TP53_expr    = tp53_expr,
  stringsAsFactors = FALSE
)

# Remove duplicate patient entries (keep first)
expr_df <- expr_df[!duplicated(expr_df$submitter_id), ]

cat("Unique patients with expression data:", nrow(expr_df), "\n")

# -----------------------------------------------------------------------------
# STEP 5: Prepare survival data
# -----------------------------------------------------------------------------

surv_df <- clinical %>%
  select(submitter_id, vital_status, days_to_death, days_to_last_follow_up,
         age_at_index, tumor_stage) %>%
  mutate(
    # Overall survival time: use days_to_death if dead, else days_to_last_follow_up
    OS_time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    # Event: 1 = dead, 0 = censored (alive)
    OS_event = ifelse(vital_status == "Dead", 1, 0),
    # Convert days to months
    OS_months = OS_time / 30.44
  ) %>%
  filter(!is.na(OS_time), OS_time > 0)

cat("Patients with valid survival data:", nrow(surv_df), "\n")

# -----------------------------------------------------------------------------
# STEP 6: Merge expression and survival data
# -----------------------------------------------------------------------------

merged_df <- inner_join(surv_df, expr_df, by = "submitter_id")

cat("Patients in merged dataset:", nrow(merged_df), "\n")

# Dichotomise TP53 expression at median (High vs Low)
median_tp53 <- median(merged_df$TP53_expr, na.rm = TRUE)

merged_df <- merged_df %>%
  mutate(
    TP53_group = ifelse(TP53_expr >= median_tp53, "High", "Low"),
    TP53_group = factor(TP53_group, levels = c("Low", "High"))
  )

cat("TP53 group distribution:\n")
print(table(merged_df$TP53_group))

# -----------------------------------------------------------------------------
# STEP 7: Kaplan-Meier Survival Analysis
# -----------------------------------------------------------------------------

surv_obj <- Surv(time = merged_df$OS_months, event = merged_df$OS_event)

km_fit <- survfit(surv_obj ~ TP53_group, data = merged_df)

# Log-rank test
log_rank <- survdiff(surv_obj ~ TP53_group, data = merged_df)
p_value  <- 1 - pchisq(log_rank$chisq, df = 1)

cat("\nLog-rank test p-value:", round(p_value, 4), "\n")

# -----------------------------------------------------------------------------
# STEP 8: Plot Kaplan-Meier curves
# -----------------------------------------------------------------------------

km_plot <- ggsurvplot(
  km_fit,
  data          = merged_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.col = "strata",
  palette       = c("#2196F3", "#F44336"),   # Blue = Low, Red = High
  legend.labs   = c("TP53 Low", "TP53 High"),
  legend.title  = "TP53 Expression",
  title         = "Overall Survival by TP53 Expression\nTCGA Ovarian Cancer (OV)",
  xlab          = "Time (Months)",
  ylab          = "Survival Probability",
  ggtheme       = theme_bw(base_size = 13),
  font.main     = c(14, "bold"),
  surv.median.line = "hv",
  risk.table.height = 0.25
)

# Save plot
ggsave(
  filename = "KM_TP53_OvarianCancer.png",
  plot     = km_plot$plot,
  width    = 10,
  height   = 7,
  dpi      = 300
)

cat("KM plot saved as KM_TP53_OvarianCancer.png\n")

# -----------------------------------------------------------------------------
# STEP 9: Cox Proportional Hazards Model (multivariate)
# (adjusting for age — adds statistical depth to your analysis)
# -----------------------------------------------------------------------------

cox_model <- coxph(
  Surv(OS_months, OS_event) ~ TP53_group + age_at_index,
  data = merged_df
)

cat("\n--- Cox Proportional Hazards Model ---\n")
print(summary(cox_model))

# Forest plot of Cox model
cox_plot <- ggforest(
  cox_model,
  data      = merged_df,
  main      = "Hazard Ratios — TCGA Ovarian Cancer",
  fontsize  = 1.0
)

ggsave(
  filename = "CoxModel_TP53_OvarianCancer.png",
  plot     = cox_plot,
  width    = 10,
  height   = 5,
  dpi      = 300
)

cat("Cox forest plot saved as CoxModel_TP53_OvarianCancer.png\n")

# -----------------------------------------------------------------------------
# STEP 10: Save merged data for reproducibility
# -----------------------------------------------------------------------------

write.csv(merged_df, "TCGA_OV_survival_data.csv", row.names = FALSE)
cat("Data saved as TCGA_OV_survival_data.csv\n")

cat("\n✅ Analysis complete! Files generated:\n")
cat("  - KM_TP53_OvarianCancer.png\n")
cat("  - CoxModel_TP53_OvarianCancer.png\n")
cat("  - TCGA_OV_survival_data.csv\n")
