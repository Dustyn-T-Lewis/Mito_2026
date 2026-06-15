#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/a_non_imputed  --  PRIMARY Mito DEP on the NON-imputed normalized matrix.
# limma handles NAs (no impute-before-test). Full proteoDA reporting (tables ->
# c_data, plots -> b_reports). ~0+group+(1|Replicate), 5 contrasts (H9C2_CONTRASTS).
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("00_input", "h9c2_design.R"))   # H9C2_DESIGN_FORMULA, H9C2_CONTRASTS, H9C2_*_THRESH
source(here("R", "dep_model.R"))

dal <- readRDS(here("02_Normalization", "c_data", "03_DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

out <- run_dep_model(dal, H9C2_DESIGN_FORMULA, H9C2_CONTRASTS,
                     out_dir = here("03_DEP", "a_non_imputed", "c_data"),
                     grouping_column = "group",
                     report_dir = here("03_DEP", "a_non_imputed", "b_reports"),
                     pval_thresh = H9C2_PVAL_THRESH, pi_thresh = H9C2_PI_THRESH)
cat(sprintf("duplicateCorrelation rho = %.3f | proteoDA tables + plots written\n", out$rho))
