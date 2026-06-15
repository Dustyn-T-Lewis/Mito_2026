#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/a_non_imputed  --  PRIMARY Mito DEP on the NON-imputed normalized matrix.
# limma handles NAs; no imputation before testing. ~0+group+(1|Replicate),
# 5 contrasts (H9C2_CONTRASTS). See R/dep_model.R.
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("00_input", "h9c2_design.R"))   # H9C2_DESIGN_FORMULA, H9C2_CONTRASTS, H9C2_*_THRESH
source(here("R", "dep_model.R"))

dal <- readRDS(here("02_Normalization", "c_data", "03_DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))
out <- run_dep_model(dal, H9C2_DESIGN_FORMULA, H9C2_CONTRASTS,
                     here("03_DEP", "a_non_imputed", "c_data"),
                     pval_thresh = H9C2_PVAL_THRESH, pi_thresh = H9C2_PI_THRESH)
cat(sprintf("duplicateCorrelation rho = %.3f\n", out$rho))
cat("\nDA summary (non-imputed):\n"); print(as.data.frame(out$summary))
