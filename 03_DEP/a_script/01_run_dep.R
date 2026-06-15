#!/usr/bin/env Rscript
# =============================================================================
# 01_run_dep.R  --  Mito Stage 03: Differential abundance (native limma)
#
# One fit on the NON-imputed normalized matrix:
#   ~ 0 + group + (1 | Replicate)   group in {Ctl, Mito, PHE, PHE_Mito}
# proteoDA routes (1|Replicate) through duplicateCorrelation (paired plate/passage;
# Smyth 2005); eBayes(robust=TRUE). 5 contrasts (H9C2_CONTRASTS): CTLvPHE (disease),
# CTLvMITO (intervention), PHEvPHE_MITO (rescue), Interaction, MITOvPHE_MITO.
# Significance: p<0.10, BH-FDR, Pi-score<0.05 (Xiao 2014).
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble); library(purrr)
})
set.seed(42)
source(here("00_input", "h9c2_design.R"))   # H9C2_DESIGN_FORMULA, H9C2_CONTRASTS, H9C2_*_THRESH
source(here("R", "pi_score.R"))

data_dir <- here("03_DEP", "c_data")
dir.create(file.path(data_dir, "per_contrast"), recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "03_DAList_normalized.rds"))
cat(sprintf("Loaded normalized DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

dal <- add_design(dal, H9C2_DESIGN_FORMULA)
dal <- add_contrasts(dal, contrasts_vector = H9C2_CONTRASTS)
dal <- fit_limma_model(dal)
cat(sprintf("Replicate-block duplicateCorrelation rho = %.3f\n", dal$eBayes_fit$correlation %||% NA_real_))

dal <- extract_DA_results(dal, pval_thresh = H9C2_PVAL_THRESH, lfc_thresh = 0, adj_method = "BH")
ann <- as_tibble(dal$annotation) |> select(uniprot_id, gene, protein, description)
res <- compute_pi_scores(dal$results, pi_thresh = H9C2_PI_THRESH) |> map(~ left_join(.x, ann, by = "uniprot_id"))

iwalk(res, ~ write_csv(.x, file.path(data_dir, "per_contrast", paste0(.y, ".csv"))))
write_csv(bind_rows(res), file.path(data_dir, "03_combined_results.csv"))
summary_tbl <- map_dfr(res, ~ tibble(
  contrast = .x$contrast[1], n = nrow(.x),
  sig_p10 = sum(.x$P.Value < 0.10, na.rm = TRUE),
  sig_fdr10 = sum(.x$adj.P.Val < 0.10, na.rm = TRUE),
  sig_fdr05 = sum(.x$adj.P.Val < 0.05, na.rm = TRUE),
  up_pi = sum(.x$sig_pi == 1), down_pi = sum(.x$sig_pi == -1)))
write_csv(summary_tbl, file.path(data_dir, "02_DA_summary.csv"))
saveRDS(dal, file.path(data_dir, "01_DAList_dep.rds"))
cat("\nDA summary:\n"); print(as.data.frame(summary_tbl))
cat(sprintf("\nDone -> %s/\n", data_dir))
