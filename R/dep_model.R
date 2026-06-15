# =============================================================================
# dep_model.R  --  shared limma DEP runner (used by 03_DEP/{a_non_imputed,b_imputed})
#
# Fits one proteoDA/limma model and writes per-contrast + combined + summary +
# the DAList. Identical for non-imputed and imputed inputs, so both DEP arms call
# this and differ only in the matrix they pass in.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble); library(purrr)
})
source(here("R", "pi_score.R"))

run_dep_model <- function(dal, formula, contrasts_vec, out_dir,
                          pval_thresh = 0.10, pi_thresh = 0.05) {
  dir.create(file.path(out_dir, "per_contrast"), recursive = TRUE, showWarnings = FALSE)
  dal <- add_design(dal, formula)
  dal <- add_contrasts(dal, contrasts_vector = contrasts_vec)
  dal <- fit_limma_model(dal)
  rho <- dal$eBayes_fit$correlation %||% NA_real_
  dal <- extract_DA_results(dal, pval_thresh = pval_thresh, lfc_thresh = 0, adj_method = "BH")

  ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- compute_pi_scores(dal$results, pi_thresh = pi_thresh) |> map(~ left_join(.x, ann, by = "uniprot_id"))
  iwalk(res, ~ write_csv(.x, file.path(out_dir, "per_contrast", paste0(.y, ".csv"))))
  write_csv(bind_rows(res), file.path(out_dir, "03_combined_results.csv"))
  summ <- map_dfr(res, ~ tibble(
    contrast = .x$contrast[1], n = nrow(.x),
    sig_p10 = sum(.x$P.Value < 0.10, na.rm = TRUE),
    sig_fdr10 = sum(.x$adj.P.Val < 0.10, na.rm = TRUE),
    sig_fdr05 = sum(.x$adj.P.Val < 0.05, na.rm = TRUE),
    up_pi = sum(.x$sig_pi == 1), down_pi = sum(.x$sig_pi == -1)))
  write_csv(summ, file.path(out_dir, "02_DA_summary.csv"))
  saveRDS(dal, file.path(out_dir, "01_DAList_dep.rds"))
  list(rho = rho, summary = summ, results = res)
}
