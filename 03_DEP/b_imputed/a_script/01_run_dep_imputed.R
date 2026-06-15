#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/b_imputed  --  COMPARISON Mito DEP on imputed matrices.
# imp4p = CANONICAL (full proteoDA tables + plots); mscoreutils = comparison
# (tables only). logFC compared to the non-imputed primary.
# CAVEAT: exploratory only -- a_non_imputed is primary (impute-before-test can
# inflate false positives).
# =============================================================================

suppressPackageStartupMessages({ library(here); library(readr); library(dplyr); library(purrr); library(tibble) })
source(here("00_input", "h9c2_design.R"))
source(here("R", "dep_model.R"))

CANONICAL <- "imp4p"
nd <- here("02_Normalization", "c_data")
methods <- c(imp4p = "DAList_imputed_imp4p.rds", mscoreutils = "DAList_imputed_mscoreutils.rds")

runs <- imap(methods, function(rds, m) {
  dal <- readRDS(file.path(nd, rds))
  cat(sprintf("\n[%s] imputed DEP: %d x %d%s\n", m, nrow(dal$data), ncol(dal$data),
              if (m == CANONICAL) " (canonical)" else ""))
  run_dep_model(dal, H9C2_DESIGN_FORMULA, H9C2_CONTRASTS,
                out_dir = here("03_DEP", "b_imputed", "c_data", m),
                grouping_column = "group",
                report_dir = if (m == CANONICAL) here("03_DEP", "b_imputed", "b_reports", m) else NULL,
                pval_thresh = H9C2_PVAL_THRESH, pi_thresh = H9C2_PI_THRESH)
})

ni_file <- here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv")
if (file.exists(ni_file)) {
  ni <- read_csv(ni_file, show_col_types = FALSE) |> select(uniprot_id, contrast, logFC_ni = logFC)
  cmp <- imap_dfr(runs, function(out, m) {
    bind_rows(out$results) |> select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |> group_by(contrast) |>
      summarise(method = m,
                spearman_vs_nonimputed = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
                .groups = "drop")
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "logFC_vs_nonimputed.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman):\n"); print(as.data.frame(cmp))
}
