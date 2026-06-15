# =============================================================================
# dep_model.R  --  shared DEP runner built on proteoDA's native capabilities.
#
# One call fits the model and writes proteoDA's own outputs:
#   write_limma_tables()  -> c_data: summary + per-contrast CSVs + combined + xlsx
#   write_limma_plots()   -> b_reports: volcano / MD / p-value-hist per contrast
# Pi-score (Xiao 2014) is the only thing proteoDA doesn't emit, so we add a single
# tidy `combined_results_pi.csv`. Used by 03_DEP/{a_non_imputed,b_imputed}; the two
# arms differ only in the matrix passed in.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble); library(purrr)
})
source(here("R", "pi_score.R"))

run_dep_model <- function(dal, formula, contrasts_vec, out_dir, grouping_column,
                          report_dir = NULL, pval_thresh = 0.10, pi_thresh = 0.05) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  dal <- add_design(dal, formula)
  dal <- add_contrasts(dal, contrasts_vector = contrasts_vec)
  dal <- fit_limma_model(dal)
  dal <- extract_DA_results(dal, pval_thresh = pval_thresh, lfc_thresh = 0, adj_method = "BH")

  # proteoDA native result tables -> keep the consolidated multi-sheet results.xlsx
  # and drop the loose CSVs it also emits (redundant with the workbook).
  write_limma_tables(dal, output_dir = out_dir, overwrite = TRUE,
                     annot_cols = c("uniprot_id", "gene", "protein", "description"))
  unlink(c(file.path(out_dir, "combined_results.csv"), file.path(out_dir, "DA_summary.csv"),
           file.path(out_dir, "per_contrast_results")), recursive = TRUE)

  # proteoDA native diagnostic plots (volcano / MD / p-value histogram)
  if (!is.null(report_dir)) {
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    write_limma_plots(dal, grouping_column = grouping_column,
                      table_columns = c("uniprot_id", "gene"),
                      output_dir = report_dir, overwrite = TRUE)
  }

  # Pi-score augmentation (not part of proteoDA tables)
  ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- compute_pi_scores(dal$results, pi_thresh = pi_thresh) |> map(~ left_join(.x, ann, by = "uniprot_id"))
  write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))

  list(dal = dal, rho = dal$eBayes_fit$correlation %||% NA_real_, results = res)
}
