#!/usr/bin/env Rscript
# Mito COMPARISON DEP on the imputed matrices (exploratory; a_non_imputed stays the reported analysis).
# imp4p = canonical (full proteoDA tables + plots); mscoreutils = comparison (tables only).
# CAVEAT: impute-before-test can inflate false positives -> logFC is checked against the non-imputed fit.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)

CANONICAL <- "imp4p"
clear_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE) }
clear_dir(here("03_DEP", "b_imputed", "c_data")); clear_dir(here("03_DEP", "b_imputed", "b_reports"))
nd <- here("02_Normalization", "imputation", "c_data")
methods <- c(imp4p = "DAList_imputed_imp4p.rds", mscoreutils = "DAList_imputed_mscoreutils.rds")

#### DEP per imputed matrix ####
# Identical limma workflow to the primary arm, run once per imputed matrix.

runs <- imap(methods, function(rds, m) {
  dal <- readRDS(file.path(nd, rds))
  cat(sprintf("\n[%s] imputed DEP: %d x %d%s\n", m, nrow(dal$data), ncol(dal$data),
              if (m == CANONICAL) " (canonical)" else ""))
  out_dir <- here("03_DEP", "b_imputed", "c_data", m)

  dal <- add_design(dal, "~ 0 + group + (1 | Replicate)")
  dal <- add_contrasts(dal, contrasts_vector = c(
    "CTLvPHE       = PHE - Ctl",
    "CTLvMITO      = Mito - Ctl",
    "PHEvPHE_MITO  = PHE_Mito - PHE",
    "Interaction   = (PHE_Mito - PHE) - (Mito - Ctl)",
    "MITOvPHE_MITO = PHE_Mito - Mito"))
  dal <- fit_limma_model(dal)
  dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

  write_limma_tables(dal, output_dir = out_dir, overwrite = TRUE,
                     annot_cols = c("uniprot_id", "gene", "protein", "description"))
  unlink(c(file.path(out_dir, "combined_results.csv"), file.path(out_dir, "DA_summary.csv"),
           file.path(out_dir, "per_contrast_results")), recursive = TRUE)
  if (m == CANONICAL)
    write_limma_plots(dal, grouping_column = "group", table_columns = c("uniprot_id", "gene"),
                      output_dir = here("03_DEP", "b_imputed", "b_reports", m), overwrite = TRUE)

  ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- imap(dal$results, function(r, cname)
    as_tibble(r, rownames = "uniprot_id") |>
      mutate(pi_score = P.Value ^ abs(logFC),
             sig_pi = case_when(pi_score < 0.05 & logFC > 0 ~ 1L,
                                pi_score < 0.05 & logFC < 0 ~ -1L, TRUE ~ 0L),
             contrast = cname) |>
      left_join(ann, by = "uniprot_id"))
  write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))
  res
})

#### logFC concordance vs non-imputed ####
# Spearman per contrast against the primary fit; high rho = imputation didn't distort effects.

ni_file <- here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv")
if (file.exists(ni_file)) {
  ni <- read_csv(ni_file, show_col_types = FALSE) |> select(uniprot_id, contrast, logFC_ni = logFC)
  cmp <- imap_dfr(runs, function(res, m) {
    bind_rows(res) |> select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |> group_by(contrast) |>
      summarise(method = m,
                spearman_vs_nonimputed = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
                .groups = "drop")
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "logFC_vs_nonimputed.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman):\n"); print(as.data.frame(cmp))
}
