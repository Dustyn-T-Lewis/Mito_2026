#!/usr/bin/env Rscript
# Mito PRIMARY DEP on the NON-imputed normalized matrix (limma handles per-protein NAs).
#   ~ 0 + group + (1 | Replicate), group {Ctl, Mito, PHE, PHE_Mito}; 5 factorial contrasts + recovery.
#   Significance: Pi-score (Xiao 2014, Pi = P.Value^|logFC|) < 0.05 + BH-FDR.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)

out_dir <- here("03_DEP", "a_non_imputed", "c_data")
report_dir <- here("03_DEP", "a_non_imputed", "b_reports")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(out_dir)
clear_dir(report_dir)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

#### Fit limma model ####
# (1 | Replicate) is a paired block estimated via duplicateCorrelation. The consensus rho is
# low (~0.017); proteoDA passes the block to lmFit regardless, so gls.series recomputes rho and
# the fit stays blocked GLS + robust eBayes. rho echoed in the printout below.

dal <- add_design(dal, "~ 0 + group + (1 | Replicate)")
dal <- add_contrasts(dal, contrasts_vector = c(
  "CTLvPHE       = PHE - Ctl", # Disease: PHE stress remodelling
  "CTLvMITO      = Mito - Ctl", # Intervention: transplant alone
  "PHEvPHE_MITO  = PHE_Mito - PHE", # Rescue: transplant under stress
  "Interaction   = (PHE_Mito - PHE) - (Mito - Ctl)", # Interaction: PHE-dependence of transplant
  "MITOvPHE_MITO = PHE_Mito - Mito", # Secondary: PHE effect in transplanted cells
  "CTLvPHE_MITO  = PHE_Mito - Ctl" # Recovery: residual disease signature after rescue
))
dal <- fit_limma_model(dal)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

#### Tables and plots ####
# Native proteoDA outputs; drop the loose CSVs and keep the consolidated workbook.

write_limma_tables(dal,
  output_dir = out_dir, overwrite = TRUE,
  annot_cols = c("uniprot_id", "gene", "protein", "description")
)
unlink(c(
  file.path(out_dir, "combined_results.csv"), file.path(out_dir, "DA_summary.csv"),
  file.path(out_dir, "per_contrast_results")
), recursive = TRUE)
write_limma_plots(dal,
  grouping_column = "group", table_columns = c("uniprot_id", "gene"),
  output_dir = report_dir, overwrite = TRUE
)

#### Pi-score ####
# Xiao 2014 pi on the p-scale: Pi = P.Value^|logFC| = 10^(-pi_Xiao); sig_pi = +1 up / -1 down / 0 ns at Pi < 0.05.

ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
res <- imap(dal$results, function(r, cname) {
  as_tibble(r, rownames = "uniprot_id") |>
    mutate(
      pi_score = P.Value^abs(logFC),
      sig_pi = case_when(
        pi_score < 0.05 & logFC > 0 ~ 1L,
        pi_score < 0.05 & logFC < 0 ~ -1L, TRUE ~ 0L
      ),
      contrast = cname
    ) |>
    left_join(ann, by = "uniprot_id")
})
write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))

cat(sprintf(
  "duplicateCorrelation rho = %.3f | proteoDA tables + plots + pi written\n",
  dal$eBayes_fit$correlation %||% NA_real_
))
