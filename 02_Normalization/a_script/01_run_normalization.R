#!/usr/bin/env Rscript
# =============================================================================
# 01_run_normalization.R  --  Mito Stage 02: Normalization (cycloess)
# Native normalize_data(cycloess) -> canonical NON-imputed matrix for
# 03_DEP/a_non_imputed. Imputation in 02_impute_mscoreutils.R + 03_impute_imp4p.R.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble)
})
set.seed(42)
report_dir <- here("02_Normalization", "b_reports"); data_dir <- here("02_Normalization", "c_data")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE); dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("01_Filtering", "c_data", "01_DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

write_norm_report(dal, grouping_column = "group", output_dir = report_dir,
                  filename = "01_norm_comparison.pdf", overwrite = TRUE)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal, color_column = "group", output_dir = report_dir,
                filename = "02_qc_post.pdf", overwrite = TRUE)

write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)), file.path(data_dir, "02_normalized.csv"))
saveRDS(dal, file.path(data_dir, "03_DAList_normalized.rds"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Normalized (cycloess): %d x %d | %.1f%% missing -> %s/\n",
            nrow(dal$data), ncol(dal$data), mean(is.na(dal$data)) * 100, data_dir))
