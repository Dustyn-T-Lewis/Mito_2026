#!/usr/bin/env Rscript
# =============================================================================
# 01_run_normalization.R  --  Mito Stage 02: Normalization (cycloess)
# Native normalize_data(cycloess) -> canonical NON-imputed DAList. b_reports =
# proteoDA QC PDFs. c_data = RDS handoff only (imputation/ adds imputed DALists).
# =============================================================================

suppressPackageStartupMessages({ library(proteoDA); library(here) })
set.seed(42)
report_dir <- here("02_Normalization", "b_reports"); data_dir <- here("02_Normalization", "c_data")
clear_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE) }
clear_dir(report_dir); clear_dir(data_dir)

dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

write_norm_report(dal, grouping_column = "group", output_dir = report_dir,
                  filename = "norm_comparison.pdf", overwrite = TRUE)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal, color_column = "group", output_dir = report_dir,
                filename = "qc_normalized.pdf", overwrite = TRUE)

saveRDS(dal, file.path(data_dir, "DAList_normalized.rds"))   # clean handoff to Stage 03 / imputation
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Normalized (cycloess): %d x %d | %.1f%% missing -> DAList_normalized.rds\n",
            nrow(dal$data), ncol(dal$data), mean(is.na(dal$data)) * 100))
