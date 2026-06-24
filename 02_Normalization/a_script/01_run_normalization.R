#!/usr/bin/env Rscript
# Mito Stage 02: cycloess normalization -> canonical NON-imputed DAList for Stage 03
# (limma handles per-protein NAs). The imputation/ scripts add exploratory imputed DALists.

pacman::p_load(proteoDA, here)
set.seed(42)
report_dir <- here("02_Normalization", "b_reports")
data_dir <- here("02_Normalization", "c_data")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(report_dir)
clear_dir(data_dir)

# Stage 01 handoff: filtered, un-normalized intensities.
dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# write_norm_report compares all candidate methods on the un-normalized data, so it
# must run BEFORE normalize_data; write_qc_report then inspects the chosen result.
write_norm_report(dal,
  grouping_column = "group", output_dir = report_dir,
  filename = "norm_comparison.pdf", overwrite = TRUE
)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal,
  color_column = "group", output_dir = report_dir,
  filename = "qc_normalized.pdf", overwrite = TRUE
)

# Non-imputed normalized DAList -> Stage 03 / imputation arms.
saveRDS(dal, file.path(data_dir, "DAList_normalized.rds"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf(
  "Normalized (cycloess): %d x %d | %.1f%% missing -> DAList_normalized.rds\n",
  nrow(dal$data), ncol(dal$data), mean(is.na(dal$data)) * 100
))
