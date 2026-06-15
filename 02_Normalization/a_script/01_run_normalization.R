#!/usr/bin/env Rscript
# =============================================================================
# 01_run_normalization.R  --  Mito Stage 02: Normalization (+ imputation)
# cycloess (native normalize_data) + imputed matrix for figures two ways:
# Perseus (native, primary) + missForest (comparison). DEP uses NON-imputed.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble); library(ggplot2)
})
set.seed(42)
report_dir <- here("02_Normalization", "b_reports"); data_dir <- here("02_Normalization", "c_data")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE); dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("01_Filtering", "c_data", "01_DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

write_norm_report(dal, grouping_column = "group", output_dir = report_dir,
                  filename = "01_norm_comparison.pdf", overwrite = TRUE)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal, color_column = "group", output_dir = report_dir, filename = "02_qc_post.pdf", overwrite = TRUE)

write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)), file.path(data_dir, "02_normalized.csv"))
saveRDS(dal, file.path(data_dir, "03_DAList_normalized.rds"))
was_na <- is.na(dal$data)
cat(sprintf("Normalized: %d x %d | %.1f%% missing\n", nrow(dal$data), ncol(dal$data), mean(was_na) * 100))

# Imputation A: Perseus (native, primary figure matrix)
dal_perseus <- perseus_impute(dal, shift = 1.8, width = 0.3, seed = 42)
saveRDS(dal_perseus, file.path(data_dir, "04_DAList_imputed_perseus.rds"))

# Imputation B: missForest (comparison)
mat <- dal$data[order(rownames(dal$data)), ]; set.seed(42)
mf <- missForest::missForest(t(mat), maxiter = 10, ntree = 100, verbose = FALSE)
mat_mf <- t(mf$ximp); dimnames(mat_mf) <- dimnames(mat)
dal_mf <- dal; dal_mf$data <- mat_mf[rownames(dal$data), ]
saveRDS(dal_mf, file.path(data_dir, "05_DAList_imputed_missforest.rds"))
cat(sprintf("missForest OOB error: %.4f\n", as.numeric(mf$OOBerror[1])))

idx <- which(was_na)
cmp <- tibble(perseus = dal_perseus$data[idx], missforest = mat_mf[rownames(dal$data), ][idx])
rho <- cor(cmp$perseus, cmp$missforest, method = "spearman")
cat(sprintf("Imputed-cell agreement (Perseus vs missForest): Spearman rho = %.3f over %d cells\n", rho, nrow(cmp)))
ggsave(file.path(report_dir, "03_impute_perseus_vs_missforest.png"),
       ggplot(cmp, aes(perseus, missforest)) + geom_point(alpha = 0.15, size = 0.5) +
         geom_abline(slope = 1, colour = "red") +
         labs(title = sprintf("Imputed: Perseus vs missForest (rho=%.3f)", rho)) + theme_minimal(base_size = 11),
       width = 6, height = 6, dpi = 130)
write_csv(cmp, file.path(data_dir, "06_imputation_comparison.csv"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done -> %s/ (Perseus = primary figure matrix; missForest = comparison)\n", data_dir))
