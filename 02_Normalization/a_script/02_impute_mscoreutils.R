#!/usr/bin/env Rscript
# =============================================================================
# 02_impute_mscoreutils.R  --  Mito Stage 02: hybrid imputation (MsCoreUtils)
# MsCoreUtils::impute_matrix(method="mixed", randna=<MAR features>, mar="knn",
# mnar="QRILC") -- Lazar 2016 hybrid; randna from R/mar_mnar.R consensus.
# Feeds 03_DEP/b_imputed + figures. DEP/a_non_imputed uses non-imputed matrix.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(MsCoreUtils); library(dplyr); library(tibble); library(readr)
})
set.seed(42)
source(here("R", "mar_mnar.R"))
data_dir <- here("02_Normalization", "c_data")

dal <- readRDS(file.path(data_dir, "03_DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf("Loaded normalized matrix: %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

cls <- classify_mar_mnar(mat); randna <- cls$classification != "MNAR"
cat(sprintf("MAR/MNAR split: %d MAR(+Complete) / %d MNAR features\n", sum(randna), sum(!randna)))
imp <- impute_matrix(mat, method = "mixed", randna = randna, mar = "knn", mnar = "QRILC")
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "MsCoreUtils mixed", mar = "knn", mnar = "QRILC")
saveRDS(dal, file.path(data_dir, "04_DAList_imputed_mscoreutils.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(imp)), file.path(data_dir, "04_imputed_mscoreutils.csv"))
cat(sprintf("Done: hybrid (knn/QRILC) imputed %d cells -> 04_DAList_imputed_mscoreutils.rds\n", sum(is.na(mat))))
