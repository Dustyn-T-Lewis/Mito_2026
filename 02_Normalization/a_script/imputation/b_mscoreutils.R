#!/usr/bin/env Rscript
# =============================================================================
# imputation/b_mscoreutils.R  --  STANDALONE imputation option B: MsCoreUtils hybrid
#
# Mechanism-aware hybrid using each tool AS DESIGNED, independent of our code:
#   1. imputeLCMD::model.Selector() -- the designed MAR/MNAR classifier
#      (Lazar et al. 2016): per protein, 1 = MAR, 0 = MNAR.
#   2. MsCoreUtils::impute_matrix(method = "mixed", randna = <MAR>, mar = "knn",
#      mnar = "QRILC") -- applies kNN to the MAR subset, QRILC to the MNAR subset.
# Output: 02_Normalization/c_data/DAList_imputed_mscoreutils.rds
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(MsCoreUtils); library(imputeLCMD)
  library(dplyr); library(tibble); library(readr)
})
set.seed(42)
data_dir <- here("02_Normalization", "c_data")

dal <- readRDS(file.path(data_dir, "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf("[mscoreutils] %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

# Designed MAR/MNAR detection (imputeLCMD), NOT our consensus
ms     <- model.Selector(mat)
randna <- as.logical(ms[[1]])                    # TRUE = MAR feature
cat(sprintf("[mscoreutils] model.Selector split: %d MAR / %d MNAR features\n", sum(randna), sum(!randna)))

imp <- impute_matrix(mat, method = "mixed", randna = randna, mar = "knn", mnar = "QRILC")
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "MsCoreUtils mixed (imputeLCMD model.Selector)",
                       mar = "knn", mnar = "QRILC", n_mar = sum(randna), n_mnar = sum(!randna))
saveRDS(dal, file.path(data_dir, "DAList_imputed_mscoreutils.rds"))
cat(sprintf("[mscoreutils] done: hybrid (knn/QRILC) imputed %d cells -> DAList_imputed_mscoreutils.rds\n", sum(is.na(mat))))
