#!/usr/bin/env Rscript
# STANDALONE imputation option B: MsCoreUtils hybrid (exploratory arm; not the canonical DEP input).
#
# Mechanism-aware hybrid, each tool used as designed:
#   1. imputeLCMD::model.Selector() classifies each protein MAR (1) vs MNAR (0) (Lazar 2016).
#   2. MsCoreUtils::impute_matrix(method = "mixed", ...) applies kNN to the MAR subset and
#      QRILC (left-censored) to the MNAR subset.

pacman::p_load(proteoDA, here, MsCoreUtils, imputeLCMD)
source(here("02_Normalization", "imputation", "a_script", "_impute_helpers.R"))
set.seed(42) # QRILC draws from a truncated Gaussian (stochastic)

norm <- load_normalized_matrix()
mat <- norm$mat
cat(sprintf("[mscoreutils] %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

# model.Selector gives the MAR mask; impute_matrix routes MAR->knn, MNAR->QRILC.
ms <- model.Selector(mat)
randna <- as.logical(ms[[1]]) # TRUE = MAR feature
cat(sprintf("[mscoreutils] model.Selector split: %d MAR / %d MNAR features\n", sum(randna), sum(!randna)))

imp <- impute_matrix(mat, method = "mixed", randna = randna, mar = "knn", mnar = "QRILC")

save_imputed_dalist(
  norm$dal, imp,
  list(
    method = "MsCoreUtils mixed (imputeLCMD model.Selector)",
    mar = "knn", mnar = "QRILC", n_mar = sum(randna), n_mnar = sum(!randna)
  ), "mscoreutils"
)
cat(sprintf("[mscoreutils] done: hybrid (knn/QRILC) imputed %d cells -> DAList_imputed_mscoreutils.rds\n", sum(is.na(mat))))
