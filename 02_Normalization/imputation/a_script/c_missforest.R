#!/usr/bin/env Rscript
# STANDALONE imputation option C: missForest — the imputed arm the figures use (the DEP input stays non-imputed).
#
# missForest (Stekhoven & Buhlmann 2012) is a non-parametric random-forest imputer. It is a
# MAR/MCAR method: each protein's missing values are predicted from the multivariate structure
# of the OTHER samples, with no left-censoring model. The Stage-02 missingness here is largely
# left-censored MNAR (low-abundance dropout: Spearman[abundance, missingness] ~= -0.54), which
# missForest does not model. Even so it empirically tracks the non-imputed effect sizes most
# closely (Spearman 0.98-0.99 in 03_DEP) and barely adds DEPs, whereas imp4p and the MsCoreUtils
# knn/QRILC hybrid inflate DEP counts from missing-only proteins; so missForest is the arm the
# figures use, with the other two kept as comparisons.
#
# Orientation: the normalized matrix is proteins (rows) x samples (columns). missForest treats
# COLUMNS as the variables to impute, so run on the matrix as-is — each sample column's missing
# values are learned from the other 23 sample columns across all proteins (4806 training rows,
# 23 predictors: well-powered, exploits replicate correlation). Single-threaded for a
# reproducible seed.

pacman::p_load(proteoDA, here, missForest)
set.seed(42) # missForest is stochastic (random forests)
norm_dir <- here("02_Normalization", "c_data") # read stage-02 normalized matrix
data_dir <- here("02_Normalization", "imputation", "c_data") # write imputed DAList here
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(file.path(norm_dir, "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf("[missforest] %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

mf <- missForest(mat, maxiter = 10, ntree = 100, verbose = FALSE)
imp <- mf$ximp
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(
  method = "missForest::missForest", maxiter = 10, ntree = 100,
  OOB_NRMSE = unname(mf$OOBerror["NRMSE"])
)
saveRDS(dal, file.path(data_dir, "DAList_imputed_missforest.rds"))
cat(sprintf(
  "[missforest] done: imputed %d cells (OOB NRMSE = %.4f) -> DAList_imputed_missforest.rds\n",
  sum(is.na(mat)), unname(mf$OOBerror["NRMSE"])
))
