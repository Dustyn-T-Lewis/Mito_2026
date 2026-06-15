#!/usr/bin/env Rscript
# =============================================================================
# imputation/a_imp4p.R  --  STANDALONE imputation option A: imp4p
#
# imp4p (Giai Gianetto et al.) is purpose-built for label-free proteomics with a
# MIXTURE of MCAR/MNAR missingness. It estimates the mechanism ITSELF (internally:
# estim.bound -> estim.mix -> prob.mcar -> mi.mix); we supply only the data and
# the experimental-group factor. No external MAR/MNAR classifier is used.
#   imp4p::impute.mi(tab, conditions, methodMCAR = "mle", methodMNAR = "igcda")
# Output: 02_Normalization/c_data/DAList_imputed_imp4p.rds  (figures + 03_DEP/b_imputed)
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(imp4p); library(dplyr); library(tibble); library(readr)
})
set.seed(42)
data_dir <- here("02_Normalization", "c_data")

dal  <- readRDS(file.path(data_dir, "DAList_normalized.rds"))
mat  <- as.matrix(dal$data)
cond <- factor(dal$metadata$group[match(colnames(mat), dal$metadata$Col_ID)])
cat(sprintf("[imp4p] %d x %d | %.1f%% missing | conditions: %s\n",
            nrow(mat), ncol(mat), mean(is.na(mat)) * 100, paste(levels(cond), collapse = "/")))

imp <- as.matrix(impute.mi(tab = mat, conditions = cond,
                           methodMCAR = "mle", methodMNAR = "igcda", progress.bar = FALSE))
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "imp4p::impute.mi", methodMCAR = "mle", methodMNAR = "igcda")
saveRDS(dal, file.path(data_dir, "DAList_imputed_imp4p.rds"))
cat(sprintf("[imp4p] done: imputed %d cells -> DAList_imputed_imp4p.rds\n", sum(is.na(mat))))
