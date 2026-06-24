#!/usr/bin/env Rscript
# STANDALONE imputation option A: imp4p (exploratory arm; not the canonical DEP input).
#
# imp4p (Giai Gianetto et al.) is purpose-built for label-free proteomics with a MIXTURE
# of MCAR/MNAR missingness. It estimates the mechanism ITSELF (estim.bound -> estim.mix ->
# prob.mcar -> mi.mix); we supply only the data and the experimental-group factor, so no
# external MAR/MNAR classifier is used. Output feeds figures + 03_DEP/b_imputed.

pacman::p_load(proteoDA, here, imp4p)
set.seed(42) # impute.mi is stochastic (multiple imputation)
norm_dir <- here("02_Normalization", "c_data") # read stage-02 normalized matrix
data_dir <- here("02_Normalization", "imputation", "c_data") # write imputed DAList here
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(file.path(norm_dir, "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cond <- factor(dal$metadata$group[match(colnames(mat), dal$metadata$Col_ID)])
cat(sprintf(
  "[imp4p] %d x %d | %.1f%% missing | conditions: %s\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100, paste(levels(cond), collapse = "/")
))

# mle for the MCAR part, igcda for the MNAR part; imp4p decides the per-protein mix.

imp <- as.matrix(impute.mi(
  tab = mat, conditions = cond,
  methodMCAR = "mle", methodMNAR = "igcda", progress.bar = FALSE
))
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "imp4p::impute.mi", methodMCAR = "mle", methodMNAR = "igcda")
saveRDS(dal, file.path(data_dir, "DAList_imputed_imp4p.rds"))
cat(sprintf("[imp4p] done: imputed %d cells -> DAList_imputed_imp4p.rds\n", sum(is.na(mat))))
