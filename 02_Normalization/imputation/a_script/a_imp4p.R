#!/usr/bin/env Rscript
# STANDALONE imputation option A: imp4p (exploratory arm; not the canonical DEP input).
#
# imp4p (Giai Gianetto et al.) is purpose-built for label-free proteomics with a MIXTURE
# of MCAR/MNAR missingness. It estimates the mechanism ITSELF (estim.bound -> estim.mix ->
# prob.mcar -> mi.mix); we supply only the data and the experimental-group factor, so no
# external MAR/MNAR classifier is used. Output feeds the 03_DEP/b_imputed comparison;
# the figures use the missForest arm, which is the more concordant imputer here.

pacman::p_load(proteoDA, here, imp4p)
source(here("02_Normalization", "imputation", "a_script", "_impute_helpers.R"))
set.seed(42) # impute.mi is stochastic (multiple imputation)

norm <- load_normalized_matrix()
mat <- norm$mat
cond <- factor(norm$dal$metadata$group[match(colnames(mat), norm$dal$metadata$Col_ID)])
cat(sprintf(
  "[imp4p] %d x %d | %.1f%% missing | conditions: %s\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100, paste(levels(cond), collapse = "/")
))

# mle for the MCAR part, igcda for the MNAR part; imp4p decides the per-protein mix.
imp <- as.matrix(impute.mi(
  tab = mat, conditions = cond,
  methodMCAR = "mle", methodMNAR = "igcda", progress.bar = FALSE
))

save_imputed_dalist(
  norm$dal, imp,
  list(method = "imp4p::impute.mi", methodMCAR = "mle", methodMNAR = "igcda"), "imp4p"
)
cat(sprintf("[imp4p] done: imputed %d cells -> DAList_imputed_imp4p.rds\n", sum(is.na(mat))))
