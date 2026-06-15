#!/usr/bin/env Rscript
# =============================================================================
# 03_impute_imp4p.R  --  Mito Stage 02: mixture multiple-imputation (imp4p)
# imp4p::impute.mi(tab, conditions=group, methodMCAR="mle", methodMNAR="igcda")
# -- mixture MI for label-free proteomics. Feeds 03_DEP/b_imputed + figures.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(imp4p); library(dplyr); library(tibble); library(readr)
})
set.seed(42)
data_dir <- here("02_Normalization", "c_data")

dal <- readRDS(file.path(data_dir, "03_DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cond <- factor(dal$metadata$group[match(colnames(mat), dal$metadata$Col_ID)])
cat(sprintf("Loaded normalized matrix: %d x %d | %.1f%% missing | conditions: %s\n",
            nrow(mat), ncol(mat), mean(is.na(mat)) * 100, paste(levels(cond), collapse = "/")))

imp <- impute.mi(tab = mat, conditions = cond, methodMCAR = "mle", methodMNAR = "igcda", progress.bar = FALSE)
imp <- as.matrix(imp); dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "imp4p impute.mi", methodMCAR = "mle", methodMNAR = "igcda")
saveRDS(dal, file.path(data_dir, "05_DAList_imputed_imp4p.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(imp)), file.path(data_dir, "05_imputed_imp4p.csv"))
cat(sprintf("Done: imp4p mixture-MI imputed %d cells -> 05_DAList_imputed_imp4p.rds\n", sum(is.na(mat))))
