#!/usr/bin/env Rscript
# =============================================================================
# imputation/c_our_approach.R  --  STANDALONE imputation option C: our hybrid
#
# Our mechanism-aware hybrid, fully interpretable:
#   1. R/mar_mnar.R -- our 3-method consensus (k-means + global-logistic +
#      left-tail; >=2/3 votes = MNAR) classifies each protein.
#   2. MAR (+ Complete) proteins  -> missForest  (benchmarked top performer here;
#      borrows from correlated proteins).
#   3. MNAR proteins              -> imputeLCMD::impute.MinProb (down-shifted,
#      left-censored draw).
#   Each missing cell is filled by the method matching its protein's class.
# Output: 02_Normalization/c_data/DAList_imputed_our.rds
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(missForest); library(imputeLCMD)
  library(dplyr); library(tibble); library(readr)
})
set.seed(42)
source(here("R", "mar_mnar.R"))
data_dir <- here("02_Normalization", "c_data")

dal <- readRDS(file.path(data_dir, "03_DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf("[our] %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

cls <- classify_mar_mnar(mat)
mnar_row <- cls$classification == "MNAR"
cat(sprintf("[our] consensus split: %d MAR/Complete / %d MNAR features\n", sum(!mnar_row), sum(mnar_row)))

set.seed(42)
mf       <- missForest(t(mat), maxiter = 10, ntree = 100, verbose = FALSE)
mat_mar  <- t(mf$ximp); dimnames(mat_mar) <- dimnames(mat)
mat_mnar <- as.matrix(impute.MinProb(mat, q = 0.01, tune.sigma = 1)); dimnames(mat_mnar) <- dimnames(mat)

imp <- mat
na  <- is.na(mat)
mnar_cell <- na & matrix(mnar_row, nrow(mat), ncol(mat))
mar_cell  <- na & !matrix(mnar_row, nrow(mat), ncol(mat))
imp[mnar_cell] <- mat_mnar[mnar_cell]
imp[mar_cell]  <- mat_mar[mar_cell]
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "our hybrid (consensus -> missForest[MAR] + MinProb[MNAR])",
                       n_mar = sum(!mnar_row), n_mnar = sum(mnar_row), oob = as.numeric(mf$OOBerror[1]))
saveRDS(dal, file.path(data_dir, "DAList_imputed_our.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(imp)), file.path(data_dir, "imputed_our.csv"))
cat(sprintf("[our] done: missForest(MAR)+MinProb(MNAR) imputed %d cells | OOB=%.4f -> DAList_imputed_our.rds\n",
            sum(na), as.numeric(mf$OOBerror[1])))
