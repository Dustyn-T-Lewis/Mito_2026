# Stage 02 setup — sourced by 01_impute.R and 02_generate_reports.R.
# Defines the shared project config + stage paths. NOT run directly.

suppressPackageStartupMessages(library(here))

# Shared project design config.
source(here::here("00_input", "h9c2_design.R"))

# Stage paths
NORM_DAT <- here::here("01_normalization", "c_data")
DAT      <- here::here("02_Imputation", "c_data")
RPT      <- here::here("02_Imputation", "b_reports")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# Inputs from Stage 01
NORM_CSV    <- file.path(NORM_DAT, "02_normalized.csv")        # numeric matrix
NORM_DALIST <- file.path(NORM_DAT, "03_DAList_normalized.rds") # proteoDA object

# Optional imputation-benchmark output. The core imputation script runs fine
# without it; if present it is appended to the supplement workbook as a sheet.
BENCH_RANKING <- file.path(DAT, "benchmark", "04_composite_ranking.csv")
