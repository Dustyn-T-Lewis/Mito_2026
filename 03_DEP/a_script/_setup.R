# Stage 03 setup — sourced by 01_run_dep.R, 02_generate_reports.R,
# 03_run_robustness.R. Defines the shared project config + stage paths.
# NOT run directly.

suppressPackageStartupMessages(library(here))

# Shared project design config.
source(here::here("00_input", "h9c2_design.R"))

# Stage paths
NORM_DAT <- here::here("01_normalization", "c_data")
IMP_DAT  <- here::here("02_Imputation", "c_data")
DAT      <- here::here("03_DEP", "c_data")
RPT      <- here::here("03_DEP", "b_reports")
PDA      <- file.path(RPT, "01_proteoDA")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(PDA, recursive = TRUE, showWarnings = FALSE)

# Inputs from Stages 01-02
NORM_CSV    <- file.path(NORM_DAT, "02_normalized.csv")        # non-imputed matrix
NORM_DALIST <- file.path(NORM_DAT, "03_DAList_normalized.rds")
IMP_DALIST  <- file.path(IMP_DAT, "01_DAList_imputed.rds")     # for sensitivity check

# Output workbook (built by 01_run_dep.R, extended by 02/03)
XLSX <- file.path(DAT, "03_DEP_results.xlsx")
