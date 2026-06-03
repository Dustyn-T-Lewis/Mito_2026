#!/usr/bin/env Rscript
# F01 orchestrator — QC + proteome-wide overview.
# Builds the MAIN composite (dual-PCA + DA counts + p/Π distributions + overlap)
# and the two SUPP QC pages (normalization, imputation). Reads existing pipeline
# outputs only; never re-runs 01-03.

source(here::here("05_Figures", "F01_QC_overview", "a_script", "01_main_panels.R"))
source(here::here("05_Figures", "F01_QC_overview", "a_script", "02_supp_panels.R"))

message("F01 complete")
