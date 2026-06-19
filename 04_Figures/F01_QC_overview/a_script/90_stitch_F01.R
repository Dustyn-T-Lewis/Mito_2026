#!/usr/bin/env Rscript
# F01 orchestrator — QC + proteome-wide overview.
# Builds the MAIN composite (dual-PCA + DA counts + p/Π distributions + overlap).
# Reads existing pipeline outputs only; never re-runs 01-03.

source(here::here("04_Figures", "F01_QC_overview", "a_script", "01_main_panels.R"))

message("F01 complete")
