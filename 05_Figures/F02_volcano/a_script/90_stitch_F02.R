#!/usr/bin/env Rscript
# F02 orchestrator — asymmetric all-contrast volcano + enrichment-ring composite.
# Reads existing pipeline outputs + frozen rat fGSEA cache; never re-runs 01-03.

source(here::here("05_Figures", "F02_volcano", "a_script", "01_main_panels.R"))

message("F02 complete")
