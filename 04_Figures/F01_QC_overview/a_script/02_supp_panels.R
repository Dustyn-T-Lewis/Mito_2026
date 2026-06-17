#!/usr/bin/env Rscript
# F01 SUPP — DROPPED in the proteoDA-native rebuild.
# Every panel (A-N) read report intermediates that no longer exist:
#   P05$norm_interm / P05$imp_interm (normalization + imputation report RDS) and
#   the results.xlsx "DA_summary" sheet. None are produced by the current 01-03
#   pipeline, and the imputation benchmark CSV is absent, so these QC pages
#   (SUPP_F01_normalization A-G, SUPP_F01_imputation H-N) cannot be reproduced.
# The MAIN composite (01_main_panels.R) is unaffected.

message("F01 SUPP QC pages skipped: report-intermediate inputs (norm/imp RDS, DA_summary sheet, imputation benchmark) do not exist in the rebuilt pipeline.")
