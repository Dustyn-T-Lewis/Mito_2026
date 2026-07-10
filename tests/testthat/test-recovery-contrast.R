test_that("non-imputed DE emits the Recovery contrast, well-formed", {
  csv <- here::here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv")
  skip_if_not(file.exists(csv), "run 01_run_dep.R first")
  res <- readr::read_csv(csv, show_col_types = FALSE)

  expect_true("CTLvPHE_MITO" %in% res$contrast)
  rec <- dplyr::filter(res, contrast == "CTLvPHE_MITO")
  expect_gt(nrow(rec), 4000)
  expect_true(all(c("logFC", "t", "P.Value", "adj.P.Val", "pi_score", "sig_pi") %in% names(rec)))
  expect_false(all(is.na(rec$logFC)))
  expect_true(all(rec$pi_score >= 0 & rec$pi_score <= 1, na.rm = TRUE))
})

test_that("imputed missForest arm also emits the Recovery contrast", {
  csv <- here::here("03_DEP", "b_imputed", "c_data", "missforest", "combined_results_pi.csv")
  skip_if_not(file.exists(csv), "run 01_run_dep_imputed.R first")
  res <- readr::read_csv(csv, show_col_types = FALSE)
  expect_true("CTLvPHE_MITO" %in% res$contrast)
})
