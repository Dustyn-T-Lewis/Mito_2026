test_that("fgsea cache carries the Recovery contrast", {
  csv <- here::here("04_Figures", "shared", "c_data", "fgsea_tstat_all_h9c2.csv")
  skip_if_not(file.exists(csv), "run 02_build_fgsea_cache.R first")
  fg <- readr::read_csv(csv, show_col_types = FALSE)
  expect_true("CTLvPHE_MITO" %in% fg$contrast)
})

test_that("primary-contrast Pi counts are unchanged by adding Recovery", {
  skip_if_not(file.exists("/tmp/primary_summary_pre.csv"))
  post <- readr::read_csv(
    here::here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(contrast %in% c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction", "MITOvPHE_MITO")) |>
    dplyr::group_by(contrast) |>
    dplyr::summarise(n = dplyr::n(), pi = sum(sig_pi != 0), .groups = "drop") |>
    dplyr::mutate(n = as.double(n), pi = as.double(pi))
  pre <- readr::read_csv("/tmp/primary_summary_pre.csv", show_col_types = FALSE) |>
    tibble::as_tibble()
  expect_equal(post, pre)
})
