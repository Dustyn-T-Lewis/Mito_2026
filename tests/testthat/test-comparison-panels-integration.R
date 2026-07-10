source(here::here("04_Figures", "functions", "F01-F03_style_palettes_theme.R"))
source(here::here("04_Figures", "functions", "F01-F03_data_paths_and_loaders.R"))
source(here::here("04_Figures", "functions", "F01-F03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("all four panels build for the Disease x Rescue pair", {
  skip_if_not(file.exists(P05$comb), "run the DE stage first")
  skip_if_not(file.exists(here::here("04_Figures", "shared", "c_data", "fgsea_tstat_all_h9c2.csv")))
  comb <- load_combined_wide()
  cache <- load_fgsea_cache()
  expect_s3_class(build_nes_scatter(cache, "CTLvPHE", "PHEvPHE_MITO")$plot, "ggplot")
  expect_s3_class(build_rrho2(comb, "CTLvPHE", "PHEvPHE_MITO")$plot, "ggplot")
  expect_s3_class(build_quadrant_ora(comb, "CTLvPHE", "PHEvPHE_MITO")$plot, "gg")
})
