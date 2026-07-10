source(here::here("04_Figures", "functions", "F01-F03_style_palettes_theme.R"))
source(here::here("04_Figures", "functions", "F01-F03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("build_nes_scatter pivots the cache and computes reversal rho", {
  cache <- tibble::tibble(
    pathway = rep(c("HALLMARK_A", "HALLMARK_B", "HALLMARK_C"), 2),
    database = "Hallmark",
    contrast = rep(c("CTLvPHE", "PHEvPHE_MITO"), each = 3),
    NES = c(2, 1, -2, -2, -1, 2), # perfectly reversed
    padj = c(0.01, 0.2, 0.01, 0.01, 0.2, 0.01),
    size = 50
  )
  out <- build_nes_scatter(cache, "CTLvPHE", "PHEvPHE_MITO")
  expect_s3_class(out$plot, "ggplot")
  expect_named(out$table, c("pathway", "database", "NES_x", "NES_y", "padj_x", "padj_y"), ignore.order = TRUE)
  expect_equal(nrow(out$table), 3)
  expect_lt(attr(out$table, "spearman_rho"), -0.9) # reversal
})
