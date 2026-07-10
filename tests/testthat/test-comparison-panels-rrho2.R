source(here::here("04_Figures", "functions", "F01-F03_style_palettes_theme.R"))
source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("build_rrho2 detects anti-correlated ranks as off-diagonal reversal", {
  withr::local_seed(42)
  n <- 400
  tx <- rnorm(n)
  comb <- tibble::tibble(
    uniprot_id = sprintf("P%03d", seq_len(n)), gene = sprintf("G%03d", seq_len(n)),
    t_CTLvPHE = tx,
    t_PHEvPHE_MITO = -tx + rnorm(n, sd = 0.2) # near-perfect reversal
  )
  out <- build_rrho2(comb, "CTLvPHE", "PHEvPHE_MITO")
  expect_s3_class(out$plot, "ggplot")
  expect_true(all(c("uu", "ud", "du", "dd") %in% out$table$quadrant))
  off <- max(out$table$max_neglog10p[out$table$quadrant %in% c("ud", "du")])
  diag <- max(out$table$max_neglog10p[out$table$quadrant %in% c("uu", "dd")])
  expect_gt(off, diag) # reversal dominates the concordant corners
})
