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

test_that("build_rrho2 attributes the reversal peak to its true quadrant under a skewed sign balance", {
  withr::local_seed(42)
  n <- 400
  n_up_x <- 320 # disease t mostly positive (80/20 split, far from the 50/50 midpoint)
  tx <- c(rnorm(n_up_x, mean = 3), rnorm(n - n_up_x, mean = -3))
  ty <- -tx + rnorm(n, sd = 0.2) # near-perfect reversal, so rescue t skews mostly negative
  shuffle <- sample(n)
  comb <- tibble::tibble(
    uniprot_id = sprintf("P%03d", seq_len(n)), gene = sprintf("G%03d", seq_len(n)),
    t_CTLvPHE = tx[shuffle],
    t_PHEvPHE_MITO = ty[shuffle]
  )
  out <- build_rrho2(comb, "CTLvPHE", "PHEvPHE_MITO")
  peak_quadrant <- out$table$quadrant[which.max(out$table$max_neglog10p)]
  # true up/down boundary sits at 80% of the ranks, not the naive 50% index
  # midpoint; the real concentration of overlap belongs to x-up/y-down ("ud")
  expect_identical(peak_quadrant, "ud")
})
