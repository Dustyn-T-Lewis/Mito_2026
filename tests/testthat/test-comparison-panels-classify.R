source(here::here("04_Figures", "functions", "F01-F03_style_palettes_theme.R"))
source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("classify_quadrant assigns signed quadrants and gates on significance", {
  comb <- tibble::tibble(
    uniprot_id = c("P1", "P2", "P3", "P4", "P5", "P6"),
    gene = c("Ga", "Gb", "Gc", "Gd", "Ge", "Gf"),
    logFC_CTLvPHE = c(1, -1, 1, -1, 0.9, 1),
    logFC_PHEvPHE_MITO = c(-1, 1, 1, -1, -0.8, -1),
    pi_score_CTLvPHE = c(0.01, 0.01, 0.01, 0.01, 0.20, NA),
    pi_score_PHEvPHE_MITO = c(0.01, 0.01, 0.01, 0.01, 0.20, NA)
  )
  q <- classify_quadrant(comb, "CTLvPHE", "PHEvPHE_MITO", pi_thresh = 0.05)
  expect_equal(q$quadrant[q$uniprot_id == "P1"], "x_up_y_dn") # reversed
  expect_equal(q$quadrant[q$uniprot_id == "P2"], "x_dn_y_up") # reversed
  expect_equal(q$quadrant[q$uniprot_id == "P3"], "x_up_y_up") # concordant
  expect_equal(q$quadrant[q$uniprot_id == "P4"], "x_dn_y_dn") # concordant
  expect_equal(q$quadrant[q$uniprot_id == "P5"], "ns") # neither significant
  expect_equal(q$quadrant[q$uniprot_id == "P6"], "ns") # NA pi_score treated as ns
})
