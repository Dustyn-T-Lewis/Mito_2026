source(here::here("04_Figures", "functions", "F01-F03_style_palettes_theme.R"))
source(here::here("04_Figures", "functions", "F01-F03_data_paths_and_loaders.R"))
source(here::here("04_Figures", "functions", "F01-F03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("build_quadrant_ora enriches a quadrant seeded with a pathway's genes", {
  set_a <- paste0("Gene", 1:20) # a fake pathway
  bg <- paste0("Bg", 1:300)
  comb <- tibble::tibble(
    uniprot_id = c(set_a, bg), gene = c(set_a, bg),
    logFC_CTLvPHE = c(rep(1, 20), rep(0, 300)),
    logFC_PHEvPHE_MITO = c(rep(-1, 20), rep(0, 300)),
    pi_score_CTLvPHE = c(rep(0.001, 20), rep(0.9, 300)),
    pi_score_PHEvPHE_MITO = c(rep(0.001, 20), rep(0.9, 300))
  )
  pw <- list(FAKE_SET = set_a)
  out <- build_quadrant_ora(
    comb, "CTLvPHE", "PHEvPHE_MITO",
    pw_collection = pw, pw_by_db = list(Hallmark = "FAKE_SET"),
    universe = comb$gene
  )
  hit <- dplyr::filter(out$table, quadrant == "x_up_y_dn", pathway == "FAKE_SET")
  expect_gt(nrow(hit), 0)
  expect_lt(min(hit$padj), 0.05)
})
