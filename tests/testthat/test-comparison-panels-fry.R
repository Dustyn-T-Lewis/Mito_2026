source(here::here("04_Figures", "functions", "F03-F04_comparison_panels.R"))

test_that("running_es rises on clustered hits and ends near zero", {
  # hits concentrated at the top of the ranking -> ES peaks positive early
  in_set <- c(rep(TRUE, 10), rep(FALSE, 90))
  stat <- sort(rnorm(100), decreasing = TRUE)
  es <- running_es(stat, in_set)
  expect_length(es, 100)
  expect_gt(max(es), 0.3) # meaningful enrichment peak
  expect_lt(abs(es[length(es)]), 1e-8) # walk returns to zero at the end
})
