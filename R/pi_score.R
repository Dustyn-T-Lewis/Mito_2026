# =============================================================================
# pi_score.R  --  Pi-score significance (Xiao 2014, Bioinformatics 30:801)
#
# Pi = P.Value ^ |logFC|  (Eq.2 form) -- ranks by significance AND effect size.
# sig_pi: +1 up / -1 down / 0 not-significant at pi_thresh.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(purrr) })

# results_list: named list of per-contrast data.frames (proteoDA dal$results),
# each with rownames = uniprot_id and columns logFC, P.Value, adj.P.Val, ...
compute_pi_scores <- function(results_list, pi_thresh = 0.05) {
  imap(results_list, function(res, cname) {
    as_tibble(res, rownames = "uniprot_id") |>
      mutate(
        pi_score = P.Value ^ abs(logFC),
        sig_pi = case_when(
          pi_score < pi_thresh & logFC > 0 ~ 1L,
          pi_score < pi_thresh & logFC < 0 ~ -1L,
          TRUE                             ~ 0L
        ),
        contrast = cname
      )
  })
}
