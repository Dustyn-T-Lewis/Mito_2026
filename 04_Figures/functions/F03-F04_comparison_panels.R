# Native contrast-pair comparison panels for F03 (concordance) and F04 (reversal):
# quadrant ORA, fry barcode, NES scatter, RRHO2. Each builder takes
# a wide DE table and an old-name contrast pair and returns list(plot, table).

classify_quadrant <- function(comb, ctr_x, ctr_y, pi_thresh = H9C2_PI_THRESH) {
  lfc_x <- comb[[paste0("logFC_", ctr_x)]]
  lfc_y <- comb[[paste0("logFC_", ctr_y)]]
  sig_x <- comb[[paste0("pi_score_", ctr_x)]] < pi_thresh
  sig_y <- comb[[paste0("pi_score_", ctr_y)]] < pi_thresh
  sig_x[is.na(sig_x)] <- FALSE
  sig_y[is.na(sig_y)] <- FALSE
  quadrant <- dplyr::case_when(
    !(sig_x | sig_y) ~ "ns",
    lfc_x >= 0 & lfc_y >= 0 ~ "x_up_y_up",
    lfc_x >= 0 & lfc_y < 0 ~ "x_up_y_dn",
    lfc_x < 0 & lfc_y >= 0 ~ "x_dn_y_up",
    TRUE ~ "x_dn_y_dn"
  )
  tibble::tibble(
    uniprot_id = comb$uniprot_id, gene = comb$gene,
    lfc_x = lfc_x, lfc_y = lfc_y, sig_x = sig_x, sig_y = sig_y, quadrant = quadrant
  )
}
