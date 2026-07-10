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

build_rrho2 <- function(comb, ctr_x, ctr_y, labels = c(ctr_x, ctr_y)) {
  col_x <- paste0("t_", ctr_x)
  col_y <- paste0("t_", ctr_y)
  d <- comb[!is.na(comb[[col_x]]) & !is.na(comb[[col_y]]), ]
  l1 <- data.frame(Genes = d$uniprot_id, score = d[[col_x]])
  l2 <- data.frame(Genes = d$uniprot_id, score = d[[col_y]])

  obj <- RRHO2::RRHO2_initialize(l1, l2, labels = labels, log10.ind = TRUE)
  hypermat <- obj$hypermat
  n_row <- nrow(hypermat)
  n_col <- ncol(hypermat)
  mid_row <- n_row %/% 2
  mid_col <- n_col %/% 2

  corner_max <- function(rows, cols) max(hypermat[rows, cols], na.rm = TRUE)
  top_rows <- seq_len(mid_row)
  bot_rows <- (mid_row + 1):n_row
  top_cols <- seq_len(mid_col)
  bot_cols <- (mid_col + 1):n_col

  ord1 <- l1$Genes[order(l1$score, decreasing = TRUE)]
  ord2 <- l2$Genes[order(l2$score, decreasing = TRUE)]
  n_shared_genes <- nrow(l1)
  mid_genes <- n_shared_genes %/% 2
  up1 <- ord1[seq_len(mid_genes)]
  dn1 <- ord1[(mid_genes + 1):n_shared_genes]
  up2 <- ord2[seq_len(mid_genes)]
  dn2 <- ord2[(mid_genes + 1):n_shared_genes]

  table <- tibble::tibble(
    quadrant = c("uu", "ud", "du", "dd"),
    max_neglog10p = c(
      corner_max(top_rows, top_cols), corner_max(top_rows, bot_cols),
      corner_max(bot_rows, top_cols), corner_max(bot_rows, bot_cols)
    ),
    n_shared = c(
      length(intersect(up1, up2)), length(intersect(up1, dn2)),
      length(intersect(dn1, up2)), length(intersect(dn1, dn2))
    )
  )

  hyper_long <- tibble::tibble(
    row = rep(seq_len(n_row), times = n_col),
    col = rep(seq_len(n_col), each = n_row),
    value = as.vector(hypermat)
  )

  plot <- ggplot2::ggplot(hyper_long, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = c("#2166AC", "white", "#B2182B"), na.value = "grey90",
      name = expression(-log[10](italic(p)))
    ) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::annotate(
      "text",
      x = c(n_col * 0.05, n_col * 0.95, n_col * 0.05, n_col * 0.95),
      y = c(n_row * 0.05, n_row * 0.05, n_row * 0.95, n_row * 0.95),
      label = c("Exacerbated", "Reversed", "Reversed", "Exacerbated"),
      fontface = "bold", size = 2, color = "grey20"
    ) +
    ggplot2::labs(
      x = bquote(.(labels[2]) ~ "rank" ~ (Up %->% Down)),
      y = bquote(.(labels[1]) ~ "rank" ~ (Up %->% Down)),
      subtitle = sprintf(
        paste0(
          "%d shared proteins · warm off-diagonal = transplant reverses ",
          "disease · No MTC (Cahill 2018)"
        ),
        n_shared_genes
      )
    ) +
    ggplot2::coord_equal() +
    FIG_THEME +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  list(plot = plot, table = table)
}
