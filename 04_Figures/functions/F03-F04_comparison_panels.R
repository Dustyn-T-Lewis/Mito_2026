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

# RRHO2::RRHO2_initialize splits hypermat on the sign-based boundary
# (count of positive scores at each step index) plus an asymmetric NA strip
# (lenStrip = round(len * boundary)), not on the index midpoint. Re-derive
# the same boundary here so quadrant corners line up with RRHO2's own.
rrho2_quadrant_table <- function(obj, l1, l2) {
  hypermat <- obj$hypermat
  n1 <- nrow(l1)
  n2 <- nrow(l2)
  stepsize <- ceiling(min(sqrt(c(n1, n2))))
  score1 <- sort(l1$score, decreasing = TRUE)
  score2 <- sort(l2$score, decreasing = TRUE)
  step1 <- seq(1, n1, stepsize)
  step2 <- seq(1, n2, stepsize)
  len1 <- length(step1)
  len2 <- length(step2)
  boundary1 <- sum(score1[step1] > 0)
  boundary2 <- sum(score2[step2] > 0)
  strip1 <- round(len1 * 0.1)
  strip2 <- round(len2 * 0.1)

  safe_seq <- function(from, to) if (from > to) integer(0) else seq(from, to)
  corner_max <- function(rows, cols) {
    if (length(rows) == 0 || length(cols) == 0) {
      return(NA_real_)
    }
    block <- hypermat[rows, cols]
    if (all(is.na(block))) {
      return(NA_real_) # all-NA corner (empty quadrant): no cell to report
    }
    max(block, na.rm = TRUE)
  }

  up_rows <- safe_seq(1, boundary1)
  up_cols <- safe_seq(1, boundary2)
  dn_rows <- strip1 + safe_seq(boundary1 + 1, len1)
  dn_cols <- strip2 + safe_seq(boundary2 + 1, len2)

  tibble::tibble(
    quadrant = c("uu", "ud", "du", "dd"),
    max_neglog10p = c(
      corner_max(up_rows, up_cols), corner_max(up_rows, dn_cols),
      corner_max(dn_rows, up_cols), corner_max(dn_rows, dn_cols)
    ),
    n_shared = c(
      length(obj$genelist_uu$gene_list_overlap_uu),
      length(obj$genelist_ud$gene_list_overlap_ud),
      length(obj$genelist_du$gene_list_overlap_du),
      length(obj$genelist_dd$gene_list_overlap_dd)
    )
  )
}

rrho2_render_plot <- function(hypermat, labels, n_total) {
  n_row <- nrow(hypermat)
  n_col <- ncol(hypermat)
  hyper_long <- tibble::tibble(
    row = rep(seq_len(n_row), times = n_col),
    col = rep(seq_len(n_col), each = n_row),
    value = as.vector(hypermat)
  )

  ggplot2::ggplot(hyper_long, ggplot2::aes(col, row, fill = value)) +
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
        n_total
      )
    ) +
    ggplot2::coord_equal() +
    FIG_THEME +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

build_rrho2 <- function(comb, ctr_x, ctr_y, labels = c(ctr_x, ctr_y)) {
  col_x <- paste0("t_", ctr_x)
  col_y <- paste0("t_", ctr_y)
  d <- comb[!is.na(comb[[col_x]]) & !is.na(comb[[col_y]]), ]
  l1 <- data.frame(Genes = d$uniprot_id, score = d[[col_x]])
  l2 <- data.frame(Genes = d$uniprot_id, score = d[[col_y]])

  obj <- RRHO2::RRHO2_initialize(l1, l2, labels = labels, log10.ind = TRUE)
  quad_table <- rrho2_quadrant_table(obj, l1, l2)
  plot <- rrho2_render_plot(obj$hypermat, labels, nrow(l1))

  list(plot = plot, table = quad_table)
}

# Quadrant background tints for the NES scatter: reversal (ref_slope = -1)
# warms the off-diagonal (blue, disease reversed by transplant) and the
# diagonal (red, exacerbated); concordance (ref_slope = +1) warms only the
# diagonal (green, concordant).
nes_quadrant_tints <- function(ref_slope) {
  diag_fill <- if (ref_slope < 0) "#D6604D" else "#4DAF4A"
  quadrants <- tibble::tibble(
    xmin = c(0, -Inf, -Inf, 0),
    xmax = c(Inf, 0, 0, Inf),
    ymin = c(0, 0, -Inf, -Inf),
    ymax = c(Inf, Inf, 0, 0),
    fill = c(diag_fill, "#4393C3", diag_fill, "#4393C3")
  )
  if (ref_slope < 0) {
    return(quadrants)
  }
  quadrants[c(1, 3), ]
}

build_nes_scatter <- function(fgsea_cache, ctr_x, ctr_y, ref_slope = -1) {
  wide <- fgsea_cache |>
    dplyr::filter(contrast %in% c(ctr_x, ctr_y), database %in% CANONICAL_DBS) |>
    dplyr::select(pathway, database, contrast, NES, padj) |>
    tidyr::pivot_wider(names_from = contrast, values_from = c(NES, padj)) |>
    dplyr::rename(
      NES_x = !!paste0("NES_", ctr_x), NES_y = !!paste0("NES_", ctr_y),
      padj_x = !!paste0("padj_", ctr_x), padj_y = !!paste0("padj_", ctr_y)
    ) |>
    tidyr::drop_na(NES_x, NES_y)

  attr(wide, "spearman_rho") <- stats::cor(
    wide$NES_x, wide$NES_y,
    method = "spearman"
  )

  min_padj <- pmin(wide$padj_x, wide$padj_y, na.rm = TRUE)
  sig <- !is.na(min_padj) & min_padj < 0.05
  off_diagonal <- (wide$NES_x * wide$NES_y) < 0
  pct_reversed <- if (any(sig)) 100 * mean(off_diagonal[sig]) else NA_real_

  callouts <- wide |>
    dplyr::mutate(
      diff = abs(NES_x - NES_y),
      label = clean_display_label(pathway)
    ) |>
    dplyr::slice_max(diff, n = 8, with_ties = FALSE)

  plot <- ggplot2::ggplot(wide, ggplot2::aes(NES_x, NES_y)) +
    ggplot2::geom_rect(
      data = nes_quadrant_tints(ref_slope),
      ggplot2::aes(
        xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
      ),
      inherit.aes = FALSE, alpha = 0.15
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_abline(
      slope = ref_slope, intercept = 0, linetype = "dashed", color = "grey40"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = database, size = -log10(pmin(padj_x, padj_y))),
      alpha = 0.7
    ) +
    ggplot2::scale_color_manual(values = ORA_DB_COLORS) +
    (if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(
        data = callouts, ggplot2::aes(label = label),
        size = 2, max.overlaps = Inf, segment.size = 0.2
      )
    } else {
      ggplot2::geom_text(
        data = callouts, ggplot2::aes(label = label),
        size = 2
      )
    }) +
    ggplot2::labs(
      x = sprintf("%s NES", ctr_x), y = sprintf("%s NES", ctr_y),
      size = expression(-log[10](italic(padj))),
      subtitle = sprintf(
        "ρ = %.2f · %.0f%% reversed",
        attr(wide, "spearman_rho"), pct_reversed
      )
    ) +
    FIG_THEME

  list(plot = plot, table = wide)
}
