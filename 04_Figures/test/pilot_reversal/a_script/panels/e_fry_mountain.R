# Fry mountain (ported from YvO panel_C_fry), generalised over contrast pairs.
# Proteins are ranked by one contrast's t-statistic; a signature set defined from
# another contrast is drawn as a GSEA-style running enrichment score with a member
# barcode and the fry direction/p. Leading-edge proteins feed the flanking ORA.

running_es <- function(t_vals, in_set) {
  n <- length(t_vals)
  n_hit <- sum(in_set)
  if (n_hit == 0) {
    return(rep(0, n))
  }
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -1 / (n - n_hit)))
}

sig_stars <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

.pretty_pw <- function(x) {
  x <- sub("^(HALLMARK|REACTOME|KEGG|GOSLIM|MITOCARTA|GOBP|GOCC|GOMF)_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x)))
}

.es_barcode <- function(t_rank, in_col, es_col, fry_row, title, color, n_all) {
  marks <- t_rank[t_rank[[in_col]], ]
  line_color <- if (isTRUE(fry_row$fdr < 0.05)) color else scales::alpha(color, 0.4)
  lab <- sprintf("fry %s, %s (n = %d)", fry_row$dir, fmt_p(fry_row$fdr), fry_row$n)
  es <- ggplot2::ggplot(t_rank, ggplot2::aes(rank, .data[[es_col]])) +
    ggplot2::geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    ggplot2::geom_line(color = line_color, linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = lab, hjust = 1.05, vjust = 1.5, size = 1.9, fontface = "bold", color = "grey20") +
    ggplot2::labs(title = title, x = NULL, y = "ES") +
    ggplot2::scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 0, 0, "mm")
    )
  bc <- ggplot2::ggplot(marks, ggplot2::aes(rank, xend = rank, y = 0, yend = 1)) +
    ggplot2::geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    ggplot2::scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "grey97"),
      plot.margin = ggplot2::margin(0, 1, 0, 2, "mm")
    )
  list(es = es, bc = bc)
}

# top 5 pathways always shown; non-significant bars shaded lighter
.ora_bars <- function(ora_df, title, color) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "no pathways", size = 2, color = "grey60"))
  }
  bars <- ora_df |>
    dplyr::slice_head(n = 5) |>
    dplyr::mutate(
      nlp = -log10(pmax(padj, 1e-20)), star = sig_stars(padj), sig = padj < 0.05,
      fill = ifelse(sig, scales::alpha(color, 0.9), scales::alpha(color, 0.35)),
      y = rev(dplyr::row_number()), label = .pretty_pw(pathway)
    )
  xmax <- max(bars$nlp, na.rm = TRUE)
  ggplot2::ggplot(bars, ggplot2::aes(y = y)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = nlp, ymin = y - 0.42, ymax = y + 0.42),
      fill = bars$fill, color = "black", linewidth = 0.25
    ) +
    ggplot2::geom_text(ggplot2::aes(x = xmax * 0.03, y = y, label = label),
      hjust = 0, size = 1.7, fontface = "bold", color = ifelse(bars$sig, "white", "grey15"), lineheight = 0.85
    ) +
    ggplot2::geom_text(ggplot2::aes(x = nlp + xmax * 0.02, label = star),
      hjust = 0, size = 2.2, fontface = "bold", color = "black"
    ) +
    ggplot2::labs(title = title, x = expression(-log[10](p[adj])), y = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, xmax * 1.18), breaks = scales::pretty_breaks(3), expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, nrow(bars) + 0.6), expand = c(0, 0)) +
    FIG_THEME +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

build_mountain <- function(t_rank, fry_up, fry_dn, ora_up, ora_dn, n_all, cfg) {
  up <- .es_barcode(t_rank, "in_up", "es_up", fry_up, cfg$up_title, unname(DIR_COLORS["Up"]), n_all)
  dn <- .es_barcode(t_rank, "in_down", "es_down", fry_dn, cfg$dn_title, unname(DIR_COLORS["Down"]), n_all)
  t_area <- ggplot2::ggplot(t_rank, ggplot2::aes(rank, t_rank_stat)) +
    ggplot2::geom_area(fill = scales::alpha("grey50", 0.25), color = "grey45", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::labs(x = cfg$rank_label, y = NULL) +
    ggplot2::scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 1, 0, "mm"))

  design <- c(
    patchwork::area(1, 1, 1, 1), patchwork::area(2, 1, 2, 1),
    patchwork::area(3, 1, 3, 1), patchwork::area(4, 1, 4, 1),
    patchwork::area(5, 1, 5, 1), patchwork::area(1, 2, 2, 2),
    patchwork::area(3, 2, 4, 2)
  )
  up$es + up$bc + dn$es + dn$bc + t_area +
    .ora_bars(ora_up, cfg$ora_up_title, unname(DIR_COLORS["Up"])) +
    .ora_bars(ora_dn, cfg$ora_dn_title, unname(DIR_COLORS["Down"])) +
    patchwork::plot_layout(design = design, heights = c(2, 0.25, 2, 0.25, 0.7), widths = c(1.7, 1.6)) +
    patchwork::plot_annotation(
      title = cfg$title, subtitle = cfg$subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = FIG_TITLE_SIZE, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
      )
    )
}
