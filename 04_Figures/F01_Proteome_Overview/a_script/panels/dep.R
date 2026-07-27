# DEP count bars for F01 panel B, plus the shared dep_data() loader that both
# this panel and the effect-size panel (dep_effect.R) build from.

# Loads the per-contrast DE tables once and derives the threshold-count table,
# so the composite hands one object to both builders instead of leaking globals.
# Narrative order: Disease -> Transplant -> Rescue -> Interaction.
dep_data <- function() {
  core <- H9C2_CONTRAST_ORDER
  ctr_lab <- stats::setNames(contrast_brief(core), core)
  ctr_levels <- unname(ctr_lab)
  thr_levels <- c(
    "p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), paste0("Π < ", H9C2_PI_THRESH)
  )

  comb_long <- load_dep_long()
  results <- stats::setNames(
    lapply(core, \(c) as.data.frame(dplyr::filter(comb_long, contrast == c))),
    core
  )
  # Denominator is proteins, matching the per-protein rows the counts below
  # tally; gene symbols would be the wrong unit (~30 proteins fall back to an
  # accession for a name).
  n_total <- dplyr::n_distinct(comb_long$uniprot_id)

  counts <- dplyr::bind_rows(lapply(core, function(ctr) {
    r <- results[[ctr]]
    n_p <- sum(!is.na(r$P.Value) & r$P.Value < 0.05, na.rm = TRUE)
    n_fdr <- sum(!is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE)
    n_pi <- sum(!is.na(r$pi_score) & r$pi_score < H9C2_PI_THRESH, na.rm = TRUE)
    tibble::tibble(
      contrast = ctr_lab[ctr], threshold = thr_levels, n = c(n_p, n_fdr, n_pi)
    )
  })) |>
    dplyr::mutate(
      contrast = factor(contrast, levels = rev(ctr_levels)), # rev for coord_flip
      threshold = factor(threshold, levels = thr_levels),
      pct = 100 * n / n_total,
      fill_key = paste(as.character(contrast), threshold, sep = "___")
    ) |>
    dplyr::filter(n > 0)

  list(
    results = results, counts = counts, n_total = n_total,
    core = core, ctr_lab = ctr_lab, ctr_levels = ctr_levels,
    thr_levels = thr_levels
  )
}

# Panel B: DEP counts (horizontal nested-threshold)
build_dep_count_panel <- function(dep) {
  frac_df <- dep$counts
  CORE <- dep$core
  CTR_LAB <- dep$ctr_lab
  ctr_levels <- dep$ctr_levels
  THR_LEVELS <- dep$thr_levels

  # Contrast color map keyed on the rev-levelled factor labels
  SET_COLS <- stats::setNames(unname(CONTRAST_COLORS[CORE]), unname(CTR_LAB))
  # Nested fills: p = 18% alpha, FDR = 45% alpha, Pi = solid
  FRAC_FILL <- c()
  for (cn in names(SET_COLS)) {
    col <- unname(SET_COLS[cn])
    FRAC_FILL[paste(cn, THR_LEVELS[1], sep = "___")] <- grDevices::adjustcolor(col, alpha.f = 0.18)
    FRAC_FILL[paste(cn, THR_LEVELS[2], sep = "___")] <- grDevices::adjustcolor(col, alpha.f = 0.45)
    FRAC_FILL[paste(cn, THR_LEVELS[3], sep = "___")] <- col
  }

  # Background rect strip per contrast
  panel_bg <- tibble::tibble(
    contrast = factor(unname(CTR_LAB), levels = rev(ctr_levels)),
    fill     = unname(CONTRAST_COLORS[CORE])
  )

  x_lim_top <- 23

  # label anchored to the Π<0.05 bar's own tip: Disease/Interaction sit just past
  # it, Transplant/Rescue sit inside it (plenty of room there). Both show the plain
  # Π-score count, the subtitle already says what it is.
  lab_df <- dplyr::filter(frac_df, threshold == THR_LEVELS[3]) |>
    dplyr::mutate(
      inside = as.character(contrast) %in% unname(CTR_LAB[c("CTLvMITO", "PHEvPHE_MITO")]),
      # Disease's bar is the tightest fit of the two outside contrasts, so it gets a
      # touch more clearance than Interaction
      hjust_val = dplyr::case_when(
        as.character(contrast) == unname(CTR_LAB["CTLvPHE"]) ~ -0.45,
        !inside ~ -0.2,
        TRUE ~ 1.15
      )
    )

  ggplot2::ggplot(frac_df, ggplot2::aes(contrast, pct, fill = fill_key)) +
    ggplot2::geom_rect(
      data = panel_bg,
      ggplot2::aes(
        xmin = as.integer(contrast) - 0.5,
        xmax = as.integer(contrast) + 0.5,
        ymin = -Inf, ymax = Inf,
        fill = I(fill)
      ),
      alpha = 0.16, color = "grey70", linewidth = 0.2, inherit.aes = FALSE
    ) +
    ggplot2::geom_col(
      position  = "identity",
      width     = 0.96,
      color     = "black",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(x = contrast, y = pct, label = n, hjust = hjust_val),
      size = 1.8, fontface = "bold", color = "black", inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = FRAC_FILL) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 5, 10, 20),
      limits = c(0, x_lim_top),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.42)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "DEP counts",
      subtitle = "counts at independent p / FDR / Π thresholds",
      x = NULL, y = "% of proteome"
    ) +
    FIG_THEME +
    ggplot2::theme(
      legend.position = "none",
      plot.subtitle = ggplot2::element_text(
        size = FIG_SUBTITLE_SIZE, face = "italic", color = "grey40"
      ),
      axis.text.y = ggplot2::element_text(
        face = "bold", size = FIG_AXIS_TEXT,
        color = "grey15"
      ),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 0, 1, 2)
    )
}
