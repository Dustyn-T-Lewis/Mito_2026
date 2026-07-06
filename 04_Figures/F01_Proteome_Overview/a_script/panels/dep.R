# DEP count bars + effect-size histograms for F01 panels B and C.
# Exposes build_dep_count_panel(), build_dep_effect_panel(), dep_count_data().

# Narrative order: Disease -> Transplant -> Rescue -> Interaction
CORE <- H9C2_CONTRAST_ORDER
CTR_LAB <- stats::setNames(contrast_brief(CORE), CORE)
ctr_levels <- unname(CTR_LAB) # top-to-bottom narrative order

THR_LEVELS <- c(
  "p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), paste0("Π < ", H9C2_PI_THRESH)
)

# Load data once at source time so both builders share it.
comb_long <- readr::read_csv(P05$comb, show_col_types = FALSE)
dep_results <- stats::setNames(
  lapply(CORE, \(c) as.data.frame(dplyr::filter(comb_long, contrast == c))),
  CORE
)

comb <- load_combined_wide()
n_total <- length(unique(comb$gene[!is.na(comb$gene)]))

# counts tibble (YvO frac_df style)
.build_counts_df <- function() {
  dplyr::bind_rows(lapply(CORE, function(ctr) {
    r <- dep_results[[ctr]]
    n_p <- sum(!is.na(r$P.Value) & r$P.Value < 0.05, na.rm = TRUE)
    n_fdr <- sum(!is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR, na.rm = TRUE)
    n_pi <- sum(!is.na(r$pi_score) & r$pi_score < H9C2_PI_THRESH, na.rm = TRUE)
    tibble::tibble(
      contrast  = CTR_LAB[ctr],
      threshold = THR_LEVELS,
      n         = c(n_p, n_fdr, n_pi)
    )
  })) |>
    dplyr::mutate(
      contrast  = factor(contrast, levels = rev(ctr_levels)), # rev for coord_flip
      threshold = factor(threshold, levels = THR_LEVELS),
      pct       = 100 * n / n_total,
      fill_key  = paste(as.character(contrast), threshold, sep = "___")
    ) |>
    dplyr::filter(n > 0)
}

counts_df_cache <- .build_counts_df()

# public: raw data for workbook
dep_count_data <- function() counts_df_cache

# Panel B: DEP counts (horizontal nested-threshold)
build_dep_count_panel <- function() {
  frac_df <- counts_df_cache

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

# Panel C: effect-size companion (contrast-colored)
build_dep_effect_panel <- function() {
  # All logFC (no window filter) for the median stat
  lfc_long_all <- dplyr::bind_rows(lapply(CORE, \(c)
  tibble::tibble(contrast = CTR_LAB[c], logFC = dep_results[[c]]$logFC))) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::mutate(contrast = factor(contrast, levels = ctr_levels))

  # Median over ALL proteins; label clarifies scope
  lfc_stats <- lfc_long_all |>
    dplyr::summarise(med_abs = median(abs(logFC), na.rm = TRUE), .by = contrast) |>
    dplyr::mutate(lab = sprintf("median |log2FC| %.2f", med_abs))

  # Histogram restricted to ±0.8 window
  lfc_plot <- dplyr::filter(lfc_long_all, abs(logFC) <= 0.8)

  hbw <- 1.6 / 40

  ggplot2::ggplot(lfc_plot, ggplot2::aes(logFC)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = contrast),
      breaks = seq(-0.8, 0.8, by = hbw),
      color = "white",
      linewidth = 0.1,
      alpha = 0.85
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count) * hbw),
      color = "grey20", linewidth = 0.4
    ) +
    ggplot2::geom_label(
      data = lfc_stats,
      ggplot2::aes(x = -0.78, y = Inf, label = lab),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.0,
      size = scale_text(BASE_STAT, 50) + 0.3,
      fontface = "bold", color = "grey15",
      fill = scales::alpha("white", 0.8), label.size = 0,
      label.padding = ggplot2::unit(0.3, "mm")
    ) +
    ggplot2::facet_wrap(~contrast, ncol = 1, scales = "free_y", strip.position = "top") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(unname(CONTRAST_COLORS[CORE]), unname(CTR_LAB)),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = c(-0.8, 0, 0.8)) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(
      title = "Effect size",
      subtitle = "log2FC per contrast (median over all proteins)",
      x = expression(bold(log[2] ~ FC)),
      y = NULL
    ) +
    FIG_THEME +
    ggplot2::theme(
      strip.text          = ggplot2::element_blank(),
      strip.background    = ggplot2::element_blank(),
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      axis.text.x         = ggplot2::element_text(size = FIG_AXIS_TEXT),
      panel.grid          = ggplot2::element_blank(),
      panel.spacing.y     = ggplot2::unit(1.5, "pt"),
      plot.margin         = ggplot2::margin(5, 2, 1, 0)
    )
}
