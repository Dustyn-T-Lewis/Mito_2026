# 04_Figures_v2/02_DEP_bars/a_script/_build.R
# Builder functions for F02 DEP bars + effect-size panels.
# Exposes: build_dep_count_panel(), build_dep_effect_panel(), dep_count_data().
# Sourced by 01_main_panels.R — never run directly.

# Narrative order: Disease first (H9C2_CONTRAST_ORDER)
CORE <- H9C2_CONTRAST_ORDER
ctr_levels <- role_label(CORE)

THR_LEVELS <- c("p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), "Π < 0.05")
ALPHAS <- c(0.30, 0.60, 1.00)

# Load data once at source time so both builders share it.
comb_long <- readr::read_csv(P05$comb, show_col_types = FALSE)
dep_results <- stats::setNames(
  lapply(CORE, \(c) as.data.frame(dplyr::filter(comb_long, contrast == c))),
  CORE
)

comb <- load_combined_wide()
n_total <- length(unique(comb$gene[!is.na(comb$gene)]))

# ---- counts tibble -----------------------------------------------------------
sig_flag <- function(r, thr) {
  switch(thr,
    "p"   = !is.na(r$P.Value) & r$P.Value < 0.05,
    "fdr" = !is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR,
    "pi"  = !is.na(r$pi_score) & r$pi_score < H9C2_PI_THRESH
  )
}

.build_counts_df <- function() {
  dplyr::bind_rows(lapply(CORE, function(c) {
    r <- dep_results[[c]]
    dplyr::bind_rows(lapply(c(p = "p", fdr = "fdr", pi = "pi"), function(thr) {
      s <- sig_flag(r, thr)
      up <- sum(s & r$logFC > 0, na.rm = TRUE)
      dn <- sum(s & r$logFC < 0, na.rm = TRUE)
      tibble::tibble(
        contrast  = role_label(c),
        threshold = thr,
        direction = c("Up", "Down"),
        n         = c(up, dn)
      )
    }), .id = NULL)
  })) |>
    dplyr::mutate(
      threshold = dplyr::recode(threshold,
        p = THR_LEVELS[1], fdr = THR_LEVELS[2], pi = THR_LEVELS[3]
      ),
      contrast = factor(contrast, levels = ctr_levels),
      threshold = factor(threshold, levels = THR_LEVELS),
      direction = factor(direction, levels = c("Up", "Down")),
      pct = 100 * n / n_total
    )
}

counts_df_cache <- .build_counts_df()

# ---- public: raw data for workbook ------------------------------------------
dep_count_data <- function() counts_df_cache

# ---- Panel A: DEP counts -----------------------------------------------------
build_dep_count_panel <- function() {
  counts_df <- counts_df_cache

  y_max <- max(counts_df$pct) * 1.12

  # Threshold faint->solid key (neutral grey swatches, inset legend)
  thr_key <- tibble::tibble(
    y    = 3:1,
    lab  = THR_LEVELS,
    fill = vapply(ALPHAS, \(a) grDevices::adjustcolor("grey25", alpha.f = a), character(1))
  )
  p_thrkey <- ggplot2::ggplot(thr_key) +
    ggplot2::geom_point(ggplot2::aes(0, y),
      shape = 22, size = 2,
      fill = thr_key$fill, color = "grey40", stroke = 0.3
    ) +
    ggplot2::geom_text(ggplot2::aes(0.25, y, label = lab),
      hjust = 0, size = 1.5, color = "grey20"
    ) +
    ggplot2::scale_x_continuous(limits = c(-0.2, 2.4)) +
    ggplot2::scale_y_continuous(limits = c(0.4, 3.6)) +
    ggplot2::theme_void()

  p_counts <- ggplot2::ggplot(
    counts_df,
    ggplot2::aes(
      x     = contrast,
      y     = pct,
      fill  = direction,
      alpha = threshold,
      group = interaction(direction, threshold)
    )
  ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge2(width = 0.9, padding = 0.1),
      color     = "grey20",
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = head(seq_along(ctr_levels), -1) + 0.5,
      color = "grey85", linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(
      values = DIR_COLORS[c("Up", "Down")],
      name   = NULL
    ) +
    ggplot2::scale_alpha_manual(
      values = stats::setNames(ALPHAS, levels(counts_df$threshold)),
      name   = NULL,
      guide  = "none"
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = "DEP counts by contrast",
      subtitle = "Up / Down dodged; p / FDR / Π faint→solid (independent thresholds)",
      x = NULL, y = "% of proteome", tag = "A"
    ) +
    FIG_THEME +
    ggplot2::theme(
      legend.position = c(0.99, 0.97),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.7), color = NA
      ),
      legend.key.size = ggplot2::unit(2.5, "mm"),
      plot.subtitle = ggplot2::element_text(
        size = FIG_SUBTITLE_SIZE, face = "italic", color = "grey40"
      ),
      axis.text.x = ggplot2::element_text(face = "bold", size = FIG_AXIS_TEXT),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 3, 1, 2)
    ) +
    patchwork::inset_element(p_thrkey,
      left = 0.62, right = 0.99,
      top = 0.78, bottom = 0.52
    )

  p_counts
}

# ---- Panel B: effect-size companion ------------------------------------------
build_dep_effect_panel <- function() {
  # Plot window: ±1; median computed on ALL logFC (no window filter in summarise)
  lfc_long <- dplyr::bind_rows(lapply(CORE, \(c)
  tibble::tibble(contrast = role_label(c), logFC = dep_results[[c]]$logFC))) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::mutate(contrast = factor(contrast, levels = ctr_levels))

  # Median over ALL proteins (unfiltered), then restrict plotting window
  lfc_stats <- lfc_long |>
    dplyr::summarise(med_abs = median(abs(logFC), na.rm = TRUE), .by = contrast) |>
    dplyr::mutate(lab = sprintf("median |log2FC| %.2f (all)", med_abs))

  lfc_plot <- lfc_long |> dplyr::filter(abs(logFC) <= 1)

  hbw <- 2 / 44

  p_eff <- ggplot2::ggplot(lfc_plot, ggplot2::aes(logFC)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
    ggplot2::geom_histogram(
      breaks    = seq(-1, 1, by = hbw),
      fill      = "grey75",
      color     = "white",
      linewidth = 0.1
    ) +
    ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(count) * hbw),
      color = "grey20", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = lfc_stats,
      ggplot2::aes(x = -0.95, y = Inf, label = lab),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.4,
      size = scale_text(BASE_STAT, 50) + 0.3,
      fontface = "bold",
      color = "grey25"
    ) +
    ggplot2::facet_wrap(~contrast,
      ncol = 1, scales = "free_y",
      strip.position = "right"
    ) +
    ggplot2::scale_x_continuous(breaks = c(-1, 0, 1)) +
    ggplot2::scale_y_continuous(breaks = NULL) +
    ggplot2::labs(
      title = "Effect size",
      subtitle = "signed log2FC (±1 shown); median over all proteins",
      x = expression(bold(log[2] ~ FC)), y = NULL, tag = "B"
    ) +
    FIG_THEME +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(
        face = "bold", size = FIG_STRIP_SIZE - 0.5, angle = 0
      ),
      strip.background = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(
        size = FIG_SUBTITLE_SIZE, face = "italic", color = "grey40"
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = FIG_AXIS_TEXT),
      panel.grid = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(1.5, "pt"),
      plot.margin = ggplot2::margin(5, 2, 1, 1)
    )

  p_eff
}
