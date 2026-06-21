# Library + style imports for the volcano-ring suite are loaded once by the
# stub 05_volcano_ring_plot_builder.R; do not re-source them here.

make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title              = NULL,
                              contrast_title     = NULL,
                              contrast_subtitle  = NULL,
                              title_size         = 22,   # scaled from F01's 12pt @ 215mm to F03's 380mm canvas
                              subtitle_size      = NULL,
                              n_terms            = 12,
                              gap_degrees        = 3,
                              start_offset       = 0,
                              databases          = c("Hallmark", "GO Slim"),
                              volcano_radius     = 3.5,
                              tick_r0            = 4.4,
                              tick_r1            = 4.8,
                              arc_r0             = 4.8,
                              arc_r1             = 5.6,
                              label_r            = 7.0,
                              label_gap          = NULL,
                              fc_thresh          = log2(1.5),
                              p_thresh           = 0.05,
                              up_color           = DIR_COLORS["Up"],
                              down_color         = DIR_COLORS["Down"],
                              ns_color           = DIR_COLORS["NS"],
                              point_size         = 0.6,
                              point_alpha        = 0.5,
                              label_size         = 3.0,
                              count_label_size    = 2.8,
                              count_label_padding = 3,
                              count_border_width  = 0.4,
                              count_y_mult        = 1.0,
                              count_x_mult        = 0.5,
                              label_padding       = 2,
                              min_angle_gap       = 18,
                              nudge_outward       = 0.8,
                              ring_data_override = NULL,
                              bg_color           = NULL,
                              bg_alpha           = 0.12,
                              mito_pattern       = NULL,
                              show_legend        = TRUE) {

  if (is.null(subtitle_size)) subtitle_size <- title_size * 0.65

  if (!is.null(ring_data_override)) {
    ring_data <- ring_data_override
  } else {
    ring_data <- prepare_ring_data(
      go_df = go_df, contrast = contrast, n_terms = n_terms,
      gap_degrees = gap_degrees, start_offset = start_offset,
      databases = databases
    )
  }

  tick_data <- build_tick_data(
    ring_data = ring_data, de_df = de_df, contrast = contrast,
    tick_r0 = tick_r0, tick_r1 = tick_r1
  )

  volcano_layers <- build_volcano_layers(
    de_df = de_df, contrast = contrast, volcano_radius = volcano_radius,
    fc_thresh = fc_thresh, p_thresh = p_thresh,
    up_color = up_color, down_color = down_color, ns_color = ns_color,
    point_size = point_size, point_alpha = point_alpha,
    count_label_size = count_label_size,
    count_label_padding = count_label_padding,
    count_border_width = count_border_width,
    count_y_mult = count_y_mult,
    count_x_mult = count_x_mult
  )

  ring_layers <- build_ring_layers(
    ring_data = ring_data, tick_data = tick_data,
    tick_r0 = tick_r0, tick_r1 = tick_r1,
    arc_r0 = arc_r0, arc_r1 = arc_r1,
    up_color = up_color, down_color = down_color,
    mito_pattern = mito_pattern
  )

  label_layers <- build_label_layer(
    ring_data = ring_data, label_r = label_r, label_size = label_size,
    label_gap = label_gap, label_padding = label_padding,
    min_angle_gap = min_angle_gap, nudge_outward = nudge_outward,
    up_color = up_color, down_color = down_color
  )

  max_label_r <- attr(label_layers, "max_label_r")
  if (is.null(max_label_r)) max_label_r <- label_r

  # disc fill behind ring (contrast color, up to tick ring)
  bg_layer <- if (!is.null(bg_color)) {
    bg_circle <- data.frame(
      x = tick_r0 * cos(seq(0, 2 * pi, length.out = 200)),
      y = tick_r0 * sin(seq(0, 2 * pi, length.out = 200))
    )
    geom_polygon(data = bg_circle, aes(x = x, y = y),
                 fill = bg_color, alpha = bg_alpha, color = NA,
                 inherit.aes = FALSE)
  }

  legend_pos <- if (show_legend) "right" else "none"

  title_lab    <- contrast_title %||% title
  subtitle_lab <- contrast_subtitle

  p <- ggplot() +
    bg_layer +
    ring_layers$tick_bg +
    volcano_layers +
    ring_layers$ticks +
    ring_layers$enrich_arcs +
    ring_layers$fill_scale +
    ring_layers$mito_outline +
    label_layers +
    labs(title = title_lab, subtitle = subtitle_lab) +
    coord_fixed(
      xlim = c(-(max_label_r + 0.8), max_label_r + 0.8),
      # Top/bottom: extra headroom so multi-line labels at 12 and 6 o'clock
      # don't clip into the composite title strip / panel border.
      ylim = c(-(max_label_r + 1.2), max_label_r + 1.2),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = title_size,
                                       hjust = 0.5, margin = margin(b = 0, unit = "mm")),
          plot.subtitle = element_text(face = "bold.italic", size = subtitle_size,
                                       color = "grey30", hjust = 0.5,
                                       margin = margin(b = 0.5, unit = "mm")),
          plot.tag      = element_text(face = "bold", size = 26),  # scaled from F01's 15pt @ 215mm to F03's 380mm canvas
          plot.tag.position = c(0.02, 0.99),
          plot.margin   = margin(1, 1, 1, 1, "mm"),
          legend.position = legend_pos,
          legend.title = element_text(size = 7, face = "bold", color = "grey30"),
          legend.text  = element_text(size = 6, color = "grey40"),
          legend.key.width  = unit(2, "mm"),
          legend.key.height = unit(12, "mm"),
          legend.margin = margin(l = 0, r = 0)) +
    guides(color = "none",
           fill = guide_colorbar(direction = "vertical"))

  attr(p, "ring_data")  <- ring_data
  attr(p, "tick_data")  <- tick_data
  attr(p, "n_up")       <- attr(volcano_layers, "n_up")
  attr(p, "n_down")     <- attr(volcano_layers, "n_down")

  p
}

build_nes_legend_bar <- function(text_size = 10, title_size = 11,
                                 bar_margin = margin(-9, 120, 0, 120, "mm")) {
  nes_data <- data.frame(NES = seq(-3, 3, length.out = 200), y = 1)
  ggplot(nes_data, aes(x = .data$NES, y = .data$y, fill = .data$NES)) +
    geom_raster(interpolate = TRUE) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(
      colours = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
      values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
      limits  = c(-3, 3),
      guide   = "none"
    ) +
    scale_x_continuous(
      breaks = c(-3, -1.5, 0, 1.5, 3),
      labels = c("-3", "-1.5", "0", "1.5", "3"),
      expand = c(0, 0)
    ) +
    labs(x = "NES") +
    theme_void() +
    theme(
      axis.text.x  = element_text(size = text_size, face = "bold", color = "grey25"),
      axis.title.x = element_text(size = title_size, face = "bold",
                                   color = "grey25", margin = margin(t = 1)),
      plot.margin  = bar_margin
    )
}
