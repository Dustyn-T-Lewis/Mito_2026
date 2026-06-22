# Library + style imports for the volcano-ring suite are loaded once by the
# stub 05_volcano_ring_plot_builder.R; do not re-source them here.

build_volcano_layers <- function(de_df,
                                 contrast,
                                 volcano_radius = 3.5,
                                 fc_thresh = log2(1.5),
                                 p_thresh = 0.05,
                                 up_color = DIR_COLORS["Up"],
                                 down_color = DIR_COLORS["Down"],
                                 ns_color = DIR_COLORS["NS"],
                                 point_size = 0.6,
                                 point_alpha = 0.5,
                                 count_label_size = 2.8,
                                 count_label_padding = 3,
                                 count_border_width = 0.4,
                                 count_y_mult = 1.0,
                                 count_x_mult = 0.5) {
  logfc_col <- paste0("logFC_", contrast)
  pval_col <- paste0("P.Value_", contrast)
  pi_col <- paste0("pi_score_", contrast)

  vdf <- de_df |>
    transmute(
      gene       = gene,
      logFC      = .data[[logfc_col]],
      pvalue     = .data[[pval_col]],
      pi_score   = .data[[pi_col]],
      neg_log10p = -log10(pvalue)
    ) |>
    filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p)) |>
    mutate(
      direction = case_when(
        pi_score < H9C2_PI_THRESH & abs(logFC) > fc_thresh & logFC > 0 ~ "Up",
        pi_score < H9C2_PI_THRESH & abs(logFC) > fc_thresh & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  n_up <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")

  x_data_max <- max(abs(vdf$logFC), na.rm = TRUE)
  y_data_max <- max(vdf$neg_log10p, na.rm = TRUE)

  margin <- 0.92
  vr <- volcano_radius * margin

  vdf <- vdf |>
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  vdf_ns <- vdf |> filter(direction == "NS")
  vdf_sig <- vdf |> filter(direction != "NS")

  layers <- list(
    ns_points = geom_point(
      data = vdf_ns, aes(x = x_plot, y = y_plot),
      color = ns_color, size = point_size * 0.8, alpha = point_alpha * 0.35,
      inherit.aes = FALSE
    ),
    sig_points = geom_point(
      data = vdf_sig, aes(x = x_plot, y = y_plot, color = direction),
      size = point_size * 1.4, alpha = point_alpha * 1.4, stroke = 0.3,
      inherit.aes = FALSE
    ),
    color_scale = scale_color_manual(
      values = c(Up = unname(up_color), Down = unname(down_color)),
      guide  = "none"
    ),
    x_axis_line = annotate(
      "segment",
      x = -vr * 0.42, xend = vr * 0.42, y = -vr, yend = -vr,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "both", length = unit(1.2, "mm"), type = "closed")
    ),
    x_axis_up = annotate(
      "text",
      x = vr * 0.45, y = -vr,
      label = "up", size = count_label_size * 0.9,
      color = unname(up_color), fontface = "bold.italic", hjust = 0
    ),
    x_axis_down = annotate(
      "text",
      x = -vr * 0.45, y = -vr,
      label = "down", size = count_label_size * 0.9,
      color = unname(down_color), fontface = "bold.italic", hjust = 1
    ),
    x_axis_label = annotate(
      "text",
      x = 0, y = -vr - 0.35,
      label = "log2 FC", size = count_label_size * 0.9,
      color = "grey40", fontface = "bold.italic"
    ),
    y_axis_line = annotate(
      "segment",
      x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "last", length = unit(1.2, "mm"), type = "closed")
    ),
    y_axis_label = annotate(
      "text",
      x = 0, y = vr * 1.04,
      label = expression(bold(-log[10]) ~ bolditalic(p)), size = count_label_size * 0.9,
      color = "grey40"
    ),
    n_up_box = annotate(
      "label",
      x = vr * count_x_mult, y = vr * count_y_mult,
      label = n_up, size = count_label_size * 1.25,
      color = "black", fill = alpha(up_color, 0.9), fontface = "bold",
      label.padding = unit(count_label_padding, "pt"), label.r = unit(2, "pt"),
      linewidth = count_border_width
    ),
    n_up_text = annotate(
      "text",
      x = vr * count_x_mult, y = vr * count_y_mult,
      label = n_up, size = count_label_size * 1.25,
      color = "white", fontface = "bold"
    ),
    n_down_box = annotate(
      "label",
      x = -vr * count_x_mult, y = vr * count_y_mult,
      label = n_down, size = count_label_size * 1.25,
      color = "black", fill = alpha(down_color, 0.9), fontface = "bold",
      label.padding = unit(count_label_padding, "pt"), label.r = unit(2, "pt"),
      linewidth = count_border_width
    ),
    n_down_text = annotate(
      "text",
      x = -vr * count_x_mult, y = vr * count_y_mult,
      label = n_down, size = count_label_size * 1.25,
      color = "white", fontface = "bold"
    )
  )

  attr(layers, "x_data_max") <- x_data_max
  attr(layers, "y_data_max") <- y_data_max
  attr(layers, "n_up") <- n_up
  attr(layers, "n_down") <- n_down

  layers
}
