# Library + style imports for the volcano-ring suite are loaded once by the
# stub 05_volcano_ring_plot_builder.R; do not re-source them here.

build_ring_layers <- function(ring_data,
                              tick_data,
                              tick_r0    = 4.4,
                              tick_r1    = 4.8,
                              arc_r0     = 4.8,
                              arc_r1     = 5.6,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"],
                              mito_pattern = NULL) {
  # mito_pattern: optional regex (perl-style). When provided, arcs whose
  # `pathway` name matches the regex OR starts with "MITOCARTA_" get a second
  # bold-dark outline layer drawn on top of the standard NES-filled arc.
  # Visual cue: "this term is mitochondrial biology".

  if (nrow(ring_data) == 0) return(list())

  layers <- list()

  layers$tick_bg <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
        start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15,
    inherit.aes = FALSE
  )

  if (nrow(tick_data) > 0) {
    layers$ticks <- geom_segment(
      data = tick_data,
      aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  layers$enrich_arcs <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1_var,
        start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )

  layers$fill_scale <- scale_fill_gradientn(
    colours = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
    values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
    limits  = c(-3, 3),
    oob     = scales::squish,
    name    = "NES"
  )

  # Optional second layer: bold dark outline on mito-relevant arcs.
  if (!is.null(mito_pattern) && "pathway" %in% names(ring_data)) {
    is_mito <- grepl(mito_pattern, ring_data$pathway, perl = TRUE) |
               startsWith(ring_data$pathway, "MITOCARTA_")
    mito_arcs <- ring_data[is_mito, , drop = FALSE]
    if (nrow(mito_arcs) > 0) {
      layers$mito_outline <- geom_arc_bar(
        data = mito_arcs,
        aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1_var,
            start = start_rad, end = end_rad),
        fill = NA, color = "grey10", linewidth = 0.8,
        inherit.aes = FALSE
      )
    }
  }

  layers
}

build_label_layer <- function(ring_data,
                              label_r    = 7.0,
                              label_size = 3.0,
                              label_gap  = NULL,
                              label_padding = 2,
                              min_angle_gap = 18,
                              nudge_outward = 0.8,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"]) {

  if (nrow(ring_data) == 0) return(list())

  lbl_df <- ring_data |>
    mutate(
      label_r_term = if (!is.null(label_gap)) arc_r1_var + label_gap else label_r
    )

  # nudge close labels outward
  if (nrow(lbl_df) >= 2) {
    lbl_df <- lbl_df  |> arrange(mid_deg)
    for (i in 2:nrow(lbl_df)) {
      if (abs(lbl_df$mid_deg[i] - lbl_df$mid_deg[i - 1]) < min_angle_gap) {
        lbl_df$label_r_term[i] <- lbl_df$label_r_term[i] + nudge_outward
      }
    }
    wrap_gap <- 360 - lbl_df$mid_deg[nrow(lbl_df)] + lbl_df$mid_deg[1]
    if (wrap_gap < min_angle_gap) {
      lbl_df$label_r_term[1] <- lbl_df$label_r_term[1] + 0.8
    }
  }

  lbl_df <- lbl_df |>
    mutate(
      lbl_x      = label_r_term * sin(mid_rad),
      lbl_y      = label_r_term * cos(mid_rad),
      lead_x     = (arc_r1_var + 0.1) * sin(mid_rad),
      lead_y     = (arc_r1_var + 0.1) * cos(mid_rad),
      lead_end_x = (label_r_term - 0.25) * sin(mid_rad),
      lead_end_y = (label_r_term - 0.25) * cos(mid_rad),
      side_x     = sin(mid_rad),
      lbl_hjust  = case_when(
        side_x >  0.15 ~ 0,
        side_x < -0.15 ~ 1,
        TRUE           ~ 0.5
      ),
      nudge_x = case_when(
        side_x >  0.15 ~  0.6,
        side_x < -0.15 ~ -0.6,
        TRUE           ~  0
      ),
      legend_label = clean_label
    )

  attr(lbl_df, "max_label_r") <- max(lbl_df$label_r_term)

  up_lbl   <- lbl_df  |> filter(NES > 0)
  down_lbl <- lbl_df  |> filter(NES <= 0)

  layers <- list()

  layers$leaders <- geom_segment(
    data = lbl_df,
    aes(x = lead_x, y = lead_y, xend = lead_end_x, yend = lead_end_y),
    linewidth = 0.5, color = "grey35",
    inherit.aes = FALSE
  )

  if (nrow(up_lbl) > 0) {
    layers$up_labels <- geom_label(
      data = up_lbl,
      aes(x = lbl_x + nudge_x, y = lbl_y, label = legend_label),
      hjust = 0.5, vjust = 0.5,
      fill = unname(up_color), color = "white",
      fontface = "bold", size = label_size,
      label.padding = unit(label_padding, "pt"), label.r = unit(1.5, "pt"),
      lineheight = 0.85,
      inherit.aes = FALSE
    )
  }

  if (nrow(down_lbl) > 0) {
    layers$down_labels <- geom_label(
      data = down_lbl,
      aes(x = lbl_x + nudge_x, y = lbl_y, label = legend_label),
      hjust = 0.5, vjust = 0.5,
      fill = unname(down_color), color = "white",
      fontface = "bold", size = label_size,
      label.padding = unit(label_padding, "pt"), label.r = unit(1.5, "pt"),
      lineheight = 0.85,
      inherit.aes = FALSE
    )
  }

  attr(layers, "max_label_r") <- attr(lbl_df, "max_label_r")
  layers
}

# min_size excludes small gene sets prone to tissue-irrelevant GO artifacts
# (Reimand et al. 2019 Nat Protocols S3.4); n_each = NULL passes all sig terms.
select_ring_terms <- function(go_df, contrast_name, n_each = NULL,
                              databases = c("Hallmark", "GO Slim"),
                              min_size = 15) {
  sig <- go_df |>
    filter(contrast == contrast_name, database %in% databases,
           padj < 0.05, size >= min_size) |>
    arrange(padj)

  up   <- sig  |> filter(NES > 0)
  down <- sig  |> filter(NES < 0)

  if (!is.null(n_each)) {
    up   <- up    |> slice_head(n = n_each)
    down <- down  |> slice_head(n = n_each)
  }

  bind_rows(up, down)
}

center_ring_angles <- function(ring, n_up) {
  n <- nrow(ring)
  if (n < 2 || n_up < 1) return(ring)
  up_mid <- (ring$start_deg[1] + ring$end_deg[min(n_up, n)]) / 2
  offset <- 90 - up_mid   # center Up block at 90° (right side)
  ring$start_deg <- ring$start_deg + offset
  ring$end_deg   <- ring$end_deg   + offset
  ring$mid_deg   <- ring$mid_deg   + offset
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring
}

# YvO-style ring builder. All arcs are equal width (computed so n_total arcs +
# gaps cover the circumference). Up cluster centered at 3 o'clock (90°),
# Down cluster centered at 9 o'clock (270°). Whichever direction has more
# arcs naturally "encroaches" past 12 and 6 o'clock into the other side.
# Most-significant arc at the TOP of each cluster (toward 12 o'clock), least
# at the BOTTOM (toward 6 o'clock). Single-direction case (n_up=0 or n_down=0)
# spans the full circle starting at 12 o'clock going clockwise.
# Drop-in replacement for build_ring_with_gaps; called by Mito _build_volcano_panel.R.
build_ring_180_split <- function(top_terms, contrast_name, go_df,
                                 databases  = c("Hallmark", "Reactome"),
                                 gap_intra  = 3,    # gap between arcs within a cluster
                                 gap_split  = 8,    # extra gap at 12 and 6 o'clock between Up/Down
                                 min_height = 0.05,
                                 max_height = 1.6) {
  if (nrow(top_terms) == 0) return(tibble())

  # top_terms inherits all columns from the fgsea cache (incl. leadingEdge).
  # If a caller passes a trimmed top_terms without leadingEdge, fall back to
  # looking it up from go_df.
  ring <- top_terms
  if (!"leadingEdge" %in% names(ring)) {
    le_lookup <- go_df |>
      dplyr::filter(contrast == !!contrast_name,
                    pathway %in% ring$pathway) |>
      dplyr::distinct(pathway, .keep_all = TRUE) |>
      dplyr::select(pathway, leadingEdge)
    ring <- dplyr::left_join(ring, le_lookup, by = "pathway")
  }
  ring <- ring |>
    dplyr::mutate(
      clean_label = clean_ring_label(pathway),
      gene_list   = stringr::str_split(leadingEdge, ";")
    )

  up   <- ring |> dplyr::filter(NES >  0) |> dplyr::arrange(padj)
  down <- ring |> dplyr::filter(NES <= 0) |> dplyr::arrange(padj)
  n_up   <- nrow(up)
  n_down <- nrow(down)

  scale_arc_height <- function(df) {
    if (nrow(df) == 0) return(numeric(0))
    if (nrow(df) == 1) return(5.6 + (min_height + max_height) / 2)
    neg_lp <- -log10(pmax(df$padj, .Machine$double.xmin))
    rng    <- range(neg_lp)
    if (diff(rng) <= 0) return(rep(5.6 + (min_height + max_height) / 2, nrow(df)))
    scaled <- (neg_lp - rng[1]) / diff(rng)
    5.6 + min_height + (max_height - min_height) * sqrt(scaled)
  }
  if (n_up   > 0) up$arc_r1_var   <- scale_arc_height(up)
  if (n_down > 0) down$arc_r1_var <- scale_arc_height(down)

  # Compute arc width: all arcs (Up + Down) share the same width, sized so the
  # complete set + gaps covers the circumference.
  if (n_up > 0 && n_down > 0) {
    # Bilateral case: 2 split gaps (at 12 and 6 o'clock) + intra gaps
    total_gap <- 2 * gap_split + (n_up - 1) * gap_intra + (n_down - 1) * gap_intra
    arc_w     <- (360 - total_gap) / (n_up + n_down)

    # Up cluster centered at 90°
    up_span    <- n_up * arc_w + (n_up - 1) * gap_intra
    up_start   <- 90 - up_span / 2
    up$start_deg <- up_start + (seq_len(n_up) - 1) * (arc_w + gap_intra)
    up$end_deg   <- up$start_deg + arc_w

    # Down cluster centered at 270°. Reverse so most-sig sits near top-left.
    down_span    <- n_down * arc_w + (n_down - 1) * gap_intra
    down_start   <- 270 - down_span / 2
    down_offset  <- (n_down - seq_len(n_down)) * (arc_w + gap_intra)
    down$start_deg <- down_start + down_offset
    down$end_deg   <- down$start_deg + arc_w
  } else {
    # Single-direction case.
    only   <- if (n_up > 0) up else down
    n_only <- nrow(only)
    if (n_only == 1) {
      # Compact 30° arc at the appropriate side (don't sprawl across the circle
      # for a lone term — it would suggest a much stronger signal than exists).
      cluster_center <- if (n_up > 0) 90 else 270
      arc_w <- 30
      only$start_deg <- cluster_center - arc_w / 2
      only$end_deg   <- only$start_deg + arc_w
    } else {
      # 2+ same-direction arcs fill the whole circle, most-sig at 12 o'clock.
      arc_w <- (360 - n_only * gap_intra) / n_only
      only$start_deg <- (seq_len(n_only) - 1) * (arc_w + gap_intra)
      only$end_deg   <- only$start_deg + arc_w
    }
    if (n_up > 0) up <- only else down <- only
  }

  ring <- dplyr::bind_rows(up, down)
  ring$mid_deg   <- (ring$start_deg + ring$end_deg) / 2
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring$term_idx  <- seq_len(nrow(ring))
  ring
}


build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = NULL,
                                 databases = c("Hallmark", "GO Slim")) {
  real_rows <- go_df |>
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  padj_lookup <- real_rows |>
    dplyr::select(pathway, padj) |>
    distinct(pathway, .keep_all = TRUE)
  go_subset <- real_rows |>
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3, start_offset = 0,
    databases = databases
  )

  ring$padj <- padj_lookup$padj[match(ring$pathway, padj_lookup$pathway)]

  n <- nrow(ring)
  n_up <- sum(ring$NES > 0)

  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    if (n_up > 0 && n_up < n) gaps[n_up] <- gap_split
    gaps[n] <- gap_split
    arc_budget <- 360 - sum(gaps)
    arc_widths <- rep(arc_budget / n, n)
    min_height <- 0.05; max_height <- 1.6
    neg_lp <- -log10(pmax(ring$padj, .Machine$double.xmin))
    scaled <- (neg_lp - min(neg_lp)) / (max(neg_lp) - min(neg_lp))
    ring$arc_r1_var <- 4.8 + min_height +
      (max_height - min_height) * sqrt(scaled)
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_widths[i - 1] + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_widths[i]
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
    ring <- center_ring_angles(ring, n_up)
  }
  ring
}

