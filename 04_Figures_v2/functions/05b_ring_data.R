# Library + style imports for the volcano-ring suite are loaded once by the
# stub 05_volcano_ring_plot_builder.R; do not re-source them here.

prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms = 12,
                              gap_degrees = 3,
                              start_offset = 0,
                              databases = c("Hallmark", "GO Slim")) {
  ring <- go_df |>
    filter(
      contrast == !!contrast,
      database %in% databases,
      padj < 0.05
    ) |>
    arrange(padj) |>
    slice_head(n = n_terms)

  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms for '", contrast, "'")
    return(tibble())
  }

  arc_width_deg <- (360 - n * gap_degrees) / n

  ring |>
    mutate(
      term_idx    = row_number(),
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg * pi / 180,
      mid_rad     = mid_deg * pi / 180,
      arc_r1_var  = 5.6,
      clean_label = clean_ring_label(pathway),
      gene_list   = str_split(leadingEdge, ";")
    )
}

build_tick_data <- function(ring_data,
                            de_df,
                            contrast,
                            tick_r0 = 4.4,
                            tick_r1 = 4.8) {
  if (nrow(ring_data) == 0) {
    return(tibble())
  }

  logfc_col <- paste0("logFC_", contrast)

  gene_lfc <- de_df |>
    dplyr::select(gene, lfc = all_of(logfc_col)) |>
    filter(!is.na(lfc)) |>
    distinct(gene, .keep_all = TRUE)

  pad_rad <- 0.5 * pi / 180

  map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) {
      return(tibble())
    }

    arc_start <- row$start_rad + pad_rad
    arc_end <- row$end_rad - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad

    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)

    matched_lfc <- filter(gene_lfc, gene %in% genes_in_arc)
    matched <- matched_lfc[match(genes_in_arc, matched_lfc$gene), ] |>
      filter(!is.na(gene))

    n_final <- nrow(matched)
    if (n_final == 0) {
      return(tibble())
    }
    tick_angles <- tick_angles[seq_len(n_final)]

    tibble(
      gene      = matched$gene,
      logFC     = matched$lfc,
      direction = ifelse(matched$lfc > 0, "Up", "Down"),
      angle_rad = tick_angles,
      x0        = tick_r0 * sin(tick_angles),
      y0        = tick_r0 * cos(tick_angles),
      x1        = tick_r1 * sin(tick_angles),
      y1        = tick_r1 * cos(tick_angles),
      term_idx  = row$term_idx,
      pathway   = row$pathway
    )
  })
}
