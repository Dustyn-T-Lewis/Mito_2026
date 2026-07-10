# Module hub map: every member placed by kME (hubness) vs peak-contrast log2FC.
# No cutoff decides what shows; leading hubs are labelled. module_hub_nodes feeds
# the workbook's top-hub table.
pacman::p_load(ggplot2, dplyr, tibble)

# Top hubs by kME (real rat symbols); feeds the workbook top-hub table.
module_hub_nodes <- function(m, w, n_hub = 7L) {
  kcol <- paste0("kME_", m)
  tibble(uniprot_id = rownames(w$kME), kME = w$kME[, kcol]) |>
    left_join(distinct(w$ann, uniprot_id, gene), by = "uniprot_id") |>
    filter(w$module_colors[match(uniprot_id, names(w$module_colors))] == m, is_real_symbol(gene)) |>
    arrange(desc(kME)) |>
    slice_head(n = n_hub) |>
    mutate(label = sub("[.].*$", "", gene))
}

panel_hub_map <- function(modules, w, ab, n_label = 6L) {
  mc <- w$module_colors
  ann <- distinct(w$ann, uniprot_id, gene)
  members <- names(mc)[mc %in% modules]
  d <- tibble(
    uniprot_id = members,
    module = factor(mc[members], levels = modules),
    kme = w$kME[cbind(members, paste0("kME_", mc[members]))],
    gene = ann$gene[match(members, ann$uniprot_id)]
  ) |>
    left_join(ab, by = "uniprot_id") |>
    filter(is_real_symbol(gene)) |>
    mutate(log2fc = coalesce(ab_logfc, 0), label = sub("[.].*$", "", gene))

  lim <- stats::quantile(abs(d$log2fc), 0.98, na.rm = TRUE)
  leaders <- d |>
    group_by(module) |>
    slice_max(kme, n = n_label, with_ties = FALSE) |>
    ungroup()

  ggplot(d, aes(kme, log2fc)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_point(aes(fill = log2fc),
      shape = 21, colour = "grey45", stroke = 0.12, size = 1.3, alpha = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = leaders, aes(label = label),
      size = 1.8, fontface = "bold", colour = "grey10", max.overlaps = Inf,
      segment.size = 0.12, segment.colour = "grey70", min.segment.length = 0, box.padding = 0.2
    ) +
    facet_wrap(~module, scales = "free_x", nrow = 2) +
    scale_fill_gradient2(
      low = "#1A3D6E", mid = "white", high = "#8B1A1A", midpoint = 0,
      limits = c(-lim, lim), oob = scales::squish, name = "peak-contrast log2FC",
      guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top", title.hjust = 0.5)
    ) +
    labs(
      title = "Module hub map (all members: kME vs peak-contrast log2FC)",
      x = "kME (intramodular connectivity)", y = "peak-contrast log2FC"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE, colour = "grey15"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}
