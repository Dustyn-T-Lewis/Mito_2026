# Top hubs per module by kME (real rat symbols); feeds the workbook top-hub table and
# the supplementary hub panel.
pacman::p_load(dplyr, tibble, ggplot2)

module_hub_nodes <- function(m, w, n_hub = 7L) {
  kcol <- paste0("kME_", m)
  tibble(uniprot_id = rownames(w$kME), kME = w$kME[, kcol]) |>
    left_join(distinct(w$ann, uniprot_id, gene), by = "uniprot_id") |>
    filter(w$module_colors[match(uniprot_id, names(w$module_colors))] == m, is_real_symbol(gene)) |>
    arrange(desc(kME)) |>
    slice_head(n = n_hub) |>
    mutate(label = sub("[.].*$", "", gene))
}

# Hub structure: the top proteins per module by module membership (kME), one facet per
# module, bars in the module colour. Modules run largest first.
panel_hubs <- function(w, modules, n_hub = 5L) {
  d <- bind_rows(lapply(modules, \(m) mutate(module_hub_nodes(m, w, n_hub), module = m))) |>
    mutate(
      module = factor(module, levels = modules),
      key = factor(paste(module, label, sep = "___")),
      key = reorder(key, kME)
    )
  ggplot(d, aes(kME, key, fill = as.character(module))) +
    geom_col(colour = "grey20", linewidth = 0.2, width = 0.75) +
    scale_fill_identity() +
    scale_y_discrete(labels = \(x) sub("^.*___", "", x)) +
    coord_cartesian(xlim = c(0.5, 1)) +
    facet_wrap(~module, scales = "free_y", ncol = 5) +
    labs(
      title = "Hub proteins per module",
      subtitle = sprintf("Top %d proteins by module membership (kME)", n_hub),
      x = "kME", y = NULL
    ) +
    FIG_THEME
}
