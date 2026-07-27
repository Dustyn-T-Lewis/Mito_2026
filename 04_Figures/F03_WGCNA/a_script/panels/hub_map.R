# Top hubs per module by kME (real rat symbols); feeds the workbook top-hub table.
pacman::p_load(dplyr, tibble)

module_hub_nodes <- function(m, w, n_hub = 7L) {
  kcol <- paste0("kME_", m)
  tibble(uniprot_id = rownames(w$kME), kME = w$kME[, kcol]) |>
    left_join(distinct(w$ann, uniprot_id, gene), by = "uniprot_id") |>
    filter(w$module_colors[match(uniprot_id, names(w$module_colors))] == m, is_real_symbol(gene)) |>
    arrange(desc(kME)) |>
    slice_head(n = n_hub) |>
    mutate(label = sub("[.].*$", "", gene))
}
