# The representative term for a module is its lowest-FDR ORA hit surviving EnrichmentMap
# dedup; modules with nothing under FDR<0.10 get no term rather than a weak filler. A few
# modules read better with a curated theme than their single top set (turquoise is
# OXPHOS/energetics, not the lowest-scoring nucleobase-metabolism term).
pacman::p_load(dplyr)

MODULE_TERM_OVERRIDE <- c(turquoise = "OXPHOS / energy metabolism")

representative_terms <- function(ora_dedup) {
  ora_dedup |>
    group_by(module) |>
    slice_min(padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      module,
      database,
      term = ifelse(module %in% names(MODULE_TERM_OVERRIDE),
        MODULE_TERM_OVERRIDE[module], clean_display_label(pathway)
      )
    )
}
