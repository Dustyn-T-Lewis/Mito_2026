# F05 build helper — sourced by 01_main_panels.R. Selects the per-contrast ring
# terms (significant -> EnrichmentMap dedup -> top RING_N) and draws each ring
# with enrichVolcano::volcano_ring(); pi_score drives the volcano, fgsea padj the
# arcs. Label cleaning lives in the package.

# Build one standalone ring; return list(plot, tag, ctr, role, terms, full).
build_ring <- function(ctr, tag, role) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep <- if (pi_col %in% names(dep_df)) {
    sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE)
  } else {
    0
  }

  sig_pool <- fgsea_all |>
    filter(
      contrast == ctr, database %in% POOL_DBS, !is.na(padj), padj < RING_PADJ,
      size >= RING_MIN_SZ, !pathway %in% MITO_DROP_SETS
    ) |>
    arrange(padj)
  n_pre <- nrow(sig_pool)
  if (n_pre > 1) {
    sig_pool <- deduplicate_enrichment(as.data.frame(sig_pool),
      pathways = SET_POOL, cutoff = SIM_CUT, cross_db = TRUE
    ) |> as_tibble()
  }
  top_terms <- slice_head(arrange(sig_pool, padj), n = RING_N)
  n_path <- nrow(top_terms)

  # Full tested-pathway table for the workbook; shown = drawn on the ring.
  full_tbl <- fgsea_all |>
    filter(contrast == ctr, database %in% POOL_DBS, !pathway %in% MITO_DROP_SETS) |>
    transmute(pathway, database, padj, pval, NES, size,
      shown = pathway %in% top_terms$pathway
    ) |>
    arrange(padj)

  volc <- dep_df |>
    transmute(
      gene,
      logFC = .data[[paste0("logFC_", ctr)]],
      P.Value = .data[[paste0("P.Value_", ctr)]],
      padj = .data[[pi_col]]
    ) |>
    filter(!is.na(logFC), !is.na(P.Value))

  p <- enrichVolcano::volcano_ring(
    volc, top_terms,
    gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
    term_col = "pathway", nes_col = "NES", size_col = "size",
    genes_col = "leadingEdge", genes_sep = ";",
    p_threshold = H9C2_PI_THRESH, logfc_threshold = log2(1.5),
    ring_radius = 5.5, volcano_radius = 5.2, arc_height_range = c(0.1, 3.2),
    title = contrast_brief(ctr),
    subtitle = sprintf(
      "%s | %d DEPs, %d pathways", CONTRAST_MATH_BRIEF[ctr], n_dep, n_path
    ),
    label_size = 2.7, point_size = 1.4, point_alpha = 0.6,
    count_x_mult = 0.55, count_y_mult = 0.55
  )

  ggsave(file.path(RPT_PNG, sprintf("MAIN_F05_%s_ring.png", tag)), p,
    width = 150, height = 150, units = "mm", dpi = 300
  )

  shown <- if (n_path > 0) {
    top_terms |>
      select(pathway, database, padj, NES, size, any_of(c("ES", "log2err"))) |>
      arrange(desc(NES))
  } else {
    top_terms[FALSE, , drop = FALSE]
  }

  message(sprintf(
    "  [%s/%s] %s: %d pre-dedup -> %d shown",
    tag, role, contrast_brief(ctr), n_pre, n_path
  ))
  list(plot = p, tag = tag, ctr = ctr, role = role, terms = shown, full = full_tbl)
}
