# F05 build helpers — sourced by 01_main_panels.R.
# Defines build_ring() and build_nes_legend() for standalone ring panels.
# Engine lives in functions/05_volcano_ring_plot_builder.R and 05a–05e;
# nothing in this file touches those.

# Tidy a few over-long Reactome/GO labels (keys = engine-cleaned text).
LABEL_SHORTEN <- c(
  "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin\nEndocytosis",
  "Assembly Of Collagen Fibrils & Other Multimeric Structures" = "Collagen Fibril\nAssembly",
  "Collagen Chain Trimerization" = "Collagen\nTrimerization",
  "Processing Of Capped Intron Containing Pre mRNA" = "Pre-mRNA\nProcessing",
  "Respiratory Chain Complex I (Holoenzyme), Mitochondrial" = "Respiratory\nComplex I",
  "Respiratory Chain Complex I, Mitochondrial" = "Respiratory\nComplex I",
  "Mitochondrial Ribosome, Large Subunit" = "Mitoribosome\n(Large)",
  "Mitochondrial Ribosome, Small Subunit" = "Mitoribosome\n(Small)"
)

shorten_label <- function(x) {
  key <- gsub("\n", " ", x)
  out <- unname(LABEL_SHORTEN[key])
  ifelse(is.na(out), x, out)
}

# The n most significant pathways (lowest padj); the ring splits them Up/Down by NES.
pick_top <- function(pool, n) slice_head(arrange(pool, padj), n = n)

# Build one standalone ring; return list(plot, tag, ctr, role, terms, full).
build_ring <- function(ctr, tag, role) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < H9C2_PI_THRESH, na.rm = TRUE) else 0

  sig_pool <- fgsea_all |>
    filter(
      contrast == ctr, database %in% POOL_DBS, !is.na(padj), padj < RING_PADJ,
      size >= RING_MIN_SZ, !pathway %in% MITO_DROP_SETS
    ) |>
    arrange(padj)
  n_pre <- nrow(sig_pool)
  if (n_pre > 1) {
    sig_pool <- deduplicate_enrichment(as.data.frame(sig_pool),
      pathways = SET_POOL,
      cutoff = SIM_CUT, cross_db = TRUE
    ) |> as_tibble()
  }
  top_terms <- pick_top(sig_pool, RING_N)
  n_path <- nrow(top_terms)
  n_up <- sum(top_terms$NES > 0)
  n_dn <- sum(top_terms$NES < 0)

  # Full tested-pathway table for this contrast (every pooled-DB pathway + FDR);
  # shown = drawn on the ring. Filter padj < 0.05 in the workbook for the sig set.
  full_tbl <- fgsea_all |>
    filter(contrast == ctr, database %in% POOL_DBS, !pathway %in% MITO_DROP_SETS) |>
    transmute(pathway, database, padj, pval, NES, size,
      shown = pathway %in% top_terms$pathway
    ) |>
    arrange(padj)

  ring_data <- if (n_path == 0) {
    EMPTY_RING
  } else {
    build_ring_180_split(top_terms, ctr, fgsea_all, databases = POOL_DBS)
  }
  if (!is.null(ring_data) && nrow(ring_data) > 0 && "clean_label" %in% names(ring_data)) {
    ring_data$clean_label <- shorten_label(ring_data$clean_label)
  }
  adaptive_gap <- if (!is.null(ring_data) && nrow(ring_data) > 0) {
    0.7 + 0.3 * (max(ring_data$arc_r1_var, na.rm = TRUE) - 4.8) / 1.6
  } else {
    0.7
  }

  p <- make_volcano_ring(
    de_df = dep_df, go_df = fgsea_all, contrast = ctr,
    title = NULL, contrast_title = contrast_brief(ctr),
    contrast_subtitle = sprintf(
      "%s | %d DEPs, %d pathways",
      CONTRAST_MATH_BRIEF[ctr], n_dep, n_path
    ),
    databases = POOL_DBS, ring_data_override = ring_data,
    label_size = 2.7, label_gap = adaptive_gap, title_size = 7, subtitle_size = 4.5,
    point_size = 0.5, point_alpha = 0.55,
    count_label_size = scale_text(BASE_COUNT, 89) + 0.4,
    count_y_mult = 0.75, count_x_mult = 0.85,
    bg_color = "grey95", bg_alpha = 0.5,
    show_legend = FALSE
  )

  ggsave(file.path(RPT_PNG, sprintf("MAIN_F05_%s_ring.png", tag)), p,
    width = 110, height = 110, units = "mm", dpi = 300
  )

  shown <- if (n_path > 0) {
    top_terms |>
      select(pathway, database, padj, NES, size, any_of(c("ES", "log2err"))) |>
      arrange(desc(NES))
  } else {
    top_terms[FALSE, , drop = FALSE]
  }

  message(sprintf(
    "  [%s/%s] %s: %d pre-dedup -> %d shown (%d up, %d down)",
    tag, role, contrast_brief(ctr), n_pre, n_path, n_up, n_dn
  ))
  list(plot = p, tag = tag, ctr = ctr, role = role, terms = shown, full = full_tbl)
}

# Shared NES legend strip; the composite consumes this.
build_nes_legend <- function() {
  build_nes_legend_bar(
    text_size = 9, title_size = 10,
    bar_margin = margin(2, 6, 2, 6, "mm")
  )
}
