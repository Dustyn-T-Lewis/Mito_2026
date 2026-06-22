#!/usr/bin/env Rscript
# F04 build function — significant pathways per contrast, stacked by source DB.
# Diverging: Up above zero, Down below; segments = the 5 databases (not pooled),
# ordered by overall total (largest at the base). Light fills, black borders,
# swatch key in the upper-right.
# Returns list(plot=<ggplot>, bar_df=<df>, sig_pw=<df>, DB_COLORS=<named vec>).

# Database key (Set2 pastels; GO Slim kept clearly distinct from Hallmark).
PATHWAY_DB_COLORS <- c(
  Hallmark = "#66C2A5", Reactome = "#8DA0CB", KEGG = "#E78AC3",
  MitoCarta = "#FC8D62", `GO Slim` = "#FFD92F"
)

build_pathway_bar_panel <- function() {
  CORE <- H9C2_CONTRAST_ORDER # Disease-first
  PADJ_CUT <- 0.05
  MIN_SIZE <- 10
  SIM_CUT <- 0.375

  fgsea_all <- readr::read_csv(
    here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
    show_col_types = FALSE
  )
  rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
  set_pool <- do.call(c, unname(rat_gene_sets[CANONICAL_DBS]))

  per_contrast <- function(ctr) {
    sig <- fgsea_all |>
      dplyr::filter(
        contrast == ctr, database %in% CANONICAL_DBS,
        !is.na(padj), padj < PADJ_CUT,
        size >= MIN_SIZE, !pathway %in% MITO_DROP_SETS
      ) |>
      dplyr::arrange(padj)
    if (nrow(sig) > 1) {
      sig <- deduplicate_enrichment(
        as.data.frame(sig),
        pathways = set_pool,
        cutoff = SIM_CUT, cross_db = TRUE
      ) |> tibble::as_tibble()
    }
    sig |>
      dplyr::mutate(
        direction = dplyr::if_else(NES > 0, "Up", "Down"),
        is_mito = database == "MitoCarta" |
          grepl(MITO_PATHWAY_REGEX, pathway, perl = TRUE),
        contrast = ctr
      )
  }

  sig_pw <- dplyr::bind_rows(lapply(CORE, per_contrast))
  message(sprintf(
    "F04: %d total post-dedup sig pathways across %d contrasts",
    nrow(sig_pw), length(CORE)
  ))

  # databases ordered by overall total (largest first); base of the stack
  db_order <- sig_pw |>
    dplyr::count(database, sort = TRUE) |>
    dplyr::pull(database)
  db_order <- c(db_order, setdiff(names(PATHWAY_DB_COLORS), db_order))

  # counts per contrast x direction x database; Down side drawn negative
  bar_df <- sig_pw |>
    dplyr::summarise(n = dplyr::n(), .by = c(contrast, direction, database)) |>
    tidyr::complete(
      contrast = CORE, direction = c("Up", "Down"), database = db_order,
      fill = list(n = 0L)
    ) |>
    dplyr::mutate(
      contrast = factor(contrast, levels = CORE),
      database = factor(database, levels = db_order),
      y = dplyr::if_else(direction == "Down", -as.numeric(n), as.numeric(n))
    )

  # per-contrast Up/Down totals for stack-top labels
  tot_df <- bar_df |>
    dplyr::summarise(
      up   = sum(n[direction == "Up"]),
      down = sum(n[direction == "Down"]),
      .by  = contrast
    ) |>
    dplyr::mutate(contrast = factor(contrast, levels = CORE))

  p <- ggplot2::ggplot(bar_df, ggplot2::aes(contrast, y, fill = database)) +
    ggplot2::geom_col(
      width = 0.74, color = "black", linewidth = 0.2,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey35") +
    ggplot2::geom_text(
      data = tot_df, ggplot2::aes(contrast, up, label = up),
      inherit.aes = FALSE, vjust = -0.4, size = 2.2, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(tot_df, down > 0),
      ggplot2::aes(contrast, -down, label = down),
      inherit.aes = FALSE, vjust = 1.3, size = 2.2, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = PATHWAY_DB_COLORS, breaks = db_order, name = NULL,
      guide = ggplot2::guide_legend(ncol = 1)
    ) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_continuous(
      labels = abs,
      expand = ggplot2::expansion(mult = c(0.12, 0.14))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Pathway enrichment by database",
      subtitle = "stacked by source DB · Up above / Down below zero",
      x = NULL, y = "Significant pathways"
    ) +
    FIG_THEME +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = FIG_AXIS_TEXT, lineheight = 0.85, face = "bold"
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), color = "grey70", linewidth = 0.3
      ),
      legend.key.size = ggplot2::unit(2.4, "mm"),
      legend.margin = ggplot2::margin(1, 2, 1, 1),
      plot.margin = ggplot2::margin(2, 4, 1, 2)
    )

  list(plot = p, bar_df = bar_df, sig_pw = sig_pw, DB_COLORS = PATHWAY_DB_COLORS)
}
