#!/usr/bin/env Rscript
# F04 build function — significant pathways per contrast, stacked by source DB.
# Diverging: Up pathways above zero, Down below; fill = database (5-DB key).
# Faint contrast-colour band behind each contrast group.
# Returns list(plot=<ggplot>, bar_df=<df>, sig_pw=<df>, DB_COLORS=<named vec>).
# Caller handles ggsave and workbook.

# Database key (Dark2 hues; stable order = CANONICAL_DBS).
PATHWAY_DB_COLORS <- c(
  Hallmark = "#1B9E77", Reactome = "#7570B3", KEGG = "#E7298A",
  MitoCarta = "#D95F02", `GO Slim` = "#66A61E"
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

  db_levels <- names(PATHWAY_DB_COLORS)

  # counts per contrast x direction x database; Down side drawn negative
  bar_df <- sig_pw |>
    dplyr::summarise(n = dplyr::n(), .by = c(contrast, direction, database)) |>
    tidyr::complete(
      contrast  = CORE,
      direction = c("Up", "Down"),
      database  = db_levels,
      fill      = list(n = 0L)
    ) |>
    dplyr::mutate(
      contrast = factor(contrast, levels = CORE),
      database = factor(database, levels = db_levels),
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

  # faint contrast-colour band behind each contrast's bars
  panel_bg <- tibble::tibble(
    contrast = factor(CORE, levels = CORE),
    fill     = unname(CONTRAST_COLORS[CORE])
  )

  p <- ggplot2::ggplot(bar_df, ggplot2::aes(contrast, y, fill = database)) +
    ggplot2::geom_rect(
      data = panel_bg,
      ggplot2::aes(
        xmin = as.integer(contrast) - 0.5, xmax = as.integer(contrast) + 0.5,
        ymin = -Inf, ymax = Inf, fill = I(fill)
      ),
      alpha = 0.13, inherit.aes = FALSE
    ) +
    ggplot2::geom_col(width = 0.74, color = "grey25", linewidth = 0.12) +
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
    ggplot2::scale_fill_manual(values = PATHWAY_DB_COLORS, name = NULL, drop = FALSE) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_y_continuous(
      labels = abs,
      expand = ggplot2::expansion(mult = c(0.12, 0.12))
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
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(2.4, "mm"),
      legend.margin = ggplot2::margin(t = -2),
      legend.box.spacing = ggplot2::unit(1, "pt"),
      plot.margin = ggplot2::margin(2, 4, 1, 2)
    )

  list(plot = p, bar_df = bar_df, sig_pw = sig_pw, DB_COLORS = PATHWAY_DB_COLORS)
}
