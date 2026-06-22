#!/usr/bin/env Rscript
# F04 build function — pathway enrichment, YvO Panel-F tactics.
# Per contrast: Up bar (red) and Down bar (blue), each the total significant
# pathway count (light) with the mitochondrial subset overlaid darker, both
# drawn from 0 on a sqrt axis so small bars stay readable. Faint contrast band
# behind each group; 4-key legend (Up/Down x total/mito).
# Mito = MitoCarta sets OR broad-DB pathways matching MITO_PATHWAY_REGEX.
# Returns list(plot=<ggplot>, bar_df=<df>, sig_pw=<df>, FILL_TOTAL=, FILL_MITO=).

# Total = lighter shade, mito = darker (emphasised) per direction (RdBu family).
FILL_TOTAL <- c(Up = "#F4A582", Down = "#92C5DE")
FILL_MITO <- c(Up = "#B2182B", Down = "#2166AC")
PATHWAY_KEYS <- c(
  "Up total" = "#F4A582", "Up mito" = "#B2182B",
  "Down total" = "#92C5DE", "Down mito" = "#2166AC"
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

  bar_df <- sig_pw |>
    dplyr::summarise(
      total = dplyr::n(), mito = sum(is_mito),
      .by = c(contrast, direction)
    ) |>
    tidyr::complete(
      contrast  = CORE,
      direction = c("Up", "Down"),
      fill      = list(total = 0L, mito = 0L)
    ) |>
    dplyr::mutate(
      contrast  = factor(contrast, levels = CORE),
      direction = factor(direction, levels = c("Up", "Down"))
    )

  BAR_W <- 0.34
  GAP <- 0.06
  centers <- stats::setNames(seq_along(CORE), CORE)
  bar_df <- bar_df |>
    dplyr::mutate(
      x_center = centers[as.character(contrast)] +
        ifelse(direction == "Up", -(BAR_W / 2 + GAP / 2), (BAR_W / 2 + GAP / 2)),
      key_total = paste(as.character(direction), "total"),
      key_mito = paste(as.character(direction), "mito")
    )

  # faint contrast-colour band behind each contrast's bars
  panel_bg <- tibble::tibble(
    contrast = factor(CORE, levels = CORE),
    fill     = unname(CONTRAST_COLORS[CORE])
  )
  mito_df <- dplyr::filter(bar_df, mito > 0)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = panel_bg,
      ggplot2::aes(
        xmin = as.integer(contrast) - 0.5, xmax = as.integer(contrast) + 0.5,
        ymin = 0, ymax = Inf, fill = I(fill)
      ),
      alpha = 0.16, color = "grey80", linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = bar_df,
      ggplot2::aes(
        xmin = x_center - BAR_W / 2, xmax = x_center + BAR_W / 2,
        ymin = 0, ymax = total, fill = key_total
      ),
      color = "black", linewidth = 0.3
    ) +
    ggplot2::geom_rect(
      data = mito_df,
      ggplot2::aes(
        xmin = x_center - BAR_W / 2, xmax = x_center + BAR_W / 2,
        ymin = 0, ymax = mito, fill = key_mito
      ),
      color = "black", linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(bar_df, total > 0),
      ggplot2::aes(x = x_center, y = total, label = total),
      vjust = -0.4, size = scale_text(BASE_COUNT, 60) + 0.6, fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = PATHWAY_KEYS, breaks = names(PATHWAY_KEYS), name = NULL,
      guide = ggplot2::guide_legend(nrow = 2)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(CORE),
      labels = stats::setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_sqrt(
      expand = ggplot2::expansion(mult = c(0, 0.08)),
      breaks = c(5, 25, 50, 100, 200, 300)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Pathway enrichment (Up / Down)",
      subtitle = "5-DB pool · dark = mitochondrial subset · √-scaled",
      x = NULL, y = "Significant pathways"
    ) +
    FIG_THEME +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = FIG_AXIS_TEXT, lineheight = 0.85, face = "bold"
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(2.2, "mm"),
      legend.margin = ggplot2::margin(t = -2),
      legend.box.spacing = ggplot2::unit(1, "pt"),
      plot.margin = ggplot2::margin(2, 4, 1, 2)
    )

  list(plot = p, bar_df = bar_df, sig_pw = sig_pw, FILL_TOTAL = FILL_TOTAL, FILL_MITO = FILL_MITO)
}
