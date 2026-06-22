#!/usr/bin/env Rscript
# _build.R — builder for F03 Venn.
# Returns list(venn, strip, membership, region_counts, fit_stats).
# Source this file; call build_venn_panels().

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(eulerr)
})

build_venn_panels <- function() {
  SET_CONTRASTS <- c(
    Disease    = "CTLvPHE",
    Transplant = "CTLvMITO",
    Rescue     = "PHEvPHE_MITO"
  )
  SET_COLORS <- unname(CONTRAST_COLORS[SET_CONTRASTS])
  names(SET_COLORS) <- names(SET_CONTRASTS)

  # load long DEP table, keep the 3 contrasts
  comb <- read_csv(P05$comb, show_col_types = FALSE) |>
    filter(contrast %in% SET_CONTRASTS, !is.na(pi_score)) |>
    mutate(set = names(SET_CONTRASTS)[match(contrast, SET_CONTRASTS)])

  # significant membership per set (Pi < threshold), with logFC direction
  sig <- comb |>
    filter(pi_score < H9C2_PI_THRESH) |>
    transmute(
      uniprot_id, gene, set,
      direction = if_else(logFC >= 0, "Up", "Down")
    )

  set_lists <- split(sig$uniprot_id, sig$set)[names(SET_CONTRASTS)]
  set_lists <- lapply(set_lists, function(x) if (is.null(x)) character(0) else x)
  set_sizes <- vapply(set_lists, length, integer(1))

  # per-protein membership wide table + region key
  membership <- sig |>
    pivot_wider(
      id_cols = c(uniprot_id, gene),
      names_from = set,
      values_from = direction
    ) |>
    (\(d) {
      for (s in names(SET_CONTRASTS)) if (!s %in% names(d)) d[[s]] <- NA_character_
      d
    })() |>
    mutate(
      Disease_member = !is.na(Disease),
      Transplant_member = !is.na(Transplant),
      Rescue_member = !is.na(Rescue),
      region_key = paste0(
        if_else(Disease_member, "Disease", ""),
        if_else(Transplant_member, "+Transplant", ""),
        if_else(Rescue_member, "+Rescue", "")
      ) |>
        sub("^\\+", "", x = _)
    ) |>
    rename(
      Disease_dir    = Disease,
      Transplant_dir = Transplant,
      Rescue_dir     = Rescue
    ) |>
    select(
      uniprot_id, gene,
      Disease = Disease_member,
      Transplant = Transplant_member,
      Rescue = Rescue_member,
      Disease_dir, Transplant_dir, Rescue_dir, region_key
    ) |>
    arrange(region_key, gene)

  # 7 region counts
  region_counts <- membership |>
    count(region_key, name = "n") |>
    arrange(desc(n))

  # area-proportional Venn (eulerr)
  m <- membership
  eu_fit <- euler(c(
    "Disease"                   = sum(m$Disease & !m$Transplant & !m$Rescue),
    "Transplant"                = sum(!m$Disease & m$Transplant & !m$Rescue),
    "Rescue"                    = sum(!m$Disease & !m$Transplant & m$Rescue),
    "Disease&Transplant"        = sum(m$Disease & m$Transplant & !m$Rescue),
    "Disease&Rescue"            = sum(m$Disease & !m$Transplant & m$Rescue),
    "Transplant&Rescue"         = sum(!m$Disease & m$Transplant & m$Rescue),
    "Disease&Transplant&Rescue" = sum(m$Disease & m$Transplant & m$Rescue)
  ), shape = "ellipse")

  # capture fit quality
  fit_stats <- tibble(
    metric = c("stress", "diagError"),
    value  = c(eu_fit$stress, eu_fit$diagError)
  )

  set_labels <- names(SET_CONTRASTS)

  venn_grob <- plot(
    eu_fit,
    fills = list(fill = SET_COLORS, alpha = 0.5),
    edges = list(col = SET_COLORS, lwd = 1.2),
    labels = list(labels = set_labels, fontsize = 6, fontfamily = "Helvetica", font = 2),
    quantities = list(fontsize = 6, fontfamily = "Helvetica"),
    legend = FALSE
  )

  # companion Up/Down directional strip
  strip_dat <- sig |>
    count(set, direction, name = "n") |>
    mutate(group = factor(set, levels = names(SET_CONTRASTS)))

  core_ids <- m$uniprot_id[m$Disease & m$Transplant & m$Rescue]
  if (length(core_ids)) {
    core_dir <- sig |>
      filter(uniprot_id %in% core_ids) |>
      distinct(uniprot_id, set, direction) |>
      count(direction, name = "n") |>
      mutate(group = "Core (all 3)", set = "Core (all 3)")
    strip_dat <- bind_rows(strip_dat, core_dir)
  }
  strip_levels <- c(names(SET_CONTRASTS), if (length(core_ids)) "Core (all 3)")
  strip_dat$group <- factor(strip_dat$group, levels = strip_levels)
  strip_dat$direction <- factor(strip_dat$direction, levels = c("Up", "Down"))

  # faint set-colour band behind each group (Core left neutral)
  strip_bg <- tibble(
    group = factor(strip_levels, levels = strip_levels),
    fill  = c(SET_COLORS[names(SET_CONTRASTS)], if (length(core_ids)) "grey70")
  )

  p_strip <- ggplot(strip_dat, aes(group, n, fill = direction)) +
    geom_rect(
      data = strip_bg,
      aes(
        xmin = as.integer(group) - 0.5, xmax = as.integer(group) + 0.5,
        ymin = -Inf, ymax = Inf, fill = I(fill)
      ),
      alpha = 0.13, inherit.aes = FALSE
    ) +
    geom_col(
      position = position_dodge(preserve = "single", width = 0.8),
      width = 0.7, color = "grey25", linewidth = 0.15
    ) +
    geom_text(
      aes(label = n),
      position = position_dodge(preserve = "single", width = 0.8),
      vjust = -0.3, size = 1.6, color = "grey15"
    ) +
    scale_fill_manual(values = DIR_COLORS[c("Up", "Down")], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Direction within set", x = NULL, y = "Proteins (Π < 0.05)") +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT),
      legend.position = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.key.size = unit(2.5, "mm"),
      legend.margin = margin(1, 1, 1, 1),
      plot.margin = margin(3, 2, 1, 1)
    )

  # wrap venn_grob as ggplot if ggplotify is available
  have_ggplotify <- requireNamespace("ggplotify", quietly = TRUE)
  sub_txt <- sprintf(
    "Disease: %s  |  Transplant: %s  |  Rescue: %s",
    CONTRAST_MATH_BRIEF[["CTLvPHE"]],
    CONTRAST_MATH_BRIEF[["CTLvMITO"]],
    CONTRAST_MATH_BRIEF[["PHEvPHE_MITO"]]
  )
  if (have_ggplotify) {
    venn_gg <- ggplotify::as.ggplot(venn_grob) +
      labs(title = "DEP overlap (Π < 0.05)", subtitle = sub_txt) +
      theme_void(base_family = "Helvetica") +
      theme(
        plot.title = element_text(
          face = "bold", size = FIG_TITLE_SIZE,
          margin = margin(b = 1)
        ),
        plot.subtitle = element_text(
          face = "bold.italic", size = FIG_SUBTITLE_SIZE,
          color = "grey30", margin = margin(b = 2)
        ),
        plot.margin = margin(3, 2, 1, 2)
      )
  } else {
    venn_gg <- venn_grob
  }

  message(sprintf(
    "F03 Venn | sets: Disease=%d Transplant=%d Rescue=%d | %d sig proteins, %d regions",
    set_sizes[["Disease"]], set_sizes[["Transplant"]], set_sizes[["Rescue"]],
    nrow(membership), nrow(region_counts)
  ))

  list(
    venn = venn_gg,
    strip = p_strip,
    membership = membership,
    region_counts = region_counts,
    fit_stats = fit_stats
  )
}
