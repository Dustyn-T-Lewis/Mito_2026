# Primary F03 building blocks. Each module is a row across four aligned panels: a
# right-to-left protein-count bar, its four 2x2 contrast effects (eigengene shift in SD
# units, FDR/nominal significance), the member-protein z trajectory, and its top-5 ORA
# pathways on an independent -log10 FDR axis. panel_module_card assembles the row.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

HEAT_CONTRASTS <- c("Disease", "Transplant", "Rescue", "Interaction")

fmt_pq <- function(x) ifelse(x < 0.01, "< 0.01", sprintf("= %.2f", x))

# Tints every panel in a module's row with that module's own colour, so all four
# columns of a row read as one unit. inherit.aes = FALSE + fill = I(...) keeps this
# independent of whatever fill scale (logFC, database) the panel's own geoms use.
module_bg_layer <- function(modules, alpha = 0.12) {
  bg <- data.frame(module = factor(modules, levels = modules))
  geom_rect(
    data = bg, aes(fill = I(as.character(module))),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = alpha, inherit.aes = FALSE
  )
}

panel_count_bars <- function(mod_size, modules) {
  d <- mod_size |>
    filter(module %in% modules) |>
    mutate(module = factor(module, levels = modules))
  ggplot(d, aes(n, 1, fill = as.character(module))) +
    geom_col(width = 0.4, colour = "black", linewidth = 0.2, orientation = "y") +
    shadowtext::geom_shadowtext(aes(label = n),
      hjust = 1.2, size = 1.8, fontface = "bold",
      colour = "black", bg.colour = "white", bg.r = 0.15
    ) +
    facet_wrap(~module, ncol = 1, strip.position = "left") +
    scale_fill_identity() +
    scale_x_reverse(expand = expansion(mult = c(0.42, 0.08))) +
    scale_y_continuous(limits = c(0.4, 1.6), expand = c(0, 0)) +
    labs(title = "Proteins (n)", subtitle = "member-protein count per module", x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT, colour = "grey15"),
      plot.subtitle = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey40"),
      strip.placement = "outside", strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = FIG_AXIS_TEXT + 1, hjust = 1, colour = "black"),
      axis.text.x = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey35"),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(colour = "grey60", linewidth = 0.2),
      panel.grid = element_blank(), panel.spacing.y = unit(0.8, "mm"),
      legend.position = "none", plot.margin = margin(2, 1, 2, 2)
    )
}

panel_contrast_bars <- function(mod_size, mod_stats) {
  modules <- arrange(mod_size, desc(n))$module
  lim <- mod_stats |>
    filter(contrast %in% HEAT_CONTRASTS) |>
    pull(logFC) |>
    abs() |>
    max()
  d <- mod_stats |>
    filter(contrast %in% HEAT_CONTRASTS, module %in% modules) |>
    mutate(
      module = factor(module, levels = modules),
      contrast = factor(contrast, levels = HEAT_CONTRASTS),
      sig_fdr = fdr < 0.10,
      star = ifelse(fdr < 0.05, "✱ ", ""),
      lab_out = sprintf("%sp %s (q %s)", star, fmt_pq(p), fmt_pq(fdr)),
      lab_in = sprintf("%sp %s\nq %s", star, fmt_pq(p), fmt_pq(fdr)),
      vj_out = ifelse(logFC >= 0, -0.45, 1.45)
    ) |>
    # a bar tall enough (relative to its own row) fits its p/q label stacked inside;
    # a short bar gets the label pushed just clear of the tip instead
    group_by(module) |>
    mutate(fits_inside = abs(logFC) >= 0.6 * max(abs(logFC))) |>
    ungroup()
  ggplot(d, aes(contrast, logFC, fill = logFC)) +
    module_bg_layer(modules) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_col(width = 0.64, colour = "grey30", linewidth = 0.2) +
    geom_col(data = filter(d, sig_fdr), width = 0.64, fill = NA, colour = "black", linewidth = 0.8) +
    geom_text(
      data = filter(d, !fits_inside),
      aes(label = lab_out, vjust = vj_out, fontface = ifelse(sig_fdr, "bold", "plain")),
      size = 1.7, colour = "black"
    ) +
    geom_text(
      data = filter(d, fits_inside),
      aes(label = lab_in, vjust = ifelse(logFC >= 0, 1.15, -0.15), fontface = ifelse(sig_fdr, "bold", "plain")),
      size = 1.6, colour = "white", lineheight = 0.85
    ) +
    facet_wrap(~module, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#1A3D6E", mid = "white", high = "#8B1A1A", midpoint = 0, limits = c(-lim, lim),
      name = "Δ module eigengene (SD units)",
      guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top", title.hjust = 0.5)
    ) +
    scale_x_discrete(expand = expansion(add = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.16)) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Eigengene response",
      subtitle = "p (q) · bold + outline FDR < 0.10 · ✱ FDR < 0.05",
      x = NULL, y = "Δ eigengene (SD units)"
    ) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT, colour = "grey15"),
      plot.subtitle = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey40"),
      strip.text = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT, face = "bold"),
      axis.text.y = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey45"),
      panel.grid = element_blank(), panel.spacing.y = unit(0.8, "mm"), legend.position = "bottom"
    )
}

panel_cluster_trajectory <- function(prot_traj, eig_traj, modules) {
  pt <- prot_traj |>
    filter(module %in% modules) |>
    mutate(module = factor(module, levels = modules), gx = match(Group, H9C2_GROUP_LEVELS))
  et <- eig_traj |>
    filter(module %in% modules) |>
    mutate(module = factor(module, levels = modules), gx = match(Group, H9C2_GROUP_LEVELS))
  ggplot() +
    module_bg_layer(modules) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_line(data = pt, aes(gx, zmean, group = uniprot_id), colour = "grey65", linewidth = 0.1, alpha = 0.2) +
    geom_line(data = et, aes(gx, zmean, colour = as.character(module)), linewidth = 1.1) +
    geom_point(data = et, aes(gx, zmean, colour = as.character(module)), size = 0.9) +
    facet_wrap(~module, ncol = 1, scales = "free_y") +
    scale_colour_identity() +
    scale_x_continuous(breaks = 1:4, labels = unname(GROUP_LABELS), expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = 0.04)) +
    labs(
      title = "Protein cluster (z)",
      subtitle = "z-scored abundance · grey = individual proteins",
      x = NULL, y = "standardized abundance"
    ) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT, colour = "grey15"),
      plot.subtitle = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey40"),
      strip.text = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = FIG_AXIS_TEXT),
      axis.text.y = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey45"),
      panel.grid = element_blank(), panel.spacing.y = unit(0.8, "mm"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 2)
    )
}

ora_aligned <- function(ora_top5, modules, top_n = 5L, sig = H9C2_FDR_EXPLOR) {
  # Name always sits outside the bar (black, unconstrained by bar length); p/(FDR)
  # always sits inside, right-aligned near the tip (white via shadowtext, so it reads
  # against every bar colour). Both turn bold at FDR < sig; the ✱ marks FDR < 0.05 on
  # top of that. Each module keeps its own x range (free_x): a scale shared across
  # modules made the weak modules' bars unreadable slivers next to turquoise's huge
  # OXPHOS hit, and the printed p/(FDR) already gives the exact number, so
  # cross-module bar-length comparison isn't needed.
  d <- ora_top5 |>
    filter(module %in% modules) |>
    group_by(module) |>
    slice_min(padj, n = top_n, with_ties = FALSE) |>
    mutate(k = row_number(), xmax = max(-log10(padj))) |>
    ungroup() |>
    mutate(
      module = factor(module, levels = modules),
      database = factor(database, levels = names(ORA_DB_COLORS)),
      lp = -log10(padj),
      yb = top_n + 1 - k,
      sig_fdr = padj < sig,
      star = ifelse(padj < 0.05, "✱ ", ""),
      name = clean_display_label(pathway),
      name = ifelse(nchar(name) > 34, paste0(substr(name, 1, 33), "…"), name),
      lab = sprintf("%sp %s (q %s)", star, fmt_pq(p), fmt_pq(padj)),
      # bounded so a short bar never pushes its inside-aligned label past x = 0
      gap = pmin(0.03 * xmax, 0.4 * lp)
    )
  ggplot(d, aes(lp, yb, fill = database)) +
    module_bg_layer(modules) +
    geom_col(width = 0.72, colour = "grey35", linewidth = 0.2, orientation = "y") +
    geom_col(
      data = filter(d, sig_fdr), width = 0.72, fill = NA,
      colour = "black", linewidth = 0.7, orientation = "y"
    ) +
    geom_text(
      aes(x = lp + gap, label = name, fontface = ifelse(sig_fdr, "bold", "plain")),
      hjust = 0, size = 1.9, colour = "black"
    ) +
    shadowtext::geom_shadowtext(
      aes(x = lp - gap, label = lab, fontface = ifelse(sig_fdr, "bold.italic", "italic")),
      hjust = 1, size = 1.7, colour = "white", bg.colour = "grey25", bg.r = 0.08
    ) +
    facet_wrap(~module, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = ORA_DB_COLORS, drop = FALSE, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32)), breaks = scales::breaks_pretty(3)) +
    scale_y_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(
      title = sprintf("Top %d ORA pathways", top_n),
      subtitle = "p (q) · bold + outline FDR < 0.10 · ✱ FDR < 0.05",
      x = expression(-log[10] ~ FDR), y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT, colour = "grey15"),
      plot.subtitle = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey40"),
      strip.text = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey35"),
      axis.ticks.x = element_line(colour = "grey60", linewidth = 0.2),
      axis.title.x = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey35"),
      panel.grid = element_blank(), panel.spacing.y = unit(0.8, "mm"),
      legend.position = "bottom", plot.margin = margin(2, 2, 2, 1)
    )
}

panel_module_card <- function(prot_traj, eig_traj, mod_size, mod_stats, ora_top5, title) {
  modules <- arrange(mod_size, desc(n))$module
  (panel_count_bars(mod_size, modules) +
    panel_contrast_bars(mod_size, mod_stats) +
    panel_cluster_trajectory(prot_traj, eig_traj, modules) +
    ora_aligned(ora_top5, modules) +
    plot_layout(widths = c(0.34, 0.66, 0.74, 1.7))) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 2, colour = "grey15"))
    )
}
