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

fmt_pq <- function(x) ifelse(x < 0.01, "<.01", sub("^0", "", sprintf("%.2f", x)))

panel_count_bars <- function(mod_size, modules) {
  d <- mod_size |>
    filter(module %in% modules) |>
    mutate(module = factor(module, levels = modules))
  ggplot(d, aes(n, 1, fill = as.character(module))) +
    geom_col(width = 0.4, colour = "black", linewidth = 0.2, orientation = "y") +
    shadowtext::geom_shadowtext(aes(label = n),
      hjust = 1.2, size = 1.8, fontface = "bold",
      colour = "white", bg.colour = "grey15", bg.r = 0.12
    ) +
    facet_wrap(~module, ncol = 1, strip.position = "left") +
    scale_fill_identity() +
    scale_x_reverse(expand = expansion(mult = c(0.42, 0.02))) +
    scale_y_continuous(limits = c(0.4, 1.6), expand = c(0, 0)) +
    labs(title = "Proteins", x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT + 1, colour = "grey15"),
      strip.placement = "outside", strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = FIG_AXIS_TEXT + 1, hjust = 1, colour = "grey10"),
      axis.text = element_blank(), axis.ticks = element_blank(),
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
      star = ifelse(fdr < 0.10, "✱", ifelse(p < 0.10, "·", "")),
      lab = sprintf("%s q%s p%s", star, fmt_pq(fdr), fmt_pq(p)),
      vj = ifelse(logFC >= 0, -0.3, 1.3)
    )
  ggplot(d, aes(contrast, logFC, fill = logFC)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_col(width = 0.64, colour = "grey30", linewidth = 0.2) +
    geom_col(data = filter(d, sig_fdr), width = 0.64, fill = NA, colour = "black", linewidth = 0.8) +
    geom_text(aes(label = lab, vjust = vj), size = 1.3, fontface = "bold", colour = "grey20") +
    facet_wrap(~module, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#1A3D6E", mid = "white", high = "#8B1A1A", midpoint = 0, limits = c(-lim, lim),
      name = "Δ module eigengene (SD units)",
      guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top", title.hjust = 0.5)
    ) +
    scale_x_discrete(expand = expansion(add = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.32)) +
    coord_cartesian(clip = "off") +
    labs(title = "Eigengene response", x = NULL, y = "Δ eigengene (SD units)") +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT + 1, colour = "grey15"),
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
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_line(data = pt, aes(gx, zmean, group = uniprot_id), colour = "grey65", linewidth = 0.1, alpha = 0.2) +
    geom_line(data = et, aes(gx, zmean, colour = as.character(module)), linewidth = 1.1) +
    geom_point(data = et, aes(gx, zmean, colour = as.character(module)), size = 0.9) +
    facet_wrap(~module, ncol = 1, scales = "free_y") +
    scale_colour_identity() +
    scale_x_continuous(breaks = 1:4, labels = unname(GROUP_LABELS), expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = 0.04)) +
    labs(title = "Protein cluster (z)", x = NULL, y = "standardized abundance") +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT + 1, colour = "grey15"),
      strip.text = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = FIG_AXIS_TEXT),
      axis.text.y = element_text(size = FIG_AXIS_TEXT - 1, colour = "grey45"),
      panel.grid = element_blank(), panel.spacing.y = unit(0.8, "mm"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 2)
    )
}

ora_aligned <- function(ora_top5, modules, top_n = 5L, sig = H9C2_FDR_EXPLOR) {
  # Label placement scales with bar length: name inside (white) + q/p outside (dark)
  # when the name fits; q/p inside right-aligned (white) + name outside (dark) when it
  # doesn't; both outside (dark) when the bar is too short for either. bar_chars tunes
  # the fit threshold (free_x panels share physical width but not scale).
  bar_chars <- 78
  d <- ora_top5 |>
    filter(module %in% modules) |>
    group_by(module) |>
    slice_min(padj, n = top_n, with_ties = FALSE) |>
    mutate(k = row_number(), xmax = max(-log10(padj))) |>
    ungroup() |>
    mutate(
      module = factor(module, levels = modules),
      database = factor(database, levels = names(DB_COLORS)),
      lp = -log10(padj),
      yb = top_n + 1 - k,
      sig_fdr = padj < sig,
      star = ifelse(padj < sig, "✱ ", ifelse(p < sig, "· ", "")),
      name = clean_display_label(pathway),
      name = ifelse(nchar(name) > 34, paste0(substr(name, 1, 33), "…"), name),
      lab = paste0(star, "q", fmt_pq(padj), " p", fmt_pq(p)),
      gap = 0.03 * xmax,
      cap = lp / xmax * bar_chars,
      name_fits = nchar(name) <= cap,
      sig_inside = !name_fits & nchar(lab) <= cap,
      both_outside = !name_fits & nchar(lab) > cap,
      combo = paste0(name, "   ", lab)
    )
  ggplot(d, aes(lp, yb, fill = database)) +
    geom_col(width = 0.72, colour = "grey35", linewidth = 0.2, orientation = "y") +
    geom_col(
      data = filter(d, sig_fdr), width = 0.72, fill = NA,
      colour = "black", linewidth = 0.7, orientation = "y"
    ) +
    geom_text(
      data = filter(d, name_fits), aes(x = gap, label = name),
      hjust = 0, size = 1.5, fontface = "bold", colour = "white"
    ) +
    geom_text(
      data = filter(d, name_fits), aes(x = lp + gap, label = lab),
      hjust = 0, size = 1.4, fontface = "bold", colour = "grey15"
    ) +
    geom_text(
      data = filter(d, sig_inside), aes(x = lp - gap, label = lab),
      hjust = 1, size = 1.4, fontface = "bold", colour = "white"
    ) +
    geom_text(
      data = filter(d, sig_inside), aes(x = lp + gap, label = name),
      hjust = 0, size = 1.5, fontface = "bold", colour = "grey15"
    ) +
    geom_text(
      data = filter(d, both_outside), aes(x = lp + gap, label = combo),
      hjust = 0, size = 1.5, fontface = "bold", colour = "grey15"
    ) +
    facet_wrap(~module, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = DB_COLORS, drop = FALSE, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
    scale_y_continuous(limits = c(0.4, top_n + 0.6), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(
      title = sprintf("Top %d ORA pathways  (✱ FDR<%.2f  · p<%.2f)", top_n, sig, sig),
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.title = element_text(face = "bold", size = FIG_AXIS_TEXT + 1, colour = "grey15"),
      strip.text = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
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
