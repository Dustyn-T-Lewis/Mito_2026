# Panels for the F03 module card: each module is one aligned row of its protein-count bar,
# member-response heatmap, group-mean eigengene trajectory, and top ORA pathways. The four
# columns facet on the same module order, so a row reads left to right as one module.
pacman::p_load(ggplot2, dplyr, forcats, stringr, patchwork, ggfittext, ggtext, ggnewscale)

# Lightens a module colour toward white, so one hue carries emphasis as depth:
# full strength for a featured module, washed out for the rest.
tint_colour <- function(col, amount = 0.62) {
  grDevices::rgb(t((1 - amount) * grDevices::col2rgb(col) + amount * 255), maxColorValue = 255)
}

fmt_fdr <- function(q) {
  ifelse(is.na(q), "", ifelse(q < 0.01, "<.01", sprintf("%.2f", q)))
}

# White text on dark bars, black on light, by perceived luminance — a saturated red or blue
# bar reads white, a pale tan or greenyellow bar reads black.
text_on_fill <- function(hex) {
  rgb <- grDevices::col2rgb(hex)
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(lum < 140, "white", "black")
}

# Card text sizes and line weights. The figure prints large (one card row per featured
# module), so these absolute sizes read large; point sizes for element_text, mm sizes for
# geom_text/geom_col.
S_TITLE <- 13
S_SUB <- 9
S_AXIS <- 8
S_LEG <- 7
# One shared size for every inline stat label (heatmap tiles and trajectory brackets).
S_STAT <- 2.5
S_TILE <- S_STAT
S_BRK <- S_STAT
S_HDR <- 3.6
S_ORANAME <- 7.5
S_ORAQ <- 3.3
LW_LINE <- 1.6
LW_PT <- 1.8
LW_EB <- 0.65
LW_BRK <- 0.75
LW_TILE <- 0.3

# Shared vertical extent so the count bars and the heatmap tiles land at the same height.
CARD_Y <- c(0.5, 1.85)

# Panel title/subtitle tier: a clear step below the figure title, above the axis text.
card_title_theme <- function() {
  theme(
    plot.title = element_text(face = "bold", size = S_TITLE, colour = "grey10"),
    plot.subtitle = element_text(face = "italic", size = S_SUB, colour = "grey35")
  )
}

# Drop the panel box, keep only the L of axis lines.
clean_axes <- function(y = TRUE) {
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "grey40", linewidth = 0.4),
    axis.line.y = if (y) element_line(colour = "grey40", linewidth = 0.4) else element_blank()
  )
}

# Protein-count bar, one facet per module aligned to the rows beside it. Bars run right to
# left from a shared baseline; a "module · count" header names each bar; the y-axis sits on
# the right, next to the heatmap.
panel_module_count <- function(mod_size, modules, emphasize = modules) {
  d <- mod_size |>
    filter(module %in% modules) |>
    mutate(
      module = factor(module, levels = modules),
      fill_col = ifelse(module %in% emphasize, as.character(module), vapply(as.character(module), tint_colour, character(1))),
      name_col = ifelse(module %in% emphasize, "grey20", "grey55")
    )
  ggplot(d, aes(n, 1)) +
    geom_col(aes(fill = I(fill_col)),
      width = 0.5, colour = "grey30", linewidth = 0.35, orientation = "y"
    ) +
    ggtext::geom_richtext(
      aes(x = 0, y = 1.54, label = sprintf(
        "<b style='color:%s'>%s</b>  ·  <b style='color:#333333'>%d</b>", name_col, module, n
      )),
      hjust = 1, vjust = 0, size = S_HDR, fill = NA, label.colour = NA,
      label.padding = unit(c(0, 0, 0, 0), "pt")
    ) +
    facet_wrap(~module, ncol = 1) +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_continuous(limits = CARD_Y, expand = c(0, 0), position = "right") +
    coord_cartesian(clip = "off") +
    labs(title = "Proteins", subtitle = "count per module", x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      strip.text = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = S_AXIS - 2, colour = "grey35", face = "bold"),
      panel.grid = element_blank(), panel.spacing.y = unit(2, "mm"),
      legend.position = "none", plot.margin = margin(2, 5, 2, 11)
    ) +
    clean_axes(y = TRUE) +
    card_title_theme()
}

# Member-response heatmap: two contiguous rows of tiles per module across the primary contrasts
# plus the interaction. Top row is fGSEA NES on a green–purple scale (competitive, descriptive:
# it lights up almost everywhere, showing why the competitive null fails); the tile prints the
# NES only. Bottom row is fry on a red–blue scale (self-contained gate), boxed at FDR < 0.05
# with q printed on the significant tiles alone, so the gate reads as colour, not text. Gradient
# keys sit below the column.
panel_module_fry <- function(settests, module_nes, modules, contrasts) {
  ctr_lab <- c(
    Disease = "Disease", Transplant = "Transplant", Rescue = "Rescue",
    Recovery = "Recovery", Interaction = "Interaction"
  )
  up <- unname(H9C2_PAL_DIR[["Up"]])
  down <- unname(H9C2_PAL_DIR[["Down"]])
  lv <- function(m) factor(m, levels = modules)
  cv <- function(c) factor(c, levels = contrasts)
  nes_d <- module_nes |>
    filter(module %in% modules, contrast %in% contrasts) |>
    transmute(
      module = lv(module), contrast = cv(contrast), layer = "NES",
      NES,
      txt = if_else(abs(NES) > 3.5, "white", "black"),
      lab = sprintf("%.1f", NES)
    )
  fry_d <- settests |>
    filter(module %in% modules, contrast %in% contrasts) |>
    transmute(
      module = lv(module), contrast = cv(contrast), layer = "fry",
      sig = !is.na(fry_fdr) & fry_fdr < H9C2_FDR_STD,
      signed = if_else(direction == "Up", 1, -1) * -log10(pmax(fry_fdr, 1e-4)),
      txt = if_else(abs(signed) > 0.9, "white", "black"),
      face = if_else(sig, "bold", "plain"),
      lab = paste0("p=", fmt_fdr(fry_p), "\nq=", fmt_fdr(fry_fdr))
    )
  bar <- guide_colourbar(barwidth = grid::unit(24, "mm"), barheight = grid::unit(2.8, "mm"), title.position = "top", title.hjust = 0.5)
  ggplot(mapping = aes(contrast, layer)) +
    geom_tile(data = nes_d, aes(fill = NES), width = 1, height = 1, colour = "grey25", linewidth = LW_TILE) +
    geom_text(data = nes_d, aes(label = lab, colour = I(txt)), fontface = "plain", size = S_TILE) +
    scale_fill_gradientn(
      colours = c("#1B7837", "white", "#762A83"), limits = c(-5, 5), oob = scales::squish,
      breaks = c(-5, 0, 5), name = "fGSEA NES (up purple)", guide = bar
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = fry_d, aes(fill = signed), width = 1, height = 1, colour = "grey25", linewidth = LW_TILE) +
    geom_tile(data = filter(fry_d, sig), aes(fill = signed), colour = "black", linewidth = 1.2, width = 1, height = 1) +
    geom_text(data = fry_d, aes(label = lab, colour = I(txt), fontface = face), size = S_TILE, lineheight = 0.85) +
    scale_fill_gradient2(
      low = down, mid = "white", high = up, midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
      breaks = c(-1.5, 0, 1.5), name = "fry −log₁₀ q (up red)", guide = bar
    ) +
    facet_wrap(~module, ncol = 1) +
    scale_x_discrete(labels = ctr_lab, expand = c(0, 0)) +
    scale_y_discrete(limits = c("fry", "NES"), expand = expansion(add = 0.18)) +
    coord_cartesian(clip = "off") +
    labs(title = "Member response", subtitle = "NES top · fry gate below", x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      strip.text = element_blank(),
      axis.text.y = element_text(size = S_AXIS - 1, colour = "grey25", hjust = 1, face = "bold"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = S_AXIS - 1, colour = "grey15"),
      panel.grid = element_blank(), panel.spacing.y = unit(2, "mm"),
      legend.position = "bottom", legend.box = "horizontal",
      legend.title = element_text(size = S_LEG, face = "bold"), legend.text = element_text(size = S_LEG - 1),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    clean_axes(y = FALSE) +
    card_title_theme()
}

# Which designed contrast spans which group pair (groups CTL=1, MitoTx=2, Phe=3, Phe+Mito=4).
POSTHOC_SPANS <- tibble::tribble(
  ~contrast, ~x1, ~x2,
  "Transplant", 1, 2,
  "Rescue", 3, 4,
  "Disease", 1, 3,
  "Recovery", 1, 4
)

# Greedy interval stacking: narrowest brackets first, each to the lowest level with no x-overlap.
# Keeps a lone bracket (e.g. pink's Recovery) at level 1 instead of floating.
posthoc_levels <- function(df) {
  df <- df[order(df$x2 - df$x1, df$x1), , drop = FALSE]
  placed <- list()
  df$level <- NA_integer_
  for (i in seq_len(nrow(df))) {
    lev <- 1L
    repeat {
      clash <- any(vapply(placed, function(pl) {
        pl$level == lev && df$x1[i] < pl$x2 && df$x2[i] > pl$x1
      }, logical(1)))
      if (!clash) break
      lev <- lev + 1L
    }
    df$level[i] <- lev
    placed[[length(placed) + 1L]] <- list(x1 = df$x1[i], x2 = df$x2[i], level = lev)
  }
  df
}

# Group-mean module eigengene (mean ± SE), one facet per module, all lines the same weight.
# Post-hoc brackets sit over the specific group pairs whose eigengene contrast is nominally
# significant (limma on the eigengene), starred at FDR < 0.05.
panel_module_trajectory <- function(group_eig, modules, mod_stats = NULL, emphasize = modules) {
  line_cols <- stats::setNames(
    ifelse(modules %in% emphasize, modules, vapply(modules, tint_colour, character(1))),
    modules
  )
  d <- group_eig |>
    filter(module %in% modules) |>
    mutate(
      module = factor(module, levels = modules),
      gx = match(as.character(Group), H9C2_GROUP_LEVELS),
      col = line_cols[as.character(module)]
    )
  p <- ggplot(d, aes(gx, mean_eig, colour = col)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.35) +
    geom_errorbar(aes(ymin = mean_eig - se, ymax = mean_eig + se), width = 0.18, linewidth = LW_EB) +
    geom_line(linewidth = LW_LINE) +
    geom_point(size = LW_PT) +
    facet_wrap(~module, ncol = 1, scales = "free_y") +
    scale_colour_identity() +
    scale_x_continuous(breaks = 1:4, labels = unname(GROUP_LABELS), expand = expansion(mult = 0.06)) +
    scale_y_continuous(breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.08))) +
    labs(
      title = "Module eigengene trajectory",
      subtitle = "group mean ± SE (n = 6)",
      x = NULL, y = "eigengene (a.u.)"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = S_AXIS),
      axis.text.y = element_text(size = S_AXIS - 1, colour = "grey40", face = "bold"),
      axis.title.y = element_text(size = S_AXIS, face = "bold"),
      panel.grid = element_blank(), panel.spacing.y = unit(2, "mm"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 2)
    ) +
    clean_axes() +
    card_title_theme()
  if (!is.null(mod_stats)) {
    span <- d |>
      group_by(module) |>
      summarise(top = max(mean_eig + se), rng = pmax(max(mean_eig + se) - min(mean_eig - se), 1e-6), .groups = "drop")
    brk0 <- mod_stats |>
      filter(module %in% modules, contrast %in% POSTHOC_SPANS$contrast, p < 0.05) |>
      inner_join(POSTHOC_SPANS, by = "contrast")
    brk <- if (nrow(brk0)) {
      brk0 |>
        group_split(module) |>
        lapply(posthoc_levels) |>
        bind_rows() |>
        inner_join(span, by = "module") |>
        mutate(
          module = factor(module, levels = modules),
          yb = top + (0.16 + 0.42 * (level - 1)) * rng,
          ytick = yb - 0.04 * rng,
          ylab = yb + 0.015 * rng,
          xmid = (x1 + x2) / 2,
          col = "grey15",
          p_lab = if_else(p < 0.05, sprintf("<b>p=%.3f</b>", p), sprintf("p=%.3f", p)),
          q_lab = if_else(fdr < H9C2_FDR_STD, sprintf("<b>q=%.2f ✱</b>", fdr), sprintf("q=%.2f", fdr)),
          lab = sprintf("<b>%s</b><br>%s  %s", contrast, p_lab, q_lab)
        )
    } else {
      brk0
    }
    if (nrow(brk)) {
      ceil <- brk |>
        group_by(module) |>
        summarise(y = max(yb) + 0.34 * rng[1], .groups = "drop")
      p <- p +
        geom_blank(data = ceil, aes(x = 1, y = y), inherit.aes = FALSE) +
        geom_segment(data = brk, aes(x = x1, xend = x2, y = yb, yend = yb, colour = col), inherit.aes = FALSE, linewidth = LW_BRK) +
        geom_segment(data = brk, aes(x = x1, xend = x1, y = yb, yend = ytick, colour = col), inherit.aes = FALSE, linewidth = LW_BRK) +
        geom_segment(data = brk, aes(x = x2, xend = x2, y = yb, yend = ytick, colour = col), inherit.aes = FALSE, linewidth = LW_BRK) +
        ggtext::geom_richtext(
          data = brk, aes(x = xmid, y = ylab, label = lab, colour = I(col)), inherit.aes = FALSE,
          vjust = 0, size = S_BRK, lineheight = 1.0, fill = NA, label.colour = NA,
          label.padding = unit(c(0, 0, 0, 0), "pt")
        )
    }
  }
  p
}

# Top ORA pathways per module. The name sits inside its bar when the bar has room; a near-null
# hit (FDR > 0.75, a sliver) shows its FDR then its name to the right of the bar instead of
# being squished. Database axis labels carry their own colour; the FDR is bold past the tip.
# Show the top 3 pathways, or the top 5 when all 5 clear FDR — a well-annotated module earns
# the extra rows, a thin one stays compact.
ora_aligned <- function(ora_top5, modules) {
  d <- ora_top5 |>
    filter(module %in% modules) |>
    group_by(module) |>
    filter(row_number(padj) <= if_else(sum(padj < 0.05) >= 5L, 5L, 3L)) |>
    ungroup() |>
    mutate(
      module = factor(module, levels = modules),
      database = factor(database, levels = names(ORA_DB_COLORS)),
      lp = -log10(padj),
      sig_fdr = padj < 0.05,
      fill_col = if_else(sig_fdr, as.character(module), tint_colour(as.character(module))),
      name_col = text_on_fill(fill_col),
      border_col = if_else(sig_fdr, "grey5", "grey55"),
      name = clean_display_label(pathway),
      qlab = ifelse(padj < 0.01, "< 0.01", sprintf("%.2f", padj)),
      outside = padj > 0.75
    ) |>
    arrange(module, lp) |>
    mutate(row = factor(paste(module, name, sep = ""), levels = unique(paste(module, name, sep = ""))))

  db_of <- stats::setNames(as.character(d$database), as.character(d$row))

  ggplot(d, aes(lp, row, fill = fill_col)) +
    geom_col(aes(colour = border_col), width = 0.88, linewidth = 0.5, orientation = "y") +
    scale_colour_identity() +
    ggfittext::geom_fit_text(
      data = ~ filter(.x, !outside), aes(xmin = 0, xmax = lp, label = name, colour = name_col),
      reflow = TRUE, grow = FALSE, fontface = "bold",
      size = S_ORANAME, min.size = 5,
      padding.x = grid::unit(0.5, "mm"), padding.y = grid::unit(0.3, "mm"),
      show.legend = FALSE
    ) +
    geom_text(
      data = ~ filter(.x, outside), aes(x = lp, label = sprintf("%s   %s", qlab, str_trunc(name, 44))),
      hjust = -0.15, size = 2.5, colour = "grey25"
    ) +
    geom_text(
      data = ~ filter(.x, !outside),
      aes(x = lp, label = qlab, fontface = if_else(sig_fdr, "bold", "plain")),
      hjust = -0.18, size = S_ORAQ, colour = "grey10"
    ) +
    facet_wrap(~module, ncol = 1, scales = "free") +
    scale_fill_identity() +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15)),
      breaks = function(lims) Filter(function(b) b >= 0, scales::breaks_pretty(3)(lims))
    ) +
    scale_y_discrete(
      labels = \(x) {
        db <- unname(db_of[x])
        sprintf("<span style='color:%s'>%s</span>", DB_COLORS[db], db)
      },
      expand = expansion(add = 0.5)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Top ORA pathways",
      subtitle = "BH-FDR; bold q < 0.05",
      x = expression(-log[10] ~ FDR), y = NULL
    ) +
    FIG_THEME +
    theme(
      strip.text = element_blank(),
      axis.text.y = ggtext::element_markdown(size = S_AXIS - 1, face = "bold", hjust = 1),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = S_AXIS - 2, colour = "grey35", face = "bold"),
      axis.title.x = element_text(size = S_AXIS - 1, colour = "grey35"),
      panel.grid = element_blank(), panel.spacing.y = unit(2, "mm"),
      legend.position = "none", plot.margin = margin(2, 3, 2, 1)
    ) +
    clean_axes() +
    card_title_theme()
}

# One card row per module: count bar, fry heatmap, eigengene trajectory, ORA pathways, all
# faceted on the same module order so the row aligns. Returns a bare patchwork; callers add
# the title (main vs supplement).
panel_module_card <- function(group_eig, mod_size, ora_top5, settests, module_nes, fry_contrasts,
                              mod_stats = NULL, emphasize = NULL) {
  modules <- arrange(filter(mod_size, module %in% unique(group_eig$module)), desc(n))$module
  if (is.null(emphasize)) emphasize <- modules
  wrap_plots(
    list(
      panel_module_count(mod_size, modules, emphasize),
      panel_module_fry(settests, module_nes, modules, fry_contrasts),
      panel_module_trajectory(group_eig, modules, mod_stats, emphasize),
      ora_aligned(ora_top5, modules)
    ),
    nrow = 1, widths = c(0.5, 1.0, 1.12, 2.2)
  )
}
