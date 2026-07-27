# Horizontal NES bar panel shared by the F02 claim and aim figures: pathway name
# fitted inside the bar, database on the axis in its own colour, FDR at the tip.

source(here::here("04_Figures", "functions", "shared_theme_palettes.R"))

pacman::p_load(ggplot2, dplyr, tibble, stringr, ggfittext)

fmt_q <- function(padj) {
  ifelse(padj < 0.001, formatC(padj, format = "e", digits = 0),
    formatC(padj, format = "f", digits = 3)
  )
}

# Direction and significance in one four-level fill, so a panel carries a single
# key instead of a colour legend plus an opacity legend. The n.s. tints are mixed
# toward white rather than made transparent, keeping ggfittext's contrast check on
# a solid background.
.tint <- function(col, amount = 0.62) {
  grDevices::rgb(
    t((1 - amount) * grDevices::col2rgb(col) + amount * 255),
    maxColorValue = 255
  )
}

NES_FILL_KEY <- c(
  Up = unname(H9C2_PAL_DIR[["Up"]]),
  Down = unname(H9C2_PAL_DIR[["Down"]]),
  `Up (n.s.)` = .tint(H9C2_PAL_DIR[["Up"]]),
  `Down (n.s.)` = .tint(H9C2_PAL_DIR[["Down"]])
)

# How many half-axes a panel occupies: one if its bars all run the same way, two
# if it has both. Pass these as patchwork widths so a two-sided panel gets twice
# the horizontal room of a one-sided one.
nes_sides <- function(d) sum(c(any(d$NES < 0), any(d$NES > 0)))

# geom_fit_text reflows and shrinks each name to sit inside its own bar and picks
# black or white from the fill, so no calibration is needed. The x-limit runs from
# the panel's own bar extremes with a pad reserving room for the FDR at each tip.
nes_bar_panel <- function(d, title, subtitle, xlab = "NES", base_theme = FIG_THEME,
                          name_size = 7, q_size = 2.4, show_names = TRUE) {
  if (!"display" %in% names(d)) d$display <- as.character(d$label)

  levs <- levels(d$label)
  db_of <- as.character(d$database)[match(levs, as.character(d$label))]
  axis_cols <- unname(DB_COLORS[db_of])

  d <- mutate(d,
    fill_key = factor(
      paste0(direction, if_else(sig == "TRUE", "", " (n.s.)")),
      levels = names(NES_FILL_KEY)
    ),
    q_x = NES + sign(NES) * 0.04,
    q_hjust = if_else(NES > 0, 0, 1),
    q_face = if_else(sig == "TRUE", "bold", "plain")
  )

  q_pad <- 0.13 * max(nchar(fmt_q(d$padj))) + 0.1
  lim <- c(
    if (any(d$NES < 0)) min(d$NES) - q_pad else 0,
    if (any(d$NES > 0)) max(d$NES) + q_pad else 0
  )

  # A level with no rows draws an empty key glyph, so zero-width bars carrying all
  # four levels give each key entry data to draw from and every panel shows the
  # same four-part key.
  key_rows <- tibble::tibble(
    NES = 0, label = d$label[1],
    fill_key = factor(names(NES_FILL_KEY), levels = names(NES_FILL_KEY))
  )

  ggplot(d, aes(NES, label)) +
    geom_col(aes(fill = fill_key), width = 0.74, colour = "grey20", linewidth = 0.25) +
    geom_col(
      data = key_rows, aes(fill = fill_key), width = 0,
      colour = "grey20", linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    (if (show_names) {
      ggfittext::geom_fit_text(
        aes(
          xmin = pmin(NES, 0), xmax = pmax(NES, 0),
          ymin = as.numeric(label) - 0.37, ymax = as.numeric(label) + 0.37,
          label = display, fill = fill_key
        ),
        reflow = TRUE, grow = FALSE, contrast = TRUE,
        fontface = "bold", size = name_size, min.size = 0,
        padding.x = grid::unit(0.6, "mm"), padding.y = grid::unit(0.4, "mm"),
        show.legend = FALSE
      )
    }) +
    geom_text(
      aes(x = q_x, label = fmt_q(padj), hjust = q_hjust, fontface = q_face),
      colour = "black", size = q_size, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = NES_FILL_KEY, name = NULL, drop = FALSE,
      guide = guide_legend(nrow = 1, keywidth = unit(3.2, "mm"), keyheight = unit(3.2, "mm"))
    ) +
    scale_y_discrete(labels = stats::setNames(db_of, levs)) +
    scale_x_continuous(limits = lim, expand = expansion(mult = 0.02)) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    base_theme +
    theme(
      axis.text.y = element_text(colour = axis_cols, face = "bold"),
      panel.grid.major.y = element_blank()
    )
}
