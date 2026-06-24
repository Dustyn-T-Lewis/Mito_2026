# 04_Figures/functions/08_composite_layout.R
# Pipeline step 08: stitch refined standalone panels into composites. The last
# assembly step — sources nothing the panels do not already provide.
source(here::here("04_Figures", "functions", "01_style_palettes_theme.R"))

library(patchwork)

# Panel letter is baked into the title (constant gap to the title text, so all
# panels read uniformly) rather than a free-floating tag whose offset tracks
# each panel's y-axis width.
add_tag <- function(p, tag) {
  cur <- p$labels$title
  p + labs(
    tag = NULL,
    title = paste0(tag, "  ", if (is.null(cur)) "" else cur)
  )
}

composite_caption <- function(text) {
  patchwork::plot_annotation(
    caption = text,
    theme = theme(plot.caption = element_text(
      size = 5, color = "grey35",
      hjust = 0, lineheight = 1.1
    ))
  )
}

save_composite <- function(plot, base, name, width_mm, height_mm) {
  out_png <- file.path(base, "b_reports", "main", "png")
  dir.create(out_png, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_png, paste0(name, ".png"))
  ggsave(png_path, plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 300, limitsize = FALSE
  )
  invisible(png_path)
}
