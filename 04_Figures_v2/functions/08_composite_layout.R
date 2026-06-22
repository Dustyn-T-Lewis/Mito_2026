# 04_Figures_v2/functions/08_composite_layout.R
# Pipeline step 08: stitch refined standalone panels into composites. The last
# assembly step — sources nothing the panels do not already provide.
source(here::here("04_Figures_v2", "functions", "01_style_palettes_theme.R"))

library(patchwork)

add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(
      plot.tag = element_text(face = "bold", size = BASE_TAG),
      plot.tag.position = c(0.01, 0.99)
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
  pdf_dev <- get_pdf_device()
  out_pdf <- file.path(base, "b_reports", "main", "pdf")
  out_png <- file.path(base, "b_reports", "main", "png")
  for (d in c(out_pdf, out_png)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_pdf, paste0(name, ".pdf")), plot,
    width = width_mm,
    height = height_mm, units = "mm", device = pdf_dev, limitsize = FALSE
  )
  ggsave(file.path(out_png, paste0(name, ".png")), plot,
    width = width_mm,
    height = height_mm, units = "mm", dpi = 300, limitsize = FALSE
  )
  invisible(file.path(out_png, paste0(name, ".png")))
}
