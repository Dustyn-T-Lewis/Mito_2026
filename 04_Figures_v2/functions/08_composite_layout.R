# 04_Figures_v2/functions/08_composite_layout.R
# Pipeline step 08: stitch refined standalone panels into composites. The last
# assembly step — sources nothing the panels do not already provide.
source(here::here("04_Figures_v2", "functions", "01_style_palettes_theme.R"))

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
  pdf_dev <- get_pdf_device()
  out_pdf <- file.path(base, "b_reports", "main", "pdf")
  out_png <- file.path(base, "b_reports", "main", "png")
  for (d in c(out_pdf, out_png)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(out_pdf, paste0(name, ".pdf"))
  pdf_ok <- tryCatch(
    {
      ggsave(pdf_path, plot,
        width = width_mm, height = height_mm, units = "mm",
        device = pdf_dev, limitsize = FALSE
      )
      file.exists(pdf_path) && file.info(pdf_path)$size > 1000
    },
    error = function(e) FALSE
  )

  if (!pdf_ok) {
    message("preferred pdf device failed or wrote empty file — falling back to base grDevices::pdf()")
    grDevices::pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4)
    tryCatch(print(plot), finally = grDevices::dev.off())
  }

  ggsave(file.path(out_png, paste0(name, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 300, limitsize = FALSE
  )
  invisible(file.path(out_png, paste0(name, ".png")))
}
