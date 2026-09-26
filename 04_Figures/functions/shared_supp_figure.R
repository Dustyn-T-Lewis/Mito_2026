# Supplementary panels and the S Figures the manuscript cites. Each panel is saved as a
# PNG and, if it belongs to an S Figure, as one labelled vector PDF page. Once every page
# of that S Figure exists, the pages are combined into the figure's b_reports/.

# Page order within each S Figure, keyed by panel file stem.
S_FIGURE_PAGES <- list(
  S1 = c("SUPP_F01_pca_group", "SUPP_F01_pca_transplant", "SUPP_F01_fgsea_pc1", "SUPP_F01_ora", "SUPP_F01_compartments"),
  S2 = c("SUPP_F02_enrichment_by_contrast", "SUPP_F02_disease_specificity", "SUPP_F02_reversal"),
  S3 = c("SUPP_F03_construction", "SUPP_F03_preservation", "SUPP_F03_hub_map"),
  S4 = "SUPP_F03_orthogonal_axes"
)

# `...` goes to ggsave unchanged, so each panel keeps its own PNG settings.
save_supp_panel <- function(plot, dir, stem, width, height, ...) {
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "mm", dpi = 300, ...
  )
  id <- names(Filter(\(pages) stem %in% pages, S_FIGURE_PAGES))
  if (length(id) != 1) {
    return(invisible())
  }
  pages <- S_FIGURE_PAGES[[id]]
  label <- if (length(pages) == 1) paste(id, "Figure") else paste(id, "Figure", LETTERS[match(stem, pages)])

  band <- 8
  grDevices::cairo_pdf(file.path(dir, paste0(stem, ".pdf")),
    width = width / 25.4, height = (height + band) / 25.4
  )
  grid::grid.newpage()
  layout <- grid::grid.layout(2, 1, heights = grid::unit(c(band, height), "mm"))
  grid::pushViewport(grid::viewport(layout = layout))
  grid::grid.text(label,
    x = grid::unit(3, "mm"), just = "left",
    gp = grid::gpar(fontface = "bold", fontsize = 11, fontfamily = "Helvetica"),
    vp = grid::viewport(layout.pos.row = 1)
  )
  print(plot, vp = grid::viewport(layout.pos.row = 2))
  grDevices::dev.off()

  pdfs <- file.path(dir, paste0(pages, ".pdf"))
  if (all(file.exists(pdfs))) {
    out <- file.path(dirname(dir), paste0(id, "_Figure.pdf"))
    qpdf::pdf_combine(pdfs, out)
    message("Wrote ", out, " (", length(pdfs), " pages)")
  }
  invisible()
}
