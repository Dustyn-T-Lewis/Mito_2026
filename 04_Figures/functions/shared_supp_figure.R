# Supplementary panels and the S Figures the manuscript cites. Each panel is saved as a
# PNG for browsing and as one labelled vector PDF page; once every page of an S Figure
# exists, the pages are combined into b_reports/S<n>_Figure.pdf of the figure that owns it.

# Page order within each S Figure, keyed by panel file stem.
S_FIGURE_PAGES <- list(
  S1 = c("SUPP_F01_pca_group", "SUPP_F01_pca_transplant", "SUPP_F01_fgsea_pc1", "SUPP_F01_ora", "SUPP_F01_compartments"),
  S2 = c("SUPP_F02_enrichment_by_contrast", "SUPP_F02_disease_specificity", "SUPP_F02_reversal"),
  S3 = c("SUPP_F03_construction", "SUPP_F03_preservation", "SUPP_F03_hub_map"),
  S4 = "SUPP_F03_orthogonal_axes"
)

s_page_label <- function(stem) {
  id <- names(S_FIGURE_PAGES)[vapply(S_FIGURE_PAGES, \(p) stem %in% p, logical(1))]
  if (length(id) != 1) {
    return(NA_character_)
  }
  pages <- S_FIGURE_PAGES[[id]]
  if (length(pages) == 1) {
    return(sprintf("%s Figure", id))
  }
  sprintf("%s Figure %s", id, LETTERS[match(stem, pages)])
}

# `...` goes to ggsave unchanged, so each panel keeps its own PNG settings.
save_supp_panel <- function(plot, dir, stem, width, height, ...) {
  ggplot2::ggsave(file.path(dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "mm", dpi = 300, ...
  )
  label <- s_page_label(stem)
  if (is.na(label)) {
    return(invisible())
  }
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
  build_s_figure(sub(" .*", "", label))
}

# Combines an S Figure's pages once all of them exist; the output lands beside the
# supp/ folder of the figure that owns the S Figure's first page.
build_s_figure <- function(id) {
  pages <- S_FIGURE_PAGES[[id]]
  owner <- list.dirs(here::here("04_Figures"), recursive = FALSE)
  owner <- owner[startsWith(basename(owner), sub("^SUPP_(F\\d+)_.*", "\\1", pages[1]))]
  supp <- file.path(owner, "b_reports", "supp")
  pdfs <- file.path(supp, paste0(pages, ".pdf"))
  if (!all(file.exists(pdfs))) {
    return(invisible())
  }
  out <- file.path(owner, "b_reports", sprintf("%s_Figure.pdf", id))
  qpdf::pdf_combine(pdfs, out)
  message("Wrote ", out, " (", length(pdfs), " pages)")
  invisible(out)
}
