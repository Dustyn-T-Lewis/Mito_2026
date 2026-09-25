#!/usr/bin/env Rscript
# Packages the supplementary files under the names the manuscript cites. Each S Figure
# is one multipage PDF, one panel per page, stamped with its label; each S Table is the
# figure's supplementary workbook. Run after the figure scripts, from the project root.

pacman::p_load(here, magick, purrr)

fig <- \(...) here::here("04_Figures", ...)
supp_png <- \(dir, name) fig(dir, "b_reports", "supp", paste0("SUPP_", name, ".png"))

S_FIGURES <- list(
  S1 = supp_png("F01_Proteome_Overview", c(
    "F01_pca_group", "F01_pca_transplant", "F01_fgsea_pc1", "F01_ora", "F01_compartments"
  )),
  S2 = supp_png("F02_Enrich_Volcanoes", c(
    "F02_enrichment_by_contrast", "F02_disease_specificity", "F02_reversal"
  )),
  S3 = supp_png("F03_WGCNA", c("F03_construction", "F03_preservation", "F03_hub_map")),
  S4 = supp_png("F03_WGCNA", "F03_orthogonal_axes")
)

S_TABLES <- c(
  S1 = fig("F01_Proteome_Overview", "c_data", "F01_supplementary.xlsx"),
  S2 = fig("F02_Enrich_Volcanoes", "c_data", "F02_supplementary.xlsx"),
  S3 = fig("F03_WGCNA", "c_data", "F03_supplementary.xlsx")
)

out_dir <- here::here("MANUSCRIPT_FILES")
stopifnot(file.exists(unlist(S_FIGURES)), file.exists(S_TABLES))

stamp_page <- function(path, label) {
  img <- image_read(path)
  band <- round(image_info(img)$width * 0.03)
  img |>
    image_border("white", sprintf("0x%d", band)) |>
    image_crop(sprintf("%dx%d", image_info(img)$width, image_info(img)$height + band)) |>
    image_annotate(label,
      size = round(band * 0.6), weight = 700, gravity = "northwest",
      location = sprintf("+%d+%d", band, round(band * 0.2))
    )
}

iwalk(S_FIGURES, function(paths, id) {
  n <- length(paths)
  labels <- if (n == 1) {
    sprintf("%s Figure", id)
  } else {
    sprintf("%s Figure %s (page %d of %d)", id, LETTERS[seq_len(n)], seq_len(n), n)
  }
  pages <- image_join(map2(paths, labels, stamp_page))
  image_write(pages, file.path(out_dir, sprintf("%s_Figure.pdf", id)), format = "pdf", density = 300)
})

iwalk(S_TABLES, \(path, id) file.copy(path, file.path(out_dir, sprintf("%s_Table.xlsx", id)), overwrite = TRUE))

message("Wrote ", length(S_FIGURES), " S Figure PDFs and ", length(S_TABLES), " S Table workbooks to ", out_dir)
