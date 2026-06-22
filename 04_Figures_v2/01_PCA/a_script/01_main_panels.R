#!/usr/bin/env Rscript
# F01 PCA — standalone sample PCA + PERMANOVA on the imputed matrix.
# Reads existing pipeline outputs only; never re-runs 01-03.

source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE <- here::here("04_Figures_v2", "01_PCA")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

res <- build_pca_panel()
ggsave(file.path(RPT_PDF, "MAIN_F01_pca.pdf"), res$plot,
  width = 110, height = 95,
  units = "mm", device = pdf_dev
)
ggsave(file.path(RPT_PNG, "MAIN_F01_pca.png"), res$plot,
  width = 110, height = 95,
  units = "mm", dpi = 300
)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  figure_title = "F01 — Sample PCA + PERMANOVA on the imputed protein matrix",
  sheet_specs = list(
    list(
      name = "pca_scores", df = res$scores,
      role = "Panel coordinates — the PCA scatter points",
      contents = "PC1/PC2 scores per sample (Col_ID) with Group; % variance in axis titles"
    ),
    list(
      name = "permanova", df = res$permanova,
      role = "Stats annotation block on the figure",
      contents = "adonis2 R2 and p for overall Group + Disease/Transplant/Rescue pairwise (3 uncorrected pairwise tests), plus betadisper dispersion p"
    )
  )
)
message("F01 PCA rebuilt")
