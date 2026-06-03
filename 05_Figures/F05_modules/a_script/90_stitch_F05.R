#!/usr/bin/env Rscript
# F05_modules orchestrator. Headline = knowledge-driven set-score modules
# (01_set_scores → 02_main_panels); WGCNA is the exploratory de-novo supplement
# (build network first, then 91_wgcna_supp). Then optional Box copy.

source(here::here("05_Figures", "F05_modules", "a_script", "01_set_scores.R"))
source(here::here("05_Figures", "F05_modules", "a_script", "02_main_panels.R"))
source(here::here("05_Figures", "shared", "build_wgcna_network.R"))
source(here::here("05_Figures", "F05_modules", "a_script", "91_wgcna_supp.R"))

box <- Sys.getenv("H9C2_BOX_DIR", "")
if (nzchar(box) && dir.exists(box)) {
  RPT <- here::here("05_Figures", "F05_modules", "b_reports")
  for (sub in c("figures/pdf", "figures/png"))
    dir.create(file.path(box, "02_Figures", sub), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main", "pdf", "MAIN_F05_composite.pdf"),
            file.path(box, "02_Figures", "figures", "pdf", "Figure_5_setscore_modules.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main", "png", "MAIN_F05_composite.png"),
            file.path(box, "02_Figures", "figures", "png", "Figure_5_setscore_modules.png"), overwrite = TRUE)
  message("Copied F05 outputs to Box")
}
message("F05 complete")
