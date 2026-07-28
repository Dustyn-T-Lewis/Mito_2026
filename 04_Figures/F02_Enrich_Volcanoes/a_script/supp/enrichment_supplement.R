#!/usr/bin/env Rscript
# Consolidated F02 enrichment supplement in one style: fgsea NES with the pathway
# name fitted inside the bar, database on the axis in its colour, FDR at the tip,
# one key at the bottom. Three figures, three questions:
#   enrichment_by_contrast : top themed sets in Disease / Transplant / Rescue
#   disease_specificity    : the three interaction reads (main beside interaction)
#   reversal               : one fixed set list read across the three contrasts
# The earlier aim / claim / themed panels are kept under supp/archive.

pacman::p_load(here, dplyr, tidyr, stringr, forcats, purrr, readr, ggplot2, patchwork)
fns <- here::here("04_Figures", "functions")
source(file.path(fns, "shared_data_loaders.R"))
source(file.path(fns, "shared_enrichment_ora.R"))
source(file.path(fns, "shared_nes_bars.R"))
source(file.path(fns, "shared_workbook.R"))

BASE <- here::here("04_Figures", "F02_Enrich_Volcanoes")
PNG <- file.path(BASE, "b_reports", "supp")
DAT_SUPP <- file.path(BASE, "c_data", "supp")
for (d in c(PNG, DAT_SUPP)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

FDR <- H9C2_FDR_STD
TOP_N <- 8L

# Pre-specified biological themes; a set is eligible if its name matches one. The
# regexes are the union the archived aim/claim/themed panels used, kept in one place.
THEME_RE <- c(
  Respiration = paste0(
    "OXIDATIVE_PHOSPHOR|OXPHOS|AEROBIC_RESPIRATION|RESPIRATORY_ELECTRON|CRISTAE|",
    "TCA|CITRATE_CYCLE|MITOCHONDRIAL_RIBOSOME|SLC25A|FATTY_ACID_METABOLISM|",
    "FATTY_ACID_OXIDATION|MITOCHONDRIAL_PROTEIN_IMPORT"
  ),
  `ECM / fibrosis` = paste0(
    "COLLAGEN|FIBRONECTIN|ECM|EXTRACELLULAR_MATRIX|PROTEOGLYCAN|SYNDECAN|",
    "ELASTIC_FIBRE|EPITHELIAL_MESENCHYMAL|FOCAL_ADHESION|INTEGRIN"
  ),
  Proliferation = paste0(
    "MYC_TARGETS|E2F_TARGETS|CELL_CYCLE|TELOMERE_ORGANIZATION|G2M_CHECKPOINT|",
    "MITOTIC_SPINDLE|MITOTIC_METAPHASE|S_PHASE|DNA_REPLICATION"
  ),
  Hypertrophic = "TGF_BETA|MYOGENESIS|DILATED_CARDIOMYOPATHY|STRIATED_MUSCLE",
  Endocytic = "ENDOCYT|CLATHRIN|BINDING_AND_UPTAKE_OF_LIGANDS|MACROPINOCYT|PINOCYT|PHAGOSOM"
)
# Endothelial shear-stress set; matches ECM only via "INTEGRINS" and is not fibrotic here.
OFF_THEME <- paste0(
  "REACTOME_TURBULENT_OSCILLATORY_DISTURBED_FLOW_SHEAR_STRESS_ACTIVATES_",
  "SIGNALING_BY_PIEZO1_AND_INTEGRINS_IN_ENDOTHELIAL_CELLS"
)

theme_of <- function(pathway) {
  hit <- vapply(THEME_RE, \(rx) str_detect(pathway, regex(rx, ignore_case = TRUE)), logical(1))
  if (any(hit)) names(THEME_RE)[which(hit)[1]] else NA_character_
}

fgsea_themed <- load_fgsea_cache() |>
  filter(
    database %in% CANONICAL_DBS, !pathway %in% MITO_DROP_SETS,
    !pathway %in% OFF_THEME, !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE)
  ) |>
  mutate(theme = vapply(pathway, theme_of, character(1))) |>
  filter(!is.na(theme))

set_pool <- load_set_pool(CANONICAL_DBS)

# Add the columns nes_bar_panel needs (display name, direction, significance),
# disambiguating two databases that clean to the same label.
bar_cols <- function(d) {
  d |>
    mutate(
      display = str_trunc(clean_display_label(pathway), 42),
      direction = factor(if_else(NES > 0, "Up", "Down"), levels = c("Up", "Down")),
      sig = factor(padj < FDR, levels = c(TRUE, FALSE))
    ) |>
    mutate(
      display = if (n_distinct(pathway) > 1) paste0(display, " (", database, ")") else display,
      .by = display
    )
}

sig_note <- function(d) sprintf("%d of %d sets FDR < %.2f", sum(d$sig == "TRUE", na.rm = TRUE), nrow(d), FDR)

FIG_SUB <- theme(
  plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 1),
  plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30")
)
bottom_key <- theme(legend.position = "bottom", legend.justification = "center")

## Figure 1 — enrichment by contrast -----------------------------------------
CONTRASTS <- c(Disease = "CTLvPHE", Transplant = "CTLvMITO", Rescue = "PHEvPHE_MITO")

by_contrast <- imap(CONTRASTS, function(ctr_id, ctr_name) {
  deduplicate_enrichment(filter(fgsea_themed, contrast == ctr_id), set_pool, cross_db = TRUE) |>
    slice_min(padj, n = TOP_N, with_ties = FALSE) |>
    bar_cols() |>
    arrange(theme, NES) |>
    mutate(label = factor(display, levels = unique(display)))
})
contrast_panels <- imap(by_contrast, \(d, nm) nes_bar_panel(d, nm, sig_note(d), q_size = 2.2))
contrast_panels[-1] <- lapply(contrast_panels[-1], \(p) p + guides(fill = "none"))

fig_contrast <- wrap_plots(contrast_panels, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "Pathway enrichment by contrast",
    subtitle = paste0(
      "Top ", TOP_N, " themed gene sets per contrast — Disease (Phe vs CTL), ",
      "Transplant (MitoTx vs CTL), Rescue (Phe+MitoTx vs Phe). fgsea NES, per-database ",
      "BH-FDR, EnrichmentMap-deduplicated. Sets trace to the fgsea supplementary table."
    ),
    theme = FIG_SUB
  ) &
  bottom_key
ggsave(file.path(PNG, "SUPP_F02_enrichment_by_contrast.png"), fig_contrast,
  width = 250, height = 120, units = "mm", dpi = 300, bg = "white"
)

## Figure 2 — disease specificity (interaction) ------------------------------
READS <- tibble::tribble(
  ~title, ~regex, ~main_id, ~claim, ~show_key,
  "ECM under rescue", THEME_RE[["ECM / fibrosis"]], "PHEvPHE_MITO",
  "MitoTx suppresses the established fibrotic program", TRUE,
  "ECM under transplant", THEME_RE[["ECM / fibrosis"]], "CTLvMITO",
  "MitoTx lowers ECM sets in healthy cells", FALSE,
  "Respiration under rescue", THEME_RE[["Respiration"]], "PHEvPHE_MITO",
  "The metabolic rescue is phenylephrine-specific", FALSE
)

read_row <- function(title, regex, main_id, claim, show_key) {
  main <- fgsea_themed |>
    filter(contrast == main_id, str_detect(pathway, regex(regex, ignore_case = TRUE))) |>
    slice_min(padj, n = TOP_N, with_ties = FALSE) |>
    bar_cols() |>
    mutate(label = fct_reorder(display, NES))
  disp_map <- distinct(main, pathway, display)
  inter <- fgsea_themed |>
    filter(contrast == "Interaction", pathway %in% main$pathway) |>
    left_join(disp_map, by = "pathway") |>
    mutate(
      direction = factor(if_else(NES > 0, "Up", "Down"), levels = c("Up", "Down")),
      sig = factor(padj < FDR, levels = c(TRUE, FALSE)),
      label = factor(display, levels = levels(main$label))
    )
  n_int <- sum(inter$sig == "TRUE", na.rm = TRUE)
  int_read <- if (n_int > 0) "disease-contingent" else "none detected; underpowered, not evidence of none"
  main_p <- nes_bar_panel(main, title, claim, xlab = "NES", q_size = 2.1)
  if (!show_key) main_p <- main_p + guides(fill = "none")
  inter_p <- nes_bar_panel(
    inter, "Phe-specific component",
    sprintf("(Phe+MitoTx - Phe) - (MitoTx - CTL) | %d of %d FDR < %.2f: %s", n_int, nrow(inter), FDR, int_read),
    xlab = "Interaction NES", q_size = 2.1, show_names = FALSE
  ) +
    guides(fill = "none")
  main_p + inter_p + plot_layout(widths = c(1.6, 1))
}

fig_specificity <- wrap_plots(pmap(READS, read_row), ncol = 1, guides = "collect") +
  plot_annotation(
    title = "Disease-specificity of the transplant response",
    subtitle = paste0(
      "Each read shows its main contrast beside the interaction (Phe+MitoTx - Phe) - ",
      "(MitoTx - CTL), which isolates the disease-contingent component. fgsea NES, ",
      "per-database BH-FDR; sets from the fgsea supplementary table, not deduplicated so no named set is hidden."
    ),
    theme = FIG_SUB
  ) &
  bottom_key
ggsave(file.path(PNG, "SUPP_F02_disease_specificity.png"), fig_specificity,
  width = 210, height = 210, units = "mm", dpi = 300, bg = "white"
)

## Figure 3 — reversal -------------------------------------------------------
REV_CONTRASTS <- c(Disease = "CTLvPHE", Rescue = "PHEvPHE_MITO", Transplant = "CTLvMITO")

top_sets <- function(regex, n) {
  fgsea_themed |>
    filter(contrast == "CTLvPHE", str_detect(pathway, regex(regex, ignore_case = TRUE))) |>
    slice_max(abs(NES), n = n, with_ties = FALSE) |>
    pull(pathway)
}
rev_sets <- c(
  top_sets(THEME_RE[["ECM / fibrosis"]], 4),
  top_sets(THEME_RE[["Proliferation"]], 4)
)
# One display + row order for all three columns, fixed on the Disease ranking so a
# reader tracks each set left to right.
rev_disp <- fgsea_themed |>
  filter(contrast == "CTLvPHE", pathway %in% rev_sets) |>
  bar_cols() |>
  arrange(NES) |>
  distinct(pathway, display)
rev_levels <- rev_disp$display

rev_panel <- function(ctr_id, ctr_name, show_names) {
  d <- fgsea_themed |>
    filter(contrast == ctr_id, pathway %in% rev_sets) |>
    left_join(rev_disp, by = "pathway") |>
    mutate(
      direction = factor(if_else(NES > 0, "Up", "Down"), levels = c("Up", "Down")),
      sig = factor(padj < FDR, levels = c(TRUE, FALSE)),
      label = factor(display, levels = rev_levels)
    )
  nes_bar_panel(d, ctr_name, sig_note(d), xlab = "NES", q_size = 2.1, show_names = show_names) +
    (if (!show_names) guides(fill = "none") else NULL)
}

rev_panels <- imap(REV_CONTRASTS, \(id, nm) rev_panel(id, nm, show_names = nm == "Disease"))
fig_reversal <- wrap_plots(rev_panels, nrow = 1) +
  plot_layout(guides = "collect", widths = c(1.7, 1, 1)) +
  plot_annotation(
    title = "Disease-selected pathway sets read across contrasts",
    subtitle = paste0(
      "The four ECM and four proliferation sets phenylephrine moves most under Disease, ",
      "read under Disease, Rescue and Transplant (names on the Disease column; rows shared). ",
      "One set clears FDR < 0.05 under Rescue; Disease and Rescue share the Phe group and are ",
      "coupled by design, so this apparent reversal is descriptive only. fgsea NES, per-database BH-FDR."
    ),
    theme = FIG_SUB
  ) &
  bottom_key
ggsave(file.path(PNG, "SUPP_F02_reversal.png"), fig_reversal,
  width = 250, height = 110, units = "mm", dpi = 300, bg = "white"
)

## Supplementary workbook ----------------------------------------------------
sheet_long <- function(figure, panels) {
  imap(panels, \(d, nm) transmute(d,
    figure = figure, panel = nm, pathway, database, theme, contrast,
    NES = round(NES, 3), pval = signif(pval, 3), padj = signif(padj, 3), size,
    sig = sig == "TRUE"
  )) |>
    list_rbind()
}
by_contrast_long <- sheet_long("enrichment_by_contrast", by_contrast)

append_workbook(
  file.path(BASE, "c_data", "F02_supplementary.xlsx"),
  sheet_specs = list(list(
    name = "enrichment_supp", df = as.data.frame(by_contrast_long),
    role = "Consolidated enrichment supplement: sets shown per contrast panel",
    contents = "figure, panel, pathway, database, theme, contrast, NES, p, FDR, set size, significance"
  ))
)
write_csv(by_contrast_long, file.path(DAT_SUPP, "SUPP_F02_enrichment_by_contrast.csv"))

message("F02 enrichment supplement: 3 figures -> b_reports/supp")

mirror_to_box(
  list.files(PNG, "^SUPP_F02_(enrichment_by_contrast|disease_specificity|reversal)\\.png$", full.names = TRUE),
  "02_Figures/F02_Enrich_Volcanoes/supp",
  guard = TRUE
)
mirror_to_box(file.path(BASE, "c_data", "F02_supplementary.xlsx"), "03_Supplementary", guard = TRUE)
