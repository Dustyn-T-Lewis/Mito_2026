#!/usr/bin/env Rscript
# F03 Venn — area-proportional 3-set overlap of DEP membership (Pi < 0.05) across
# Disease (CTLvPHE), Transplant (CTLvMITO), and Rescue (PHEvPHE_MITO), with a
# companion Up/Down directional strip. Reads existing pipeline outputs only.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(readr); library(ggplot2); library(eulerr)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE    <- here::here("04_Figures_v2", "03_Venn")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

# --- set definitions -------------------------------------------------------
SET_CONTRASTS <- c(Disease = "CTLvPHE", Transplant = "CTLvMITO", Rescue = "PHEvPHE_MITO")
SET_COLORS    <- unname(CONTRAST_COLORS[SET_CONTRASTS])  # red / blue / green
names(SET_COLORS) <- names(SET_CONTRASTS)

# --- load long DEP table, keep the 3 contrasts -----------------------------
comb <- read_csv(P05$comb, show_col_types = FALSE) |>
  filter(contrast %in% SET_CONTRASTS, !is.na(pi_score)) |>
  mutate(set = names(SET_CONTRASTS)[match(contrast, SET_CONTRASTS)])

# Significant membership per set (Pi < threshold), with logFC direction.
sig <- comb |>
  filter(pi_score < H9C2_PI_THRESH) |>
  transmute(uniprot_id, gene, set,
            direction = if_else(logFC >= 0, "Up", "Down"))

set_lists <- split(sig$uniprot_id, sig$set)[names(SET_CONTRASTS)]
set_lists <- lapply(set_lists, function(x) if (is.null(x)) character(0) else x)
set_sizes <- vapply(set_lists, length, integer(1))

# --- per-protein membership wide table + region key ------------------------
membership <- sig |>
  pivot_wider(id_cols = c(uniprot_id, gene), names_from = set,
              values_from = direction) |>
  # ensure all 3 set columns exist even if a set is empty
  {\(d) { for (s in names(SET_CONTRASTS)) if (!s %in% names(d)) d[[s]] <- NA_character_; d }}() |>
  mutate(
    Disease_member    = !is.na(Disease),
    Transplant_member = !is.na(Transplant),
    Rescue_member     = !is.na(Rescue),
    region_key = paste0(
      if_else(Disease_member,    "Disease",    ""),
      if_else(Transplant_member, "+Transplant", ""),
      if_else(Rescue_member,     "+Rescue",     "")) |>
      sub("^\\+", "", x = _)
  ) |>
  rename(Disease_dir = Disease, Transplant_dir = Transplant, Rescue_dir = Rescue) |>
  select(uniprot_id, gene,
         Disease = Disease_member, Transplant = Transplant_member, Rescue = Rescue_member,
         Disease_dir, Transplant_dir, Rescue_dir, region_key) |>
  arrange(region_key, gene)

# --- 7 region counts -------------------------------------------------------
region_counts <- membership |>
  count(region_key, name = "n") |>
  arrange(desc(n))

# --- area-proportional Venn (eulerr) ---------------------------------------
# eulerr fits the ellipse areas from disjoint region counts.
m <- membership
eu_fit <- euler(c(
  "Disease"                       = sum(m$Disease & !m$Transplant & !m$Rescue),
  "Transplant"                    = sum(!m$Disease & m$Transplant & !m$Rescue),
  "Rescue"                        = sum(!m$Disease & !m$Transplant & m$Rescue),
  "Disease&Transplant"            = sum(m$Disease & m$Transplant & !m$Rescue),
  "Disease&Rescue"                = sum(m$Disease & !m$Transplant & m$Rescue),
  "Transplant&Rescue"             = sum(!m$Disease & m$Transplant & m$Rescue),
  "Disease&Transplant&Rescue"     = sum(m$Disease & m$Transplant & m$Rescue)
), shape = "ellipse")

set_labels <- vapply(names(SET_CONTRASTS), function(s)
  sprintf("%s\n(%s)", s, contrast_brief(SET_CONTRASTS[[s]])), character(1))

venn_grob <- plot(
  eu_fit,
  fills  = list(fill = SET_COLORS, alpha = 0.5),
  edges  = list(col = SET_COLORS, lwd = 1.2),
  labels = list(labels = set_labels, fontsize = 6, fontfamily = "Helvetica", font = 2),
  quantities = list(fontsize = 6, fontfamily = "Helvetica"),
  legend = FALSE
)

# --- companion Up/Down directional strip -----------------------------------
# Per-set Up vs Down counts among that set's significant members, plus the
# all-three core intersection.
strip_dat <- sig |>
  count(set, direction, name = "n") |>
  mutate(group = factor(set, levels = names(SET_CONTRASTS)))

core_ids <- m$uniprot_id[m$Disease & m$Transplant & m$Rescue]
if (length(core_ids)) {
  core_dir <- sig |>
    filter(uniprot_id %in% core_ids) |>
    distinct(uniprot_id, set, direction) |>
    count(direction, name = "n") |>
    mutate(group = "Core (all 3)", set = "Core (all 3)")
  strip_dat <- bind_rows(strip_dat, core_dir)
}
strip_levels <- c(names(SET_CONTRASTS),
                  if (length(core_ids)) "Core (all 3)")
strip_dat$group <- factor(strip_dat$group, levels = strip_levels)
strip_dat$direction <- factor(strip_dat$direction, levels = c("Up", "Down"))

p_strip <- ggplot(strip_dat, aes(group, n, fill = direction)) +
  geom_col(position = position_dodge(preserve = "single", width = 0.8),
           width = 0.7, color = "grey25", linewidth = 0.15) +
  geom_text(aes(label = n),
            position = position_dodge(preserve = "single", width = 0.8),
            vjust = -0.3, size = 1.6, color = "grey15") +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Direction within set", x = NULL, y = "Proteins (Π < 0.05)") +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT),
        legend.position = "top", legend.key.size = unit(2.5, "mm"),
        plot.margin = margin(3, 2, 1, 1))

# --- combine: Venn (as ggplot) + strip via patchwork -----------------------
sub_txt <- sprintf("Disease (%s) | Transplant (%s) | Rescue (%s)",
                   contrast_brief("CTLvPHE"), contrast_brief("CTLvMITO"),
                   contrast_brief("PHEvPHE_MITO"))

have_ggplotify <- requireNamespace("ggplotify", quietly = TRUE)
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

if (have_ggplotify && have_patchwork) {
  venn_gg <- ggplotify::as.ggplot(venn_grob) +
    labs(title = "DEP overlap (Π < 0.05)", subtitle = sub_txt) +
    theme_void(base_family = "Helvetica") +
    theme(plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                       margin = margin(b = 1)),
          plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                       color = "grey30", margin = margin(b = 2)),
          plot.margin   = margin(3, 2, 1, 2))
  combined <- patchwork::wrap_plots(venn_gg, p_strip, ncol = 2, widths = c(1.6, 1))
  ggsave(file.path(RPT_PDF, "MAIN_F03_venn.pdf"), combined,
         width = 140, height = 100, units = "mm", device = pdf_dev)
  ggsave(file.path(RPT_PNG, "MAIN_F03_venn.png"), combined,
         width = 140, height = 100, units = "mm", dpi = 300)
} else {
  # Fallback: Venn as the main figure, strip as a companion file.
  pdf_dev(file.path(RPT_PDF, "MAIN_F03_venn.pdf"), width = 140 / 25.4, height = 100 / 25.4)
  grid::grid.newpage()
  grid::grid.text("DEP overlap (Π < 0.05)", y = 0.97,
                  gp = grid::gpar(fontface = "bold", fontsize = 8))
  grid::grid.text(sub_txt, y = 0.93, gp = grid::gpar(fontsize = 5, col = "grey30"))
  print(venn_grob, newpage = FALSE)
  dev.off()
  png(file.path(RPT_PNG, "MAIN_F03_venn.png"), width = 140, height = 100,
      units = "mm", res = 300)
  grid::grid.newpage()
  grid::grid.text("DEP overlap (Π < 0.05)", y = 0.97,
                  gp = grid::gpar(fontface = "bold", fontsize = 8))
  grid::grid.text(sub_txt, y = 0.93, gp = grid::gpar(fontsize = 5, col = "grey30"))
  print(venn_grob, newpage = FALSE)
  dev.off()
  ggsave(file.path(RPT_PDF, "MAIN_F03_venn_direction.pdf"), p_strip,
         width = 80, height = 70, units = "mm", device = pdf_dev)
  ggsave(file.path(RPT_PNG, "MAIN_F03_venn_direction.png"), p_strip,
         width = 80, height = 70, units = "mm", dpi = 300)
}

# --- supplementary workbook ------------------------------------------------
build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  figure_title = "F03 — DEP overlap (Π < 0.05) across Disease / Transplant / Rescue",
  sheet_specs = list(
    list(name = "membership",    df = membership,
         role     = "Per-protein set membership behind the Euler diagram",
         contents = "One row per significant protein (Π<0.05): TRUE/FALSE membership in Disease/Transplant/Rescue, the Up/Down direction in each set, and its region_key"),
    list(name = "region_counts", df = region_counts,
         role     = "The 7 Euler region areas + the directional strip counts",
         contents = "Protein count per overlap region (region_key), sorted descending")))

message(sprintf(
  "F03 Venn | sets: Disease=%d Transplant=%d Rescue=%d | %d sig proteins, %d regions",
  set_sizes[["Disease"]], set_sizes[["Transplant"]], set_sizes[["Rescue"]],
  nrow(membership), nrow(region_counts)))
print(region_counts)
