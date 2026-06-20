#!/usr/bin/env Rscript
# F04 MAIN — Panel-D-style pathway count summary across 5 DBs.
# Per contrast, Up bar and Down bar drawn from 0; mito subset overdrawn darker.
# Pool = CANONICAL_DBS (Hallmark + Reactome + KEGG + MitoCarta + GO Slim);
# padj<0.05, size>=10, MITO_DROP_SETS excluded; cross-DB Jaccard+Overlap dedup
# at 0.375 per contrast (EnrichmentMap, Merico 2010 / Reimand 2019).
# Mito flag: database == MitoCarta OR pathway matches MITO_PATHWAY_REGEX.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(scales)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE    <- here::here("04_Figures_v2", "04_Pathway_bars")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

CORE      <- c("CTLvMITO", "CTLvPHE", "PHEvPHE_MITO", "Interaction")
PADJ_CUT  <- 0.05
MIN_SIZE  <- 10
SIM_CUT   <- 0.375

fgsea_all     <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
                          show_col_types = FALSE)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
set_pool      <- do.call(c, unname(rat_gene_sets[CANONICAL_DBS]))

# ---- Per-contrast: pool 5 DBs, dedup, classify mito ---------------------------
per_contrast <- function(ctr) {
  sig <- fgsea_all |>
    filter(contrast == ctr, database %in% CANONICAL_DBS,
           !is.na(padj), padj < PADJ_CUT,
           size >= MIN_SIZE, !pathway %in% MITO_DROP_SETS) |>
    arrange(padj)
  if (nrow(sig) > 1)
    sig <- deduplicate_enrichment(as.data.frame(sig), pathways = set_pool,
                                  cutoff = SIM_CUT, cross_db = TRUE) |> as_tibble()
  sig |>
    mutate(direction = if_else(NES > 0, "Up", "Down"),
           is_mito = database == "MitoCarta" |
                     grepl(MITO_PATHWAY_REGEX, pathway, perl = TRUE),
           contrast = ctr)
}
sig_pw <- bind_rows(lapply(CORE, per_contrast))
message(sprintf("F04: %d total post-dedup sig pathways across %d contrasts",
                nrow(sig_pw), length(CORE)))

# ---- Bar data: contrast × direction → total + mito subset --------------------
bar_df <- sig_pw |>
  summarise(total = n(), mito = sum(is_mito), .by = c(contrast, direction)) |>
  complete(contrast = CORE, direction = c("Up", "Down"),
           fill = list(total = 0L, mito = 0L)) |>
  mutate(contrast = factor(contrast, levels = CORE),
         direction = factor(direction, levels = c("Up", "Down")))

# Fills: light = total, dark = mito subset.
FILL_TOTAL <- c(Up = "#F4A582", Down = "#92C5DE")
FILL_MITO  <- c(Up = "#B2182B", Down = "#2166AC")

BAR_W   <- 0.34; GAP <- 0.06
centers <- setNames(seq_along(CORE), CORE)
bar_df <- bar_df |>
  mutate(x_center  = centers[as.character(contrast)] +
                     ifelse(direction == "Up", -(BAR_W/2 + GAP/2), (BAR_W/2 + GAP/2)),
         fill_tot  = FILL_TOTAL[as.character(direction)],
         fill_mito = FILL_MITO[as.character(direction)])

bg_df <- tibble(xmin = seq_along(CORE) - 0.5, xmax = seq_along(CORE) + 0.5,
                fill = unname(CONTRAST_COLORS[CORE]))

p <- ggplot() +
  geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = bg_df$fill, alpha = 0.16, color = "grey75", linewidth = 0.2) +
  geom_rect(data = bar_df,
            aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2,
                ymin = 0, ymax = total),
            fill = bar_df$fill_tot, color = "black", linewidth = 0.3) +
  geom_rect(data = filter(bar_df, mito > 0),
            aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2,
                ymin = 0, ymax = mito),
            fill = filter(bar_df, mito > 0)$fill_mito,
            color = "black", linewidth = 0.3) +
  geom_text(data = filter(bar_df, total > 0),
            aes(x = x_center, y = total, label = total),
            vjust = -0.4, size = 2.4, fontface = "bold") +
  scale_x_continuous(breaks = seq_along(CORE),
                     labels = setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
                     expand = expansion(mult = 0)) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.08)),
               breaks = c(5, 25, 50, 100, 200, 300, 500)) +
  coord_cartesian(clip = "off") +
  labs(title    = "Pathway enrichment (Up / Down)",
       subtitle = "Hallmark + Reactome + KEGG + MitoCarta + GO Slim | dark = mito | sqrt y, cross-DB Jaccard-deduped",
       x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(axis.text.x = element_text(size = FIG_AXIS_TEXT, lineheight = 0.85, face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(2, 4, 2, 2))

# Inline 4-key legend (matches Panel D).
key_df <- tibble(y = 4:1,
                 lab  = c("Up total", "Up mito", "Down total", "Down mito"),
                 fill = c(FILL_TOTAL["Up"], FILL_MITO["Up"],
                          FILL_TOTAL["Down"], FILL_MITO["Down"]))
p_key <- ggplot(key_df) +
  geom_point(aes(0, y), shape = 22, size = 2, fill = key_df$fill,
             color = "grey30", stroke = 0.3) +
  geom_text(aes(0.22, y, label = lab), hjust = 0, size = 1.6, color = "grey20") +
  scale_x_continuous(limits = c(-0.1, 2.6)) +
  scale_y_continuous(limits = c(0.5, 4.5)) +
  theme_void()
p <- p + patchwork::inset_element(p_key, left = 0.02, right = 0.40,
                                  top = 0.99, bottom = 0.62)

FIG_W <- 120; FIG_H <- 70
ggsave(file.path(RPT_PDF, "MAIN_F04_pathway_bars.pdf"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F04_pathway_bars.png"), p,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE)

# ---- Workbook ----------------------------------------------------------------
counts_tab <- bar_df |>
  transmute(contrast = contrast_brief(as.character(contrast)),
            direction, total, mito,
            mito_fraction = ifelse(total > 0, mito / total, NA_real_)) |>
  arrange(contrast, direction)

per_ctr_tabs <- lapply(CORE, function(ctr) {
  sig_pw |> filter(contrast == ctr) |>
    transmute(pathway, database, clean_label = clean_display_label(pathway),
              NES = round(NES, 3), padj = signif(padj, 3),
              size, direction, is_mito) |>
    arrange(padj)
})

sheet_specs <- c(
  list(list(name = "dep_pathway_counts", df = counts_tab,
            role     = "F04 bar heights — total significant pathways and mito subset per contrast x direction",
            contents = "contrast, direction (Up/Down), total significant pathways, mito subset, mito_fraction")),
  Map(function(ctr, tab) list(
        name = paste0(contrast_brief(ctr), "_sig_pathways"),
        df = tab,
        role     = sprintf("Post-dedup significant pathways behind the %s bars", contrast_brief(ctr)),
        contents = "pathway, database, clean_label, NES, padj, size, direction (Up/Down), is_mito"),
      CORE, per_ctr_tabs)
)

build_workbook(file.path(DAT, "F04_supplementary.xlsx"),
               figure_title = "F04 — Panel-D pathway count summary (5-DB pool, mito subset overlaid)",
               sheet_specs = sheet_specs)

# Remove the legacy shown_pathways.csv if a prior run left it behind — no consumer now.
old_csv <- file.path(DAT, "shown_pathways.csv")
if (file.exists(old_csv)) file.remove(old_csv)

message(sprintf("F04 done | %d sig pathways shown | output %s", nrow(sig_pw), RPT_PNG))
