# F01 — Panel F: Pathway enrichment (Total + Mito two-layer stack, Up/Down).
# Per contrast, Up bar (red) and Down bar (blue); each bar split into
#   non-mito pathways (lighter) + mitochondrial pathways (darker, emphasized).
# Mito = MitoCarta sets OR broad-DB pathways matching MITO_PATHWAY_REGEX.
# Reads the frozen rat fGSEA cache (msigdbr rat orthologs + MitoCarta mouse->rat).
# Sourced from 01_main_panels.R. Inherits: PNL_PNG, PNL_PDF, DAT, pdf_dev, CORE,
#   CTR_LABEL_2L, CONTRAST_COLORS, DIR_COLORS, contrast_brief, FIG_* sizes.

source(here::here("04_Figures", "shared", "mitocarta_utils.R"))
source(here::here("04_Figures", "shared", "pathway_utils.R"))

fgsea_cache <- here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv")
stopifnot("fGSEA cache missing" = file.exists(fgsea_cache))
fgsea_raw <- read_csv(fgsea_cache, show_col_types = FALSE)

# DB set matches F02 (Hallmark + Reactome + MitoCarta) so F01 counts reconcile
# with the F02 volcano rings. GO:BP + KEGG dropped (parsimony, lower MTC burden).
BROAD_DBS <- c("Hallmark", "Reactome")   # rat orthologs (msigdbr)
ALL_DBS   <- c(BROAD_DBS, "MitoCarta")

# Significant pathways across the appropriate rat databases.
sig_raw <- fgsea_raw |>
  filter(!is.na(padj), padj < 0.05, database %in% ALL_DBS, contrast %in% CORE) |>
  mutate(direction = ifelse(NES > 0, "Up", "Down"),
         is_mito = database == "MitoCarta" |
           grepl(MITO_PATHWAY_REGEX, pathway, perl = TRUE))

# Cross-DB Jaccard dedup per contrast (collapses redundant terms) so counts are
# parsimonious and consistent with the F02 volcano rings (same DB universe).
gs_all <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
pw_all <- c(gs_all$Hallmark, gs_all$Reactome, gs_all$MitoCarta)
sig_pw <- bind_rows(lapply(CORE, function(ctr) {
  d <- as.data.frame(sig_raw[sig_raw$contrast == ctr, ])
  if (nrow(d) > 1) d <- deduplicate_enrichment(d, pathways = pw_all, jaccard_cutoff = 0.5, cross_db = TRUE)
  as_tibble(d)
}))
message(sprintf("F01 enrichment: %d sig pathways -> %d after cross-DB dedup", nrow(sig_raw), nrow(sig_pw)))

write_csv(sig_pw |> select(contrast, database, pathway, padj, NES, direction, is_mito),
          file.path(DAT, "panel_F_enrichment_sig.csv"))

# Counts: contrast x direction x {Mito, Non-mito}
count_df <- sig_pw |>
  mutate(layer = ifelse(is_mito, "Mito", "Non-mito")) |>
  count(contrast, direction, layer, name = "count")

full_grid <- expand_grid(contrast = CORE, direction = c("Up", "Down"),
                         layer = c("Non-mito", "Mito"))
count_df <- full_grid |>
  left_join(count_df, by = c("contrast", "direction", "layer")) |>
  mutate(count = tidyr::replace_na(count, 0L),
         contrast = factor(contrast, levels = CORE),
         direction = factor(direction, levels = c("Up", "Down")),
         layer = factor(layer, levels = c("Non-mito", "Mito")))

# Per contrast x direction: total sig pathways + mito subset (overlaid bars,
# both drawn from 0 so a sqrt y-axis keeps small bars visible without the
# additivity distortion that a stacked bar would suffer on a nonlinear scale).
bar_df <- count_df |>
  summarise(total = sum(count), mito = sum(count[layer == "Mito"]),
            .by = c(contrast, direction))

# Fills: total = lighter shade, mito = darker (emphasized) per direction.
FILL_TOTAL <- c(Up = "#F4A582", Down = "#92C5DE")
FILL_MITO  <- c(Up = "#B2182B", Down = "#2166AC")

BAR_W <- 0.34; GAP <- 0.06
ctr_centers <- setNames(seq_along(CORE), CORE)
bar_df <- bar_df |>
  mutate(contrast = factor(contrast, levels = CORE),
         x_center = ctr_centers[as.character(contrast)] +
           ifelse(direction == "Up", -(BAR_W/2 + GAP/2), (BAR_W/2 + GAP/2)),
         fill_total = FILL_TOTAL[as.character(direction)],
         fill_mito  = FILL_MITO[as.character(direction)])

PF_W <- 78; PF_H <- 56
y_max <- max(bar_df$total) * 1.18
bg_rects <- tibble(xmin = seq_along(CORE) - 0.5, xmax = seq_along(CORE) + 0.5,
                   fill = unname(CONTRAST_COLORS[CORE]))

pF <- ggplot() +
  geom_rect(data = bg_rects, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = bg_rects$fill, alpha = 0.16, color = "grey75", linewidth = 0.2) +
  # total bar (light) then mito bar (dark) — both from 0, sqrt axis keeps both readable
  geom_rect(data = bar_df, aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2,
                               ymin = 0, ymax = total), fill = bar_df$fill_total,
            color = "black", linewidth = 0.3) +
  geom_rect(data = bar_df |> filter(mito > 0),
            aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2, ymin = 0, ymax = mito),
            fill = (bar_df |> filter(mito > 0))$fill_mito, color = "black", linewidth = 0.3) +
  geom_text(data = bar_df |> filter(total > 0), aes(x = x_center, y = total, label = total),
            vjust = -0.4, size = scale_text(BASE_COUNT, PF_W) + 0.6, fontface = "bold") +
  scale_x_continuous(breaks = seq_along(CORE),
                     labels = setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
                     expand = expansion(mult = 0)) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.06)),
               breaks = c(5, 25, 50, 100, 200, 300)) +
  coord_cartesian(clip = "off") +
  labs(title = "Pathway enrichment (Up / Down)",
       subtitle = "Hallmark·Reactome·MitoCarta (rat) | dark = mito | √-scaled, Jaccard-deduped",
       x = NULL, y = "Significant pathways", tag = "D") +
  FIG_THEME +
  theme(axis.text.x = element_text(size = FIG_AXIS_TEXT - 1, lineheight = 0.8, face = "bold"),
        panel.grid.major.x = element_blank(), plot.margin = margin(1, 4, 1, 2))

# Stacked key (4x1) — same point/text size + spacing as panel E's key.
key_df <- tibble(
  y = 4:1,
  lab = c("Up total", "Up mito", "Down total", "Down mito"),
  fill = c(FILL_TOTAL["Up"], FILL_MITO["Up"], FILL_TOTAL["Down"], FILL_MITO["Down"]))
p_key <- ggplot(key_df) +
  geom_point(aes(0, y), shape = 22, size = 1.7, fill = key_df$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(0.2, y, label = lab), size = 1.4, hjust = 0, color = "grey20") +
  scale_x_continuous(limits = c(-0.1, 2.6)) + scale_y_continuous(limits = c(0.5, 4.5)) +
  theme_void()
pF <- pF + inset_element(p_key, left = 0.02, right = 0.42, top = 0.99, bottom = 0.64)

n_mito_sig <- sum(sig_pw$is_mito); n_sig <- nrow(sig_pw)
ggsave(file.path(PNL_PNG, "MAIN_panel_F_enrichment.png"), pF,
       width = PF_W, height = PF_H, units = "mm", dpi = 300)
message(sprintf("F01 Panel F (enrichment Total+Mito) saved | %d sig pathways, %d mito", n_sig, n_mito_sig))
