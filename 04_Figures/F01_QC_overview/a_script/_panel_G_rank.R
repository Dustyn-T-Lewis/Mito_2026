# F01 — Panel F: DEP rank location (barcode), formatted exactly like YvO F02 panel F.
# Density ribbons (Up/Down) + direction-coloured barcode ticks + descriptive
# directional callouts with leader lines to each density peak. Faceted by contrast.
# H9c2 additions: mitochondrial DEPs marked with longer black ticks + mito count
# in the subtitle. Sourced from 01_main_panels.R.

source(here::here("04_Figures", "shared", "mitocarta_utils.R"))

mito_genes <- {
  mc <- load_mitocarta_collection(drop_all = FALSE)
  unique(toupper(unlist(mc, use.names = FALSE)))
}
is_mito_gene <- function(g) toupper(g) %in% mito_genes

PG_W <- 74; PG_H <- 56
DENS_PAD <- 0.06

rank_df <- bind_rows(lapply(CORE, function(ctr) {
  t_col <- paste0("t_", ctr); pi_col <- paste0("pi_score_", ctr); lfc_col <- paste0("logFC_", ctr)
  comb |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(rank_frac = seq_len(n()) / n(),
           is_dep = !is.na(.data[[pi_col]]) & .data[[pi_col]] < H9C2_PI_THRESH,
           direction = case_when(!is_dep ~ NA_character_, .data[[lfc_col]] > 0 ~ "Up", TRUE ~ "Down"),
           is_mito = is_mito_gene(gene), contrast = ctr) |>
    select(gene, contrast, rank_frac, is_dep, direction, is_mito)
}))
rank_df$contrast <- factor(rank_df$contrast, levels = CORE)
dep_only <- rank_df |> filter(is_dep)
dep_only$direction <- factor(dep_only$direction, levels = c("Up", "Down"))

dep_counts <- dep_only |>
  summarise(n_up = sum(direction == "Up"), n_down = sum(direction == "Down"),
            n_total = n(), n_mito = sum(is_mito), .by = contrast) |>
  mutate(contrast = factor(contrast, levels = CORE))
write_csv(dep_counts, file.path(DAT, "panel_G_rank_counts.csv"))

# Density per contrast, peak-normalised to 1.
dens_list <- lapply(split(dep_only, dep_only$contrast), function(cdf) {
  lapply(split(cdf, cdf$direction, drop = TRUE), function(dd) {
    if (nrow(dd) < 2) return(NULL)
    d <- density(dd$rank_frac, adjust = 1.8, from = -DENS_PAD, to = 1 + DENS_PAD, n = 512)
    tibble(x = d$x, y = d$y, direction = dd$direction[1], contrast = dd$contrast[1])
  }) |> bind_rows()
}) |> bind_rows()
dens_list <- dens_list |> group_by(contrast) |> mutate(y_norm = y / max(y)) |> ungroup()
dens_list$direction <- factor(dens_list$direction, levels = c("Up", "Down"))
dens_list$contrast  <- factor(dens_list$contrast, levels = CORE)

TICK_DEPTH <- -0.25
ANNOT_SZ <- scale_text(BASE_STAT + 0.6, PG_W)
LABEL_NUDGE <- 0.06

peak_pos <- dens_list |> group_by(contrast, direction) |>
  slice_max(y_norm, n = 1, with_ties = FALSE) |> ungroup() |>
  select(contrast, direction, peak_x = x, peak_y = y_norm)
n_down <- dep_only |> filter(direction == "Down") |> count(contrast) |> tibble::deframe()
n_up   <- dep_only |> filter(direction == "Up")   |> count(contrast) |> tibble::deframe()

# H9c2 directional descriptions (mirror YvO panel F's plain-language callouts).
DESC_DOWN <- c(CTLvMITO = "down with Mito", CTLvPHE = "down with PHE",
               PHEvPHE_MITO = "rescued down in PHE", Interaction = "Mito stronger in healthy")
DESC_UP   <- c(CTLvMITO = "up with Mito", CTLvPHE = "up with PHE",
               PHEvPHE_MITO = "rescued up in PHE", Interaction = "Mito stronger in stressed")

bg_wash <- tibble(contrast = factor(CORE, levels = CORE), fill = unname(CONTRAST_COLORS[CORE]),
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
ad_all <- peak_pos |> filter(direction == "Down") |>
  mutate(ctr = as.character(contrast), label_x = peak_x + LABEL_NUDGE, label_y = peak_y * 0.78,
         label = paste(n_down[ctr], DESC_DOWN[ctr])) |> filter(!is.na(label), !is.na(n_down[ctr]))
au_all <- peak_pos |> filter(direction == "Up") |>
  mutate(ctr = as.character(contrast), label_x = peak_x - LABEL_NUDGE,
         label_y = if_else(ctr == "Interaction", 0.30, peak_y * 0.85),
         label = paste(n_up[ctr], DESC_UP[ctr])) |> filter(!is.na(label), !is.na(n_up[ctr]))
cd_all <- ad_all |> mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)
cu_all <- au_all |> mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)

pG <- ggplot() +
  geom_rect(data = bg_wash, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_wash$fill, alpha = 0.18, inherit.aes = FALSE) +
  geom_ribbon(data = dens_list, aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
              alpha = 0.30, outline.type = "upper") +
  geom_line(data = dens_list, aes(x = x, y = y_norm, color = direction), linewidth = 0.5) +
  geom_segment(data = dep_only, aes(x = rank_frac, xend = rank_frac, y = 0, yend = TICK_DEPTH,
                                    color = direction), linewidth = 0.35, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  { if (nrow(cd_all) > 0) geom_segment(data = cd_all, aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Down"]), alpha = 0.4, inherit.aes = FALSE) } +
  { if (nrow(ad_all) > 0) geom_label(data = ad_all, aes(x = label_x, y = label_y, label = label),
                 hjust = 0, vjust = 0.5, size = ANNOT_SZ, fill = unname(DIR_COLORS["Down"]), color = "white",
                 fontface = "bold", linewidth = 0, label.padding = unit(0.08, "lines"), inherit.aes = FALSE) } +
  { if (nrow(cu_all) > 0) geom_segment(data = cu_all, aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Up"]), alpha = 0.4, inherit.aes = FALSE) } +
  { if (nrow(au_all) > 0) geom_label(data = au_all, aes(x = label_x, y = label_y, label = label),
                 hjust = 1, vjust = 0.5, size = ANNOT_SZ, fill = unname(DIR_COLORS["Up"]), color = "white",
                 fontface = "bold", linewidth = 0, label.padding = unit(0.08, "lines"), inherit.aes = FALSE) } +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]), Down = unname(DIR_COLORS["Down"]))) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]), Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(-DENS_PAD, 1 + DENS_PAD), ylim = c(TICK_DEPTH, 1.02)) +
  facet_grid(contrast ~ ., switch = "y",
             labeller = labeller(contrast = setNames(gsub("_", "\n", contrast_brief(CORE)), CORE))) +
  labs(title = "DEP rank location",
       subtitle = sprintf("%s proteins | %d Π DEPs (t-ranked)",
                          format(length(unique(rank_df$gene)), big.mark = ","), sum(dep_counts$n_total)),
       x = "Rank position (by t-statistic)", y = NULL, tag = "F") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE + 1, face = "bold.italic", color = "grey40"),
        legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = FIG_AXIS_TEXT + 0.5),
        strip.text.y.left = element_text(face = "bold", size = FIG_AXIS_TEXT + 0.5, angle = 0,
                                         hjust = 1, vjust = 0.5, lineheight = 0.85),
        strip.background = element_blank(), strip.placement = "outside",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.y = unit(2, "pt"), plot.margin = margin(5, 4, 1, 2))

ggsave(file.path(PNL_PNG, "MAIN_panel_G_rank.png"), pG, width = PG_W, height = PG_H, units = "mm", dpi = 300)
message(sprintf("F01 Panel F (rank location, YvO format) saved | %d mito DEPs", sum(dep_counts$n_mito)))
