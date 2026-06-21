#!/usr/bin/env Rscript
# F02 MAIN — DEP bars + effect-size companion. Standalone (no composite).
# Reads existing pipeline outputs only; never re-runs 01-03.
#   Left  : DEP counts as % of proteome, per contrast, at three independent
#           thresholds (p<0.05 / FDR<0.10 / Π<0.05), Up/Down dodged, contrast-
#           shaded background (YvO-style), faint->solid by threshold.
#   Right : |log2FC| distribution per contrast — histogram + density faceted by
#           contrast, median |log2FC| annotated per facet (mirrors F01 Panel C).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(readr); library(ggplot2); library(patchwork); library(scales)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE    <- here::here("04_Figures_v2", "02_DEP_bars")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

CORE     <- c("CTLvMITO", "CTLvPHE", "PHEvPHE_MITO", "Interaction")
comb     <- load_combined_wide()
comb_long <- read_csv(P05$comb, show_col_types = FALSE)
dep_results <- setNames(lapply(CORE, \(c) as.data.frame(filter(comb_long, contrast == c))), CORE)

# proteome size = unique genes measured
all_genes <- unique(comb$gene[!is.na(comb$gene)]); n_total <- length(all_genes)

CTR_LAB    <- setNames(contrast_brief(CORE), CORE)
ctr_levels <- unname(CTR_LAB)                              # Transplant -> Interaction
THR_LEVELS <- c("p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), "Π < 0.05")

# ---- Panel 1: DEP counts (Up/Down) at three independent thresholds -----------
# Up = logFC > 0 & sig; Down = logFC < 0 & sig, at each threshold independently.
sig_flag <- function(r, thr) {
  switch(thr,
    "p"   = !is.na(r$P.Value)   & r$P.Value   < 0.05,
    "fdr" = !is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR,
    "pi"  = !is.na(r$pi_score)  & r$pi_score  < H9C2_PI_THRESH)
}
counts_df <- bind_rows(lapply(CORE, function(c) {
  r <- dep_results[[c]]
  bind_rows(lapply(c(p = "p", fdr = "fdr", pi = "pi"), function(thr) {
    s  <- sig_flag(r, thr)
    up <- sum(s & r$logFC > 0); dn <- sum(s & r$logFC < 0)
    tibble(contrast = CTR_LAB[c], threshold = thr,
           direction = c("Up", "Down"), n = c(up, dn))
  }), .id = NULL)
})) |>
  mutate(threshold = recode(threshold, p = THR_LEVELS[1],
                            fdr = THR_LEVELS[2], pi = THR_LEVELS[3]),
         contrast  = factor(contrast, levels = ctr_levels),
         threshold = factor(threshold, levels = THR_LEVELS),
         direction = factor(direction, levels = c("Up", "Down")),
         pct       = 100 * n / n_total)

# Fill key: contrast colour, faint -> solid across the three thresholds.
ALPHAS <- c(0.30, 0.60, 1.00)
FILL_KEYS <- counts_df |>
  distinct(contrast, threshold) |>
  mutate(key = paste(contrast, threshold, sep = "___"))
FILL_VALS <- setNames(
  vapply(seq_len(nrow(FILL_KEYS)), function(i) {
    ctr <- as.character(FILL_KEYS$contrast[i])
    col <- unname(CONTRAST_COLORS[CORE[match(ctr, ctr_levels)]])
    adjustcolor(col, alpha.f = ALPHAS[as.integer(FILL_KEYS$threshold[i])])
  }, character(1)),
  FILL_KEYS$key)

counts_df <- counts_df |> mutate(fill_key = paste(contrast, threshold, sep = "___"))

# accent background band per contrast group (YvO-style)
bg_df <- tibble(contrast = factor(ctr_levels, levels = ctr_levels),
                fill = unname(CONTRAST_COLORS[CORE]),
                xi = seq_along(ctr_levels))

# x position: contrast group, dodged by direction; threshold nested via fill alpha
# We lay out groups along x; within each group Up bars then Down bars, dodged by
# threshold. Build explicit x to keep p/FDR/Π adjacent and Up/Down separated.
dodge_w <- 0.92
counts_df <- counts_df |>
  mutate(dir_off = ifelse(direction == "Up", -1, 1) * 0.24,
         thr_off = (as.integer(threshold) - 2) * 0.13,
         xc = as.integer(contrast) + dir_off + thr_off)

# Up/Down direction tint marker strip under axis via secondary annotation; here
# we simply outline Down bars to distinguish from Up at a glance.
y_max <- max(counts_df$pct) * 1.12

p_counts <- ggplot(counts_df, aes(xc, pct, fill = fill_key)) +
  geom_rect(data = bg_df, aes(xmin = xi - 0.5, xmax = xi + 0.5,
            ymin = -Inf, ymax = Inf), fill = bg_df$fill, alpha = 0.14,
            inherit.aes = FALSE) +
  geom_col(aes(linetype = direction), width = 0.115, color = "grey20",
           linewidth = 0.25) +
  geom_vline(xintercept = head(seq_along(ctr_levels), -1) + 0.5,
             color = "grey85", linewidth = 0.25) +
  scale_fill_manual(values = FILL_VALS, guide = "none") +
  scale_linetype_manual(values = c(Up = "solid", Down = "22"),
                        name = NULL,
                        guide = guide_legend(override.aes = list(
                          fill = "grey60", color = "grey20"))) +
  scale_x_continuous(breaks = seq_along(ctr_levels), labels = ctr_levels,
                     expand = expansion(add = 0.06)) +
  scale_y_continuous(limits = c(0, y_max),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "DEP counts by contrast",
       subtitle = "Up / Down dodged; p / FDR / Π faint→solid (independent thresholds)",
       x = NULL, y = "% of proteome", tag = "A") +
  FIG_THEME +
  theme(legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.key.size = unit(2.5, "mm"),
        plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, face = "italic",
                                     color = "grey40"),
        axis.text.x = element_text(face = "bold", size = FIG_AXIS_TEXT),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(5, 3, 1, 2))

# Small in-panel threshold key (faint->solid swatches, neutral grey).
thr_key <- tibble(y = 3:1, lab = THR_LEVELS,
                  fill = vapply(ALPHAS, \(a) adjustcolor("grey25", alpha.f = a),
                                character(1)))
p_thrkey <- ggplot(thr_key) +
  geom_point(aes(0, y), shape = 22, size = 2, fill = thr_key$fill,
             color = "grey40", stroke = 0.3) +
  geom_text(aes(0.25, y, label = lab), hjust = 0, size = 1.5, color = "grey20") +
  scale_x_continuous(limits = c(-0.2, 2.4)) +
  scale_y_continuous(limits = c(0.4, 3.6)) + theme_void()
p_counts <- p_counts +
  inset_element(p_thrkey, left = 0.62, right = 0.99, top = 0.78, bottom = 0.52)

# ---- Panel 2: effect-size companion — |log2FC| dist per contrast -------------
lfc_long <- bind_rows(lapply(CORE, \(c)
  tibble(contrast = CTR_LAB[c], logFC = dep_results[[c]]$logFC))) |>
  filter(!is.na(logFC), abs(logFC) <= 1) |>
  mutate(contrast = factor(contrast, levels = ctr_levels))
lfc_stats <- lfc_long |>
  summarise(med_abs = median(abs(logFC)), .by = contrast) |>
  mutate(lab = sprintf("med|LFC| %.2f", med_abs))

hbw <- 2 / 44
hist_bg <- tibble(contrast = factor(ctr_levels, levels = ctr_levels),
                  fill = unname(CONTRAST_COLORS[CORE]))
p_eff <- ggplot(lfc_long, aes(logFC)) +
  geom_rect(data = hist_bg, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = hist_bg$fill, alpha = 0.14, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_histogram(aes(fill = contrast), breaks = seq(-1, 1, by = hbw),
                 color = "white", linewidth = 0.1, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * hbw), color = "grey20", linewidth = 0.4) +
  geom_text(data = lfc_stats, aes(x = -0.95, y = Inf, label = lab),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4,
            size = scale_text(BASE_STAT, 50) + 0.3, fontface = "bold",
            color = "grey25") +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_fill_manual(values = setNames(unname(CONTRAST_COLORS[CORE]), ctr_levels),
                    guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1)) + scale_y_continuous(breaks = NULL) +
  labs(title = "Effect size", subtitle = "signed log2FC per contrast",
       x = expression(bold(log[2] ~ FC)), y = NULL, tag = "B") +
  FIG_THEME +
  theme(strip.text.y = element_text(face = "bold", size = FIG_STRIP_SIZE - 0.5,
                                    angle = 0),
        strip.background = element_blank(),
        plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, face = "italic",
                                     color = "grey40"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = FIG_AXIS_TEXT),
        panel.grid = element_blank(), panel.spacing.y = unit(1.5, "pt"),
        plot.margin = margin(5, 2, 1, 1))

# ---- Compose -----------------------------------------------------------------
fig <- p_counts + p_eff + plot_layout(widths = c(1, 0.55))

FIG_W <- 178; FIG_H <- 100
ggsave(file.path(RPT_PDF, "MAIN_F02_dep_bars.pdf"), fig,
       width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F02_dep_bars.png"), fig,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE)

# ---- Supplementary workbook --------------------------------------------------
dep_tabs <- lapply(CORE, function(c) {
  dep_results[[c]] |>
    transmute(uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score) |>
    arrange(pi_score)
})
sheet_specs <- c(
  list(list(name = "dep_counts",
            df = counts_df |> select(contrast, threshold, direction, n, pct),
            role     = "Panel A bar heights — DEP counts per contrast",
            contents = "Up/Down significant-protein counts and % of proteome at each independent threshold (p<0.05, FDR, Π<0.05) per contrast")),
  Map(function(c, tab) list(name = substr(contrast_brief(c), 1, 28), df = tab,
            role     = "Underlying per-protein DE table feeding both panels",
            contents = sprintf("Every protein tested in the %s contrast: uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score (sorted by Π)", contrast_brief(c))),
      CORE, dep_tabs)
)
build_workbook(file.path(DAT, "F02_supplementary.xlsx"),
               figure_title = "F02 — DEP counts by contrast + effect-size (|log2FC|) companion",
               sheet_specs = sheet_specs)

pi_summary <- counts_df |> filter(threshold == THR_LEVELS[3]) |>
  select(contrast, direction, n) |> pivot_wider(names_from = direction, values_from = n)
message(sprintf("F02 MAIN (%dx%d mm) | proteome n=%d", FIG_W, FIG_H, n_total))
print(as.data.frame(pi_summary))
