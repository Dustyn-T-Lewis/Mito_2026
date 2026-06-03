#!/usr/bin/env Rscript
# F01 MAIN composite — QC + proteome-wide picture. Reads existing pipeline
# outputs only; never re-runs 01-03.
# Order: A PCA + B leading-edge lollipop (top); contrast-definition row;
#        C merged DEPs|effect-size + D enrichment; E UpSet + F rank location.
#   A  Sample PCA (standard group palette, matches other figures) + PERMANOVA
#   B  leading-edge proteins, faceted/stacked by contrast, GO-CC coloured  [companion]
#   C  DEPs per contrast (nested p/FDR/Π) with effect-size violins to the right
#   D  pathway enrichment (total + mito, Up/Down, √-scaled)                [companion]
#   E  contrast overlap (UpSet; dodged Up/Down, singles first — YvO logic)
#   F  DEP rank location, with mito-protein counts                        [companion]

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(readr); library(readxl); library(ggplot2); library(ggtext)
  library(patchwork); library(cowplot); library(scales); library(vegan)
})

source(here::here("05_Figures", "shared", "config.R"))
source(here::here("04_Figures", "shared", "figure_supplement_helpers.R"))

BASE    <- fig05_base("F01_QC_overview")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels"); PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

CORE  <- c("CTLvMITO", "CTLvPHE", "PHEvPHE_MITO", "Interaction")
comb <- read_csv(P05$comb, show_col_types = FALSE); dep_df <- comb
fgsea_all <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"), show_col_types = FALSE)
dep_results <- setNames(lapply(CORE, \(c) as.data.frame(read_excel(P05$dep_xlsx, sheet = c))), CORE)

# ════════════════════════════════════════════════════════════════════════════
# Panel A — Sample PCA (standard group palette) + PERMANOVA in top headroom
# ════════════════════════════════════════════════════════════════════════════
dal_imp  <- readRDS(P05$imp_rds)
imp_mat  <- as.matrix(dal_imp$data)
imp_meta <- as_tibble(dal_imp$metadata); imp_meta$Group <- factor(imp_meta$Group, levels = H9C2_GROUP_LEVELS)
imp_mat  <- imp_mat[, imp_meta$Col_ID]
pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- as.data.frame(pca$x[, 1:2]) |> mutate(Col_ID = rownames(pca$x)) |>
  left_join(imp_meta |> dplyr::select(Col_ID, Group), by = "Col_ID") |>
  mutate(Group = factor(Group, levels = H9C2_GROUP_LEVELS))
set.seed(42); dist_mat <- dist(scale(t(imp_mat)))
perm <- adonis2(dist_mat ~ Group, data = imp_meta, permutations = 9999, by = "terms")
perm_R2 <- perm["Group", "R2"]; perm_p <- perm["Group", "Pr(>F)"]
pairs_to_test <- list(Transplant = c("Ctl", "Mito"), Disease = c("Ctl", "PHE"), Rescue = c("PHE", "PHE_Mito"))
pair_res <- bind_rows(lapply(names(pairs_to_test), function(nm) {
  pr <- pairs_to_test[[nm]]; keep <- imp_meta$Group %in% pr
  sm <- imp_mat[, keep]; smeta <- imp_meta[keep, ]; smeta$Group <- droplevels(smeta$Group)
  set.seed(42); r <- adonis2(dist(scale(t(sm))) ~ Group, data = smeta, permutations = 9999, by = "terms")
  tibble(role = nm, R2 = r$R2[1], p = r$`Pr(>F)`[1]) }))
write_csv(pair_res, file.path(DAT, "panel_A_pairwise_permanova.csv")); write_csv(pca_df, file.path(DAT, "panel_A_pca.csv"))

GRP_SHP <- c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)
GRP_LAB <- c(Ctl = "Ctl", Mito = "Mito", PHE = "PHE", PHE_Mito = "PHE+Mito")
# Dispersion homogeneity check — PERMANOVA assumes equal multivariate spread.
set.seed(42)
bd_p <- permutest(betadisper(dist_mat, imp_meta$Group), permutations = 999)$tab$`Pr(>F)`[1]
if (!is.na(bd_p) && bd_p < 0.05)
  warning("Heterogeneous group dispersions (betadisper p < 0.05) — interpret PERMANOVA cautiously")
fmt_perm <- function(role, r2, p) sprintf("%s: R²=%.3f, %s", role, r2, fmt_p(p))
perm_label <- paste(c("PERMANOVA",
  fmt_perm("Group", perm_R2, perm_p), fmt_perm("Transplant (Ctl|Mito)", pair_res$R2[1], pair_res$p[1]),
  fmt_perm("Disease (Ctl|PHE)", pair_res$R2[2], pair_res$p[2]),
  fmt_perm("Rescue (PHE|PHE+Mito)", pair_res$R2[3], pair_res$p[3]),
  sprintf("dispersion (betadisper): %s", fmt_p(bd_p))), collapse = "\n")
xr <- range(pca_df$PC1); yr <- range(pca_df$PC2); ytop <- yr[2] + 0.42 * diff(yr)
pA <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.08, level = 0.80, show.legend = FALSE) +
  stat_ellipse(level = 0.80, linewidth = 0.3, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.9) +
  annotate("label", x = xr[1] - 0.02 * diff(xr), y = ytop + 0.03 * diff(yr), label = perm_label,
           hjust = 0, vjust = 1, size = 1.7, color = "grey20", fontface = "bold", lineheight = 0.95,
           fill = alpha("white", 0.85), label.size = 0, label.padding = unit(0.15, "lines")) +
  scale_color_manual(values = H9C2_PAL_GROUP, labels = GRP_LAB, name = NULL,
                     guide = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  scale_fill_manual(values = H9C2_PAL_GROUP, guide = "none") +
  scale_shape_manual(values = GRP_SHP, labels = GRP_LAB, name = NULL, guide = guide_legend(ncol = 1)) +
  coord_cartesian(ylim = c(yr[1] - 0.02 * diff(yr), ytop + 0.03 * diff(yr))) +
  labs(title = "Sample PCA",
       subtitle = sprintf("n = %d, %s proteins (imputed)", nrow(imp_meta), format(nrow(imp_mat), big.mark = ",")),
       x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]), tag = "A") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, face = "bold.italic", color = "grey30"),
        legend.position = c(0.015, 0.58), legend.justification = c(0, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.key = element_blank(), legend.key.size = unit(2.8, "mm"),
        legend.text = element_text(size = FIG_LEGEND_TEXT + 0.5), legend.spacing.y = unit(0.3, "mm"),
        axis.title.y = element_text(face = "bold", size = 5, margin = margin(r = 1)),
        axis.text.y = element_text(margin = margin(r = 0)),
        plot.margin = margin(5, 2, 1, 1))

# ════════════════════════════════════════════════════════════════════════════
# Panel B — leading-edge lollipop (companion)
# ════════════════════════════════════════════════════════════════════════════
source(here::here("05_Figures", "F01_QC_overview", "a_script", "_panel_lollipop.R"))  # -> pLOL (tag B)

# ════════════════════════════════════════════════════════════════════════════
# Panel C — DEPs per contrast (nested p/FDR/Π) + effect-size violins to the right
# ════════════════════════════════════════════════════════════════════════════
all_genes <- unique(comb$gene[!is.na(comb$gene)]); n_total <- length(all_genes)
count_sig <- function(r) { pv <- if ("pi_score" %in% names(r)) r$pi_score else NA
  c(p = sum(!is.na(r$P.Value) & r$P.Value < 0.05), fdr = sum(!is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR),
    pi = if (is.numeric(pv)) sum(!is.na(pv) & pv < H9C2_PI_THRESH) else 0L) }
sig_counts <- sapply(CORE, \(c) count_sig(dep_results[[c]]))
CTR_LAB <- setNames(contrast_brief(CORE), CORE)
ctr_levels <- unname(CTR_LAB)   # facet order top-to-bottom: Transplant -> Interaction
frac_df <- bind_rows(lapply(CORE, \(ctr) { cc <- sig_counts[, ctr]
  tibble(contrast = CTR_LAB[ctr], threshold = c("p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), "Π < 0.05"),
         n = c(cc["p"], cc["fdr"], cc["pi"])) })) |>
  mutate(contrast = factor(contrast, levels = rev(ctr_levels)),
         threshold = factor(threshold, levels = c("p < 0.05", paste0("q < ", H9C2_FDR_EXPLOR), "Π < 0.05")),
         pct = 100 * n / n_total, fill_key = paste(contrast, threshold, sep = "___")) |> filter(n > 0)
SET_COLS <- setNames(unname(CONTRAST_COLORS[CORE]), unname(CTR_LAB))
FRAC_FILL <- c(); for (cn in names(SET_COLS)) { col <- unname(SET_COLS[cn])
  FRAC_FILL[paste(cn, "p < 0.05", sep = "___")] <- adjustcolor(col, alpha.f = 0.18)
  FRAC_FILL[paste(cn, paste0("q < ", H9C2_FDR_EXPLOR), sep = "___")] <- adjustcolor(col, alpha.f = 0.45)
  FRAC_FILL[paste(cn, "Π < 0.05", sep = "___")] <- col }
THRESH_LAB <- c("p < 0.05" = "p", "Π < 0.05" = "Π"); THRESH_LAB[paste0("q < ", H9C2_FDR_EXPLOR)] <- "FDR"
label_df <- frac_df |> arrange(contrast, threshold) |>
  mutate(next_pct = lead(pct, default = 0), seg = pct - next_pct, label_y = (next_pct + pct) / 2,
         label = THRESH_LAB[as.character(threshold)],
         text_col = if_else(threshold == "p < 0.05", "grey20", "white"), .by = contrast) |> filter(seg > 1.0)
# DEP bars: single compact figure (4 close bars, colour-coded by contrast,
# legible in-bar p/FDR/Π threshold labels). No numeric counts (decluttered).
panel_C_bg <- tibble(contrast = factor(unname(CTR_LAB), levels = rev(ctr_levels)),
                     fill = unname(CONTRAST_COLORS[CORE]))
pDEP <- ggplot(frac_df, aes(contrast, pct, fill = fill_key)) +
  geom_rect(data = panel_C_bg, aes(xmin = as.integer(contrast) - 0.5, xmax = as.integer(contrast) + 0.5,
            ymin = -Inf, ymax = Inf), fill = panel_C_bg$fill, alpha = 0.16, inherit.aes = FALSE) +
  geom_col(position = "identity", width = 0.86, color = "black", linewidth = 0.3) +
  geom_text(data = label_df, aes(contrast, label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5, size = scale_text(BASE_COUNT, 60) + 1.0, fontface = "bold") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(breaks = c(0, 5, 10, 20), limits = c(0, 22), expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(title = "DEP counts", subtitle = "n at p / FDR / Π (independent thresholds)",
       x = NULL, y = "% of proteome", tag = "C") +
  FIG_THEME + theme(legend.position = "none",
        plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, face = "italic", color = "grey40"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), plot.margin = margin(5, 0, 1, 2))
# Effect size: histograms with density overlay, faceted by contrast (aligns with
# panel B). Minimal stat per facet: median |log2FC|.
lfc_long <- bind_rows(lapply(CORE, \(c) tibble(contrast = CTR_LAB[c], logFC = dep_results[[c]]$logFC))) |>
  filter(!is.na(logFC), abs(logFC) <= 1) |> mutate(contrast = factor(contrast, levels = unname(CTR_LAB)))
lfc_stats <- lfc_long |> summarise(med_abs = median(abs(logFC)), .by = contrast) |>
  mutate(lab = sprintf("med|LFC| %.2f", med_abs))
hbw <- 2 / 44
hist_bg <- tibble(contrast = factor(unname(CTR_LAB), levels = unname(CTR_LAB)),
                  fill = unname(CONTRAST_COLORS[CORE]))
pHIST <- ggplot(lfc_long, aes(logFC)) +
  geom_rect(data = hist_bg, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = hist_bg$fill, alpha = 0.16, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_histogram(aes(fill = contrast), breaks = seq(-1, 1, by = hbw), color = "white", linewidth = 0.1, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * hbw), color = "grey20", linewidth = 0.4) +
  geom_text(data = lfc_stats, aes(x = -0.95, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = 1.4, size = scale_text(BASE_STAT, 50) + 0.5, fontface = "bold", color = "grey25") +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_fill_manual(values = setNames(unname(CONTRAST_COLORS[CORE]), unname(CTR_LAB)), guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1)) + scale_y_continuous(breaks = NULL) +
  labs(title = "Effect size", x = expression(bold(log[2]~FC)), y = NULL) +
  FIG_THEME + theme(strip.text = element_blank(), strip.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = FIG_AXIS_TEXT),
        panel.grid = element_blank(), panel.spacing.y = unit(1.5, "pt"),
        plot.margin = margin(5, 2, 1, 0))
pC <- pDEP + pHIST + plot_layout(widths = c(1, 0.62))

# ════════════════════════════════════════════════════════════════════════════
# Panel D — UpSet (dodged Up/Down, singles first; discordant computed, not barred)
# ════════════════════════════════════════════════════════════════════════════
sig_long <- bind_rows(lapply(CORE, function(c) { r <- dep_results[[c]]
  sig <- if ("pi_score" %in% names(r)) !is.na(r$pi_score) & r$pi_score < H9C2_PI_THRESH else rep(FALSE, nrow(r))
  tibble(uniprot_id = r$uniprot_id[sig], contrast = c, dir = ifelse(r$logFC[sig] > 0, "Up", "Down")) }))
all_sig <- unique(sig_long$uniprot_id)
mem_mat <- vapply(CORE, function(c) all_sig %in% sig_long$uniprot_id[sig_long$contrast == c], logical(length(all_sig)))
colnames(mem_mat) <- CORE
int_key <- apply(mem_mat, 1, function(r) paste(as.integer(r), collapse = ""))
consensus <- sapply(all_sig, function(g) { d <- sig_long$dir[sig_long$uniprot_id == g]
  if (all(d == "Up")) "Up" else if (all(d == "Down")) "Down" else "Mixed" })
inter_df <- tibble(uniprot_id = all_sig, int_key = int_key, direction = consensus,
                   degree = nchar(gsub("0", "", int_key)))
# FDR subset: a protein is an "FDR DEP" if adj.P.Val < threshold in ALL the
# contrasts of its membership (consistent with its Π-based intersection).
fdr_sets <- lapply(setNames(CORE, CORE), function(c) { r <- dep_results[[c]]
  r$uniprot_id[!is.na(r$adj.P.Val) & r$adj.P.Val < H9C2_FDR_EXPLOR] })
active_bits <- strsplit(inter_df$int_key, "")
inter_df$fdr_dep <- vapply(seq_len(nrow(inter_df)), function(i) {
  acts <- CORE[active_bits[[i]] == "1"]; uid <- inter_df$uniprot_id[i]
  length(acts) > 0 && all(vapply(acts, function(ct) uid %in% fdr_sets[[ct]], logical(1)))
}, logical(1))
write_csv(inter_df, file.path(DAT, "panel_D_upset_membership.csv"))
agg <- inter_df |> filter(int_key != strrep("0", length(CORE))) |>
  summarise(up_pi = sum(direction == "Up"), up_fdr = sum(direction == "Up" & fdr_dep),
            down_pi = sum(direction == "Down"), down_fdr = sum(direction == "Down" & fdr_dep),
            mixed = sum(direction == "Mixed"), .by = int_key) |>
  mutate(degree = nchar(gsub("0", "", int_key)), total = up_pi + down_pi + mixed)
# order: single-contrast first in contrast order (Transplant, Disease, Rescue,
# Interaction), then multi-contrast (overlap) intersections by ascending size.
single_keys <- vapply(seq_along(CORE), function(i) {
  k <- rep("0", length(CORE)); k[i] <- "1"; paste(k, collapse = "") }, character(1))
single_keys <- single_keys[single_keys %in% agg$int_key]
multi_keys <- (agg |> filter(degree > 1) |> arrange(total))$int_key
ord_keys <- head(c(single_keys, multi_keys), 10); n_int <- length(ord_keys)
agg <- agg |> filter(int_key %in% ord_keys) |>
  mutate(int_key = factor(int_key, levels = ord_keys), x = as.integer(int_key))
deg_by_x <- agg |> distinct(int_key, degree) |> arrange(int_key)
single_bg <- tibble(x = which(deg_by_x$degree == 1), key = ord_keys[deg_by_x$degree == 1]) |>
  mutate(fill = unname(CONTRAST_COLORS[CORE[ vapply(key, \(k) which(strsplit(k, "")[[1]] == "1")[1], integer(1)) ]]))
# Dodged Up/Down; within each, Π bar (light) with FDR subset (dark) overlaid from 0.
off <- 0.19; bw <- 0.30
rect_df <- bind_rows(
  agg |> transmute(x, dir = "Up",   layer = "pi",  val = up_pi),
  agg |> transmute(x, dir = "Up",   layer = "fdr", val = up_fdr),
  agg |> transmute(x, dir = "Down", layer = "pi",  val = down_pi),
  agg |> transmute(x, dir = "Down", layer = "fdr", val = down_fdr)) |>
  mutate(xc = x + ifelse(dir == "Up", -off, off),
         fill = dplyr::case_when(dir == "Up" & layer == "pi" ~ "#F4A582", dir == "Up" & layer == "fdr" ~ "#B2182B",
                                 dir == "Down" & layer == "pi" ~ "#92C5DE", TRUE ~ "#2166AC"))
tops <- agg |> transmute(x, Up = up_pi, Down = down_pi) |>
  tidyr::pivot_longer(c(Up, Down), names_to = "dir", values_to = "val") |>
  mutate(xc = x + ifelse(dir == "Up", -off, off))
zero_lab <- tops |> filter(val == 0)   # floating 0s above zero-hit direction bars
y_max <- max(c(agg$up_pi, agg$down_pi)) * 1.20
p_bars <- ggplot() +
  { if (nrow(single_bg) > 0) geom_rect(data = single_bg, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf),
            fill = single_bg$fill, alpha = 0.18, inherit.aes = FALSE) } +
  geom_rect(data = rect_df |> filter(layer == "pi"), aes(xmin = xc - bw/2, xmax = xc + bw/2, ymin = 0, ymax = val),
            fill = (rect_df |> filter(layer == "pi"))$fill, color = "black", linewidth = 0.2) +
  geom_rect(data = rect_df |> filter(layer == "fdr", val > 0), aes(xmin = xc - bw/2, xmax = xc + bw/2, ymin = 0, ymax = val),
            fill = (rect_df |> filter(layer == "fdr", val > 0))$fill, color = "black", linewidth = 0.2) +
  geom_text(data = tops |> filter(val > 0), aes(xc, val, label = val), vjust = -0.35,
            size = 2.0, fontface = "bold", color = "grey20") +
  { if (nrow(zero_lab) > 0) geom_text(data = zero_lab, aes(xc, 0, label = "0"), vjust = -0.3,
            size = 2.0, fontface = "bold", color = "grey55") } +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, y_max), breaks = breaks_pretty(4)) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5)) +
  labs(x = NULL, y = "Π-DEPs", title = "Contrast overlap (UpSet)", tag = "E") +
  FIG_THEME +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.title.y = element_text(size = FIG_AXIS_TEXT - 0.5), axis.text.y = element_text(size = FIG_AXIS_TEXT - 0.5),
        plot.margin = margin(5, 4, 0, 2))
# 4-item key (Up/Down x Π/FDR), top-right.
ukey <- tibble(y = 4:1, lab = c("Up (Π)", "Up (FDR)", "Down (Π)", "Down (FDR)"),
               fill = c("#F4A582", "#B2182B", "#92C5DE", "#2166AC"))
p_ukey <- ggplot(ukey) + geom_point(aes(0, y), shape = 22, size = 1.7, fill = ukey$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(0.2, y, label = lab), size = 1.4, hjust = 0, color = "grey20") +
  scale_x_continuous(limits = c(-0.1, 2.6)) + scale_y_continuous(limits = c(0.5, 4.5)) + theme_void()
p_bars <- p_bars + inset_element(p_ukey, left = 0.60, right = 0.99, top = 0.99, bottom = 0.64)
set_levels <- rev(CORE)
dot_df <- expand_grid(x = seq_len(n_int), set = factor(set_levels, levels = set_levels)) |>
  mutate(active = vapply(seq_len(n()), function(i)
    substr(ord_keys[x[i]], which(CORE == as.character(set[i])), which(CORE == as.character(set[i]))) == "1", logical(1)))
seg_df <- dot_df |> filter(active) |> summarise(y_lo = min(as.integer(set)), y_hi = max(as.integer(set)), .by = x) |> filter(y_lo < y_hi)
stripe <- tibble(set = factor(set_levels, levels = set_levels), fill = unname(CONTRAST_COLORS[set_levels]))
p_grid <- ggplot(dot_df, aes(x, set)) +
  geom_rect(data = stripe, aes(xmin = -Inf, xmax = Inf, ymin = as.integer(set) - 0.5, ymax = as.integer(set) + 0.5),
            fill = stripe$fill, alpha = 0.18, inherit.aes = FALSE) +
  { if (nrow(seg_df) > 0) geom_segment(data = seg_df, aes(x = x, xend = x, y = y_lo, yend = y_hi), inherit.aes = FALSE, color = "grey20", linewidth = 0.5) } +
  geom_point(aes(color = active), size = 1.5) +
  scale_color_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey80"), guide = "none") +
  scale_y_discrete(labels = setNames(contrast_brief(set_levels), set_levels)) +
  scale_x_continuous(expand = expansion(add = 0)) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5)) +
  labs(x = NULL, y = NULL) + FIG_THEME +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT - 1), panel.border = element_blank(),
        plot.margin = margin(0, 4, 4, 2))
pD <- (p_bars / p_grid) + plot_layout(heights = c(2.6, 1))

# Panels E (enrichment) and F (rank) via companions
source(here::here("05_Figures", "F01_QC_overview", "a_script", "_panel_F_enrichment.R"))  # -> pF (tag E)
source(here::here("05_Figures", "F01_QC_overview", "a_script", "_panel_G_rank.R"))         # -> pG (tag F)
p_enrich <- pF; p_rank <- pG

# ════════════════════════════════════════════════════════════════════════════
# Composite (178 x 230 mm)
# ════════════════════════════════════════════════════════════════════════════
COMP_W <- 178; COMP_H <- 178
# Combined A panel = PCA + leading-edge lollipop, aligned, sharing a key strip
# (contrast key left of the GO Slim annotation key). Lollipop facet shading maps
# to the contrast key.
pLOL <- pLOL + theme(legend.position = "none")
lol_ck <- tibble(name = contrast_brief(LOLLI_CORE), fill = unname(CONTRAST_COLORS[LOLLI_CORE]),
                 y = rev(seq_along(LOLLI_CORE)))
contrast_key <- ggplot(lol_ck) +
  annotate("text", x = 0, y = max(lol_ck$y) + 1.1, label = "Contrast (facet shade)",
           hjust = 0, size = 1.9, fontface = "bold", color = "grey20") +
  geom_point(aes(0.05, y), shape = 22, size = 2.2, fill = lol_ck$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(0.20, y, label = name), hjust = 0, size = 1.6, color = "grey15") +
  scale_x_continuous(limits = c(0, 2.3)) + scale_y_continuous(limits = c(0.5, max(lol_ck$y) + 1.8)) +
  theme_void()
present_cats <- intersect(CONSOLIDATED_PATHWAY_ORDER, as.character(unique(na.omit(lolli_long$compartment))))
lol_gk <- tibble(name = present_cats, fill = unname(CONSOLIDATED_COLORS[present_cats]),
                 y = rev(seq_along(present_cats)))
goslim_key <- ggplot(lol_gk) +
  annotate("text", x = 0, y = max(lol_gk$y) + 1.1, label = "GO Slim category (points)",
           hjust = 0, size = 1.9, fontface = "bold", color = "grey20") +
  geom_point(aes(0.05, y), shape = 22, size = 2.2, fill = lol_gk$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(0.20, y, label = name), hjust = 0, size = 1.6, color = "grey15") +
  scale_x_continuous(limits = c(0, 3.2)) + scale_y_continuous(limits = c(0.5, max(lol_gk$y) + 1.8)) +
  theme_void()
keys_AB <- (contrast_key | goslim_key) + plot_layout(widths = c(1, 1.35))
pAB <- ((pA | pLOL) / wrap_elements(full = keys_AB)) + plot_layout(heights = c(1, 0.24))
row1 <- wrap_elements(full = pAB) + wrap_elements(full = pC) + plot_layout(widths = c(2.05, 1))
# Contrast key band: uniform accent-shaded boxes (swatch + 2-line name/equation),
# evenly spaced with equal gaps, no overlap.
def_df <- tibble(name = contrast_brief(CORE), def = unname(CONTRAST_MATH_BRIEF[CORE]),
                 col = unname(CONTRAST_COLORS[CORE]),
                 x0 = (seq_along(CORE) - 1) * 0.25 + 0.008)
def_bw <- 0.234
pDEF <- ggplot(def_df) +
  geom_rect(aes(xmin = x0, xmax = x0 + def_bw, ymin = 0.12, ymax = 0.88),
            fill = alpha(def_df$col, 0.20), color = def_df$col, linewidth = 0.35) +
  geom_point(aes(x = x0 + 0.014, y = 0.5), shape = 22, size = 2.1,
             fill = def_df$col, color = "grey30", stroke = 0.3) +
  geom_text(aes(x = x0 + 0.034, y = 0.5, label = sprintf("%s\n%s", name, def)),
            hjust = 0, vjust = 0.5, size = 1.55, lineheight = 0.95, color = "grey15") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  theme_void() + theme(plot.margin = margin(1, 2, 1, 2))
row2 <- wrap_elements(full = p_enrich) + wrap_elements(full = pD) + wrap_elements(full = p_rank) +
  plot_layout(widths = c(0.85, 1.55, 0.72))
composite <- (row1 / wrap_elements(full = pDEF) / row2) +
  plot_layout(heights = c(1.05, 0.10, 1.15)) +
  plot_annotation(
    title = "H9c2 mito-transplant proteome — QC & overview",
    subtitle = sprintf("%s proteins × %d wells (imputed) | 2×2 PHE×Mito (n=6/group; Interaction underpowered) | PERMANOVA Group R²=%.3f, %s | Π = P.Value^|logFC| < %.2f",
                       format(nrow(imp_mat), big.mark = ","), nrow(imp_meta), perm_R2, fmt_p(perm_p), H9C2_PI_THRESH),
    theme = theme(plot.title = element_text(face = "bold", size = composite_text_sizes(COMP_H)$title),
                  plot.subtitle = element_text(face = "italic", size = composite_text_sizes(COMP_H)$subtitle, color = "grey30"))) &
  theme(plot.tag = element_text(face = "bold", size = composite_text_sizes(COMP_H)$tag))

ggsave(file.path(RPT_PDF, "MAIN_F01_composite.pdf"), composite, width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F01_composite.png"), composite, width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  sheet_specs = list(
    list(name = "contrast_map",       df = fig05_contrast_table()),
    list(name = "panel_A_pca",        df = as.data.frame(pca_df)),
    list(name = "panel_A_permanova",  df = as.data.frame(pair_res)),
    list(name = "panel_C_dep_counts", df = as.data.frame(frac_df)),
    list(name = "panel_D_upset",      df = as.data.frame(inter_df)),
    list(name = "panel_E_enrichment", df = as.data.frame(read_csv(file.path(DAT, "panel_F_enrichment_sig.csv"), show_col_types = FALSE))),
    list(name = "panel_F_rank_counts",df = as.data.frame(read_csv(file.path(DAT, "panel_G_rank_counts.csv"), show_col_types = FALSE))),
    list(name = "panel_B_lollipop",   df = as.data.frame(read_csv(file.path(DAT, "panel_lollipop_leadingedge.csv"), show_col_types = FALSE)))))

message(sprintf("F01 MAIN composite (%dx%d mm) | PERMANOVA p=%.4f | %d unique Π DEPs", COMP_W, COMP_H, perm_p, length(all_sig)))
