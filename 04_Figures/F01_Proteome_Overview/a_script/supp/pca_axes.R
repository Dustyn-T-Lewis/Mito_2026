#!/usr/bin/env Rscript
# Supplement to F01 panel A: exploratory characterisation of the global proteome
# PCA axes. Two sample PCAs (factorial + pairwise PERMANOVA, betadisper), the
# subcellular make-up of the loadings (descriptive), and fgsea on the signed
# loadings ranked by FDR. Descriptive only — the confirmatory enrichment lives in
# F02 (DE contrasts). PC1 is variance-defined, not a designed contrast.
#
# Methods: PERMANOVA + PERMDISP (Anderson 2001; Anderson 2006); loading-based set
# analysis on the full ranking, not a top-N cut (Yamamoto 2014, BMC Bioinformatics
# 15:51); GSEA metric-agnostic (Subramanian 2005); EnrichmentMap dedup (Merico
# 2010; Reimand 2019); MitoCarta3 compartments (Rath 2021).

pacman::p_load(here, dplyr, tibble, tidyr, readr, stringr, forcats, ggplot2, patchwork, vegan, msigdbr, ggtext, scales, fgsea, ggrepel)
fns <- here("04_Figures", "functions")
source(file.path(fns, "shared_theme_palettes.R"))
source(file.path(fns, "shared_data_loaders.R"))
source(file.path(fns, "shared_enrichment_ora.R"))
source(here("04_Figures", "functions", "shared_gene_set_helpers.R"))
source(here("04_Figures", "functions", "f01_pca_stats.R"))
source(file.path(fns, "shared_nes_bars.R"))
set.seed(42)

FIG_DIR <- here("04_Figures", "F01_Proteome_Overview", "b_reports", "supp")
TAB_DIR <- here("04_Figures", "F01_Proteome_Overview", "c_data", "supp")
for (d in c(FIG_DIR, TAB_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(P05$imp_rds)
mat <- as.matrix(dal$data)
meta <- as_tibble(dal$metadata)
mat <- mat[, meta$Col_ID]
gene <- dal$annotation$gene[match(rownames(mat), dal$annotation$uniprot_id)]
group <- factor(meta$Group, levels = H9C2_GROUP_LEVELS)
label <- factor(GROUP_LABELS[as.character(group)], levels = unname(GROUP_LABELS[H9C2_GROUP_LEVELS]))
transplant <- factor(if_else(group %in% c("Mito", "PHE_Mito"), "Transplant", "No transplant"), levels = c("No transplant", "Transplant"))
phe <- factor(if_else(group %in% c("PHE", "PHE_Mito"), "PHE", "no"))
mito <- factor(if_else(group %in% c("Mito", "PHE_Mito"), "Mito", "no"))

pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
ve <- 100 * summary(pca)$importance[2, 1:2]
D <- dist(scale(t(mat)))
fac <- adonis2(D ~ phe * mito, by = "terms", permutations = 9999)
tx_perm <- adonis2(D ~ transplant, by = "terms", permutations = 9999)
pair_stats <- pca_group_stats(mat, group, list(Disease = c("Ctl", "PHE"), Transplant = c("Ctl", "Mito"), Rescue = c("PHE", "PHE_Mito")))
bd_p <- pair_stats$betadisper_p

permanova_tab <- bind_rows(
  tibble(
    model = "factorial", term = c("Phe (disease)", "Mito (transplant)", "Phe x Mito (rescue)"),
    R2 = round(fac[c("phe", "mito", "phe:mito"), "R2"], 4), p = fac[c("phe", "mito", "phe:mito"), "Pr(>F)"]
  ),
  tibble(model = "pairwise", term = pair_stats$pairwise$role, R2 = round(pair_stats$pairwise$R2, 4), p = pair_stats$pairwise$p, q = round(pair_stats$pairwise$q, 4)),
  tibble(model = "transplant", term = "Transplant vs none", R2 = round(tx_perm["transplant", "R2"], 4), p = tx_perm["transplant", "Pr(>F)"]),
  tibble(model = "dispersion", term = "betadisper", R2 = NA_real_, p = round(bd_p, 4))
)
write_csv(permanova_tab, file.path(TAB_DIR, "SUPP_F01_permanova.csv"))

## Compartments — descriptive annotation, not an enrichment test.
gene_up <- toupper(gene)
mito_genes <- toupper(readRDS(here("04_Figures", "shared", "mitocarta3_rat.rds"))$MITOCARTA_ALL)
cc <- msigdbr(species = "Rattus norvegicus", collection = "C5", subcollection = "GO:CC")
bucket <- function(re) unique(toupper(cc$gene_symbol[grepl(re, cc$gs_name)]))
# Only the compartments that move along an axis, plus Nucleus for scale. Golgi, ER
# and Cytosol sit within a point or two of each other at both ends of PC1 and add
# rows without adding signal; Extracellular is under 1% of the proteome. Everything
# not listed falls into Other, which is reported rather than hidden.
compartments <- list(
  Mitochondrion = union(mito_genes, bucket("MITOCHONDR")),
  Cytoskeleton = bucket("CYTOSKELET|ACTIN|MICROTUBULE|INTERMEDIATE_FILAMENT|MYOFIBRIL|SARCOMERE"),
  Nucleus = bucket("NUCLEUS|NUCLEOPLASM|NUCLEAR_|CHROMOSOM|CHROMATIN|NUCLEOLUS")
)
comp_levels <- c(names(compartments), "Other")
comp_v <- rep("Other", length(gene_up))
for (cn in rev(names(compartments))) comp_v[gene_up %in% compartments[[cn]]] <- cn

grp_pal <- setNames(GROUP_COLORS[H9C2_GROUP_LEVELS], GROUP_LABELS[H9C2_GROUP_LEVELS])
shp <- setNames(c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)[H9C2_GROUP_LEVELS], GROUP_LABELS[H9C2_GROUP_LEVELS])
tx_pal <- c("No transplant" = "#D6604D", "Transplant" = "#4393C3")
comp_pal <- c(
  Mitochondrion = "#009988", Cytoskeleton = "#CC3311", Nucleus = "#8FBAD6",
  Other = "grey85"
)
# The dot strip sits directly above the compartment bars and makes the same
# argument, so it borrows their colours: the end transplanted cells occupy is the
# mitochondrial end. Under tx_pal the strip spends red on "no transplant" while
# the bars beneath spend red on cytoskeleton, and one figure asks red to mean two
# things.
end_pal <- c(
  "No transplant" = unname(comp_pal["Cytoskeleton"]),
  "Transplant" = unname(comp_pal["Mitochondrion"])
)
fp <- function(p) if (p < 1e-3) "p<0.001" else sprintf("p=%.3f", p)
bold_if <- function(txt, p) if (p < 0.05) paste0("**", txt, "**") else txt

df <- tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = label, transplant = transplant)
grp_cent <- df |>
  group_by(label) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
ctl_c <- filter(grp_cent, label == "CTL")
grp_arrow <- grp_cent |>
  filter(label != "CTL") |>
  mutate(x0 = ctl_c$PC1, y0 = ctl_c$PC2)
tx_cent <- df |>
  group_by(transplant) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
tx_arrow <- tibble(x0 = tx_cent$PC1[1], y0 = tx_cent$PC2[1], PC1 = tx_cent$PC1[2], PC2 = tx_cent$PC2[2])

rbox <- function(txt) {
  geom_richtext(
    data = tibble(lab = txt), aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 2.7, lineheight = 1.15, fill = alpha("white", 0.85),
    label.colour = "grey75", label.padding = unit(c(3, 4, 3, 4), "pt"), nudge_x = 4, nudge_y = -4
  )
}
arrow_style <- arrow(length = unit(2.4, "mm"), type = "closed")
grp_txt <- paste(
  bold_if(sprintf("Phe (disease): %.0f%%, %s", 100 * fac["phe", "R2"], fp(fac["phe", "Pr(>F)"])), fac["phe", "Pr(>F)"]),
  bold_if(sprintf("Mito (transplant): %.0f%%, %s", 100 * fac["mito", "R2"], fp(fac["mito", "Pr(>F)"])), fac["mito", "Pr(>F)"]),
  bold_if(sprintf("Phe×Mito (rescue): %.0f%%, %s", 100 * fac["phe:mito", "R2"], fp(fac["phe:mito", "Pr(>F)"])), fac["phe:mito", "Pr(>F)"]),
  sep = "<br>"
)
tx_txt <- paste(
  bold_if(sprintf("Transplant: R²=%.2f, %s", tx_perm["transplant", "R2"], fp(tx_perm["transplant", "Pr(>F)"])), tx_perm["transplant", "Pr(>F)"]),
  sprintf("4-group dispersion (betadisper): %s", fp(bd_p)),
  sep = "<br>"
)

pca_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold"))
p_group <- ggplot(df, aes(PC1, PC2)) +
  stat_ellipse(aes(group = label, fill = label, colour = label), geom = "polygon", alpha = 0.18, linewidth = 0.7, level = 0.8) +
  scale_fill_manual(values = grp_pal, guide = "none") +
  geom_segment(data = grp_arrow, aes(x = x0, y = y0, xend = PC1, yend = PC2, colour = label), arrow = arrow_style, linewidth = 0.9, show.legend = FALSE) +
  geom_point(aes(colour = label, shape = label), size = 2.6, alpha = 0.95) +
  scale_colour_manual(values = grp_pal, name = NULL) +
  scale_shape_manual(values = shp, name = NULL) +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "By experimental group", subtitle = "80% ellipses; arrows = centroid shift from CTL; factorial PERMANOVA") +
  rbox(grp_txt) +
  coord_cartesian(clip = "off") +
  pca_theme
p_tx <- ggplot(df, aes(PC1, PC2)) +
  stat_ellipse(aes(group = transplant, colour = transplant, fill = transplant), geom = "polygon", alpha = 0.14, linewidth = 0.7, level = 0.8) +
  geom_segment(data = tx_arrow, aes(x = x0, y = y0, xend = PC1, yend = PC2), arrow = arrow_style, linewidth = 0.9, colour = "grey20") +
  geom_point(aes(colour = transplant, shape = label), size = 2.6, alpha = 0.95) +
  scale_colour_manual(values = tx_pal, name = NULL) +
  scale_fill_manual(values = tx_pal, guide = "none") +
  scale_shape_manual(values = shp, guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = "Transplant vs no-transplant (mito contrast)", subtitle = "points coloured by transplant, shaped by group; arrow = centroid shift") +
  rbox(tx_txt) +
  coord_cartesian(clip = "off") +
  pca_theme

# Why PC1 and not another axis. Drawing PC1 alone invites the reply that it was
# chosen because it looked right, so every axis worth plotting is scored on both
# things PC1 is claimed to be: it separates transplanted cells, and its ends
# differ in mitochondrial content. Neither property alone is special - PC5 and PC8
# are mitochondrial without tracking the transplant, PC6 tracks the transplant
# without being mitochondrial. Only PC1 does both, and that conjunction is the
# argument.
N_PC <- 8
ve_all <- 100 * summary(pca)$importance[2, ]

# A principal component's sign is arbitrary, so each axis is oriented to put the
# transplanted samples on the positive side before its ends are compared. That
# makes "the transplant-ward end" mean the same thing on every axis.
pc_scores <- pca$x[, seq_len(N_PC), drop = FALSE]
orient <- sign(vapply(
  seq_len(N_PC),
  \(k) mean(pc_scores[transplant == "Transplant", k]) - mean(pc_scores[transplant == "No transplant", k]),
  numeric(1)
))
orient[orient == 0] <- 1

end_share <- function(k, member) {
  l <- pca$rotation[, k] * orient[k]
  top <- order(abs(l), decreasing = TRUE)[seq_len(500)]
  mean(member[top][l[top] >= 0]) - mean(member[top][l[top] < 0])
}
is_mito <- gene_up %in% compartments$Mitochondrion
is_cyto <- gene_up %in% compartments$Cytoskeleton

pc_map <- tibble(
  pc = paste0("PC", seq_len(N_PC)),
  variance = ve_all[seq_len(N_PC)],
  separation = vapply(
    seq_len(N_PC),
    \(k) abs(stats::t.test(pc_scores[, k] ~ transplant)$statistic), numeric(1)
  ),
  sep_p = vapply(
    seq_len(N_PC), \(k) stats::t.test(pc_scores[, k] ~ transplant)$p.value, numeric(1)
  ),
  mito_gap = vapply(seq_len(N_PC), \(k) 100 * end_share(k, is_mito), numeric(1)),
  cyto_gap = vapply(seq_len(N_PC), \(k) 100 * end_share(k, is_cyto), numeric(1))
) |>
  mutate(is_axis = pc == "PC1")
write_csv(
  mutate(pc_map, across(where(is.numeric), \(x) round(x, 3))),
  file.path(TAB_DIR, "SUPP_F01_pc_map.csv")
)

cc_universe <- list(`top 500` = function(l) order(abs(l), decreasing = TRUE)[seq_len(500)], whole = function(l) seq_along(l))
cc_df <- bind_rows(lapply(names(cc_universe), function(un) {
  bind_rows(lapply(c("PC1", "PC2"), function(pc) {
    l <- pca$rotation[, pc]
    idx <- cc_universe[[un]](l)
    tibble(universe = un, pc = pc, side = if_else(l[idx] >= 0, "+", "−"), comp = comp_v[idx]) |>
      count(universe, pc, side, comp, .drop = FALSE) |>
      group_by(universe, pc, side) |>
      mutate(frac = n / sum(n)) |>
      ungroup()
  }))
}))
write_csv(mutate(cc_df, frac = round(frac, 4)), file.path(TAB_DIR, "SUPP_F01_compartment_shares.csv"))
# Generic end labels. Only PC1 separates transplant from control, so naming the
# ends "transplant"/"control" on every panel would assert that of PC2 as well.
end_lab <- c(`−` = "\u2212 end", `+` = "+ end")
# Two compartments and one axis. Nucleus sits at 26% against 26% and Golgi, ER
# and cytosol within a point or two, so plotting them spends four panels to show
# eight flat bars. The full four-category table stays in c_data for the question
# "what about everything else".
cc_plot_df <- cc_df |>
  filter(pc == "PC1", universe == "top 500", comp %in% c("Mitochondrion", "Cytoskeleton")) |>
  mutate(
    comp = factor(comp, levels = c("Mitochondrion", "Cytoskeleton")),
    end = factor(
      if_else(side == "+", "Transplant end", "Control end"),
      levels = c("Control end", "Transplant end")
    )
  )
# "Transplant end" is a label, not a finding: a component's sign is arbitrary, so
# orient[1] flips PC1 to put transplanted samples on the positive side and the end
# takes their name. That is only worth saying because they do separate along it,
# and this strip is where a reader can see that rather than take it on trust. It
# sits above the bars on the same axis direction, so the end a bar describes is
# the end the samples are sitting at.
pc1_strip_df <- tibble(
  pc1 = pca$x[, 1] * orient[1],
  transplant = transplant
)
# The ends follow the samples, not a symmetric axis: one control sits far out at
# -120 and a mirrored limit would park "Transplant end" in empty space.
pc1_rng <- range(pc1_strip_df$pc1)
pc1_pad <- 0.04 * diff(pc1_rng)
pc1_strip <- ggplot(pc1_strip_df, aes(pc1, transplant)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4, linetype = "22") +
  geom_point(
    aes(fill = transplant),
    shape = 21, size = 3.4, stroke = 0.3, colour = "grey25", alpha = 0.85
  ) +
  annotate("text",
    x = pc1_rng[1], y = 2.62, label = "Control end", hjust = 0,
    fontface = "bold", size = 4.2, colour = "grey25"
  ) +
  annotate("text",
    x = pc1_rng[2], y = 2.62, label = "Transplant end", hjust = 1,
    fontface = "bold", size = 4.2, colour = "grey25"
  ) +
  scale_fill_manual(values = end_pal, guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.1))) +
  coord_cartesian(xlim = c(pc1_rng[1] - pc1_pad, pc1_rng[2] + pc1_pad)) +
  labs(
    title = "PC1 score by transplant status",
    subtitle = "One dot per sample. PC1's sign is arbitrary, so it is flipped to put transplanted cells on the right.",
    x = "PC1 score", y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey30")
  )

# A slope, not paired bars: the finding is that the two compartments trade places
# along PC1, and one line crossing another states that. Faceted bars made a reader
# hold four numbers across two panels and do the subtraction themselves.
cc_end_lab <- filter(cc_plot_df, end == "Transplant end")
cc_panel <- ggplot(cc_plot_df, aes(end, frac, group = comp, colour = comp)) +
  geom_line(linewidth = 1.7) +
  geom_point(size = 4.6) +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * frac)),
    vjust = -1.15, size = 5, fontface = "bold", show.legend = FALSE
  ) +
  geom_text(
    data = cc_end_lab, aes(label = comp),
    hjust = -0.16, size = 4.6, fontface = "bold", show.legend = FALSE
  ) +
  scale_y_continuous(
    labels = percent, limits = c(0, 0.34),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(expand = expansion(add = c(0.45, 1.4))) +
  scale_colour_manual(values = comp_pal[c("Mitochondrion", "Cytoskeleton")], guide = "none") +
  labs(
    title = "Compartment membership at each end of PC1",
    subtitle = "Share of the 500 strongest PC1 loadings at each end",
    x = NULL, y = "share of proteins"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey30", lineheight = 1.2)
  )

## fgsea on signed loadings, full ranking, all 5 DBs, dedup. Bars ranked by FDR.
dep <- tibble(gene = gene, t_PC1 = pca$rotation[, 1], t_PC2 = pca$rotation[, 2])
sets <- readRDS(here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"))
pw_flat <- build_harmonized_collection()
fg_all <- bind_rows(lapply(CANONICAL_DBS, \(db) run_fgsea_cache(dep, sets[[db]], db, c("PC1", "PC2")))) |>
  filter(size >= 10, size <= 500, !pathway %in% MITO_DROP_SETS, !grepl(DISEASE_VIRAL_RE, pathway, ignore.case = TRUE)) |>
  mutate(label = clean_display_label(pathway), database = factor(database, levels = names(ORA_DB_COLORS)))
kept <- bind_rows(lapply(c("PC1", "PC2"), \(pc) deduplicate_enrichment(filter(fg_all, contrast == pc), pw_flat, cutoff = 0.375, cross_db = TRUE) |> mutate(contrast = pc)))
fg_all <- fg_all |> mutate(dedup_kept = paste(contrast, pathway) %in% paste(kept$contrast, kept$pathway))
write_csv(
  fg_all |> transmute(contrast, database, pathway, label, NES = round(NES, 3), pval = signif(pval, 3), padj = signif(padj, 3), size, dedup_kept, leadingEdge),
  file.path(TAB_DIR, "SUPP_F01_fgsea_loadings.csv")
)
write_csv(
  tibble(uniprot_id = rownames(mat), gene = gene, PC1_loading = round(pca$rotation[, 1], 5), PC2_loading = round(pca$rotation[, 2], 5), compartment = comp_v),
  file.path(TAB_DIR, "SUPP_F01_loadings.csv")
)

# Same renderer as the deck aim bars: name inside the bar, database on the axis in
# its colour, q at the tip, significant rows white and bold.
pc_bars <- function(pc, n_each) {
  d <- kept |> filter(contrast == pc)
  top <- bind_rows(
    d |> filter(NES > 0) |> slice_min(padj, n = n_each, with_ties = FALSE),
    d |> filter(NES < 0) |> slice_min(padj, n = n_each, with_ties = FALSE)
  )
  top |>
    mutate(
      display = str_trunc(label, 44),
      direction = factor(if_else(NES > 0, "Up", "Down"), levels = c("Up", "Down")),
      sig = factor(padj < H9C2_FDR_STD, levels = c(TRUE, FALSE))
    ) |>
    group_by(display) |>
    mutate(label = if (n_distinct(pathway) > 1) paste0(display, " (", database, ")") else display) |>
    ungroup() |>
    mutate(label = fct_reorder(label, NES))
}

# One PNG per panel: the deck walks these one at a time, so a baked-in composite
# would force four arguments onto one slide.
saved_panels <- character()
save_panel <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, sprintf("SUPP_F01_%s.png", name)), p,
    width = w, height = h, units = "mm", dpi = 300, bg = "white", limitsize = FALSE
  )
  saved_panels <<- c(saved_panels, name)
}
save_panel(p_group, "pca_group", 180, 150)
save_panel(p_tx, "pca_transplant", 180, 150)
save_panel(
  pc1_strip / cc_panel + patchwork::plot_layout(heights = c(1, 1.5)),
  "compartments", 190, 150
)
# PC1 alone, built to the same recipe as the F03 aim bars rather than to
# FGSEA_THEME: the deck shows this panel four slides before those, and a viewer
# reads two panels with different text weights and bar heights as two different
# kinds of evidence. Eight rows at this panel size, matching MAX_SETS there, so a
# bar is the same thickness on every slide that carries bars. The two-panel
# version above keeps its own scale for the supplement; PC2 clears nothing, so
# the deck says it rather than showing it.
# Twice the 92 mm the deck aim-bar script gives a one-sided panel, because this one runs
# both ways from zero: at 92 mm the same NES span is drawn in half the
# millimetres and the names and q-values run off the ends.
DECK_PANEL_W <- 184
DECK_SLIDE_H <- 118
DECK_THEME <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(face = "bold", size = 13, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 9, colour = "grey30", margin = margin(b = 4)),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9, colour = "grey15"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )
pc1_deck <- pc_bars("PC1", 4)
save_panel(
  nes_bar_panel(
    pc1_deck, sprintf("PC1 (%.0f%% of variance)", ve[["PC1"]]),
    "fgsea on the signed PC1 loadings, 5-DB lens, EnrichmentMap-deduplicated; + = transplant end, − = control end",
    xlab = "NES", base_theme = DECK_THEME, q_size = 2.4
  ) & theme(legend.position = "bottom", legend.justification = "center"),
  "fgsea_pc1", DECK_PANEL_W, DECK_SLIDE_H
)

message(sprintf(
  "F01 PCA-axes supplement: %d panels (%s) | %d fgsea sets (%d kept after dedup) -> %s",
  length(saved_panels), paste(saved_panels, collapse = ", "),
  nrow(fg_all), nrow(kept), FIG_DIR
))

# Optional local Box mirror (author's machine only; no-ops elsewhere).
mirror_to_box(
  list.files(FIG_DIR, "^SUPP_F01_(pca_group|pca_transplant|compartments|fgsea_pc1)\\.png$", full.names = TRUE),
  "02_Figures/F01_Proteome_Overview/supp",
  guard = TRUE
)
