# F04 Rescue — trajectory clustering of the disease (PHE) signature
# Soft (fuzzy c-means) clustering of each disease-signature protein's standardized
# profile across the rescue axis Ctl -> PHE -> PHE_Mito, coloured by the per-protein
# RESCUE class (engine band on the actual disease/rescue logFC), per-cluster ORA.
#   Engine = e1071::cmeans (fuzzy c-means Mfuzz wraps; Futschik & Carlisle 2005).
#   Descriptive only -- inference stays with the DEP contrasts / fry.
library(here)
setwd(here::here())                          # Mito root via .git (no .Rproj)
source("04_Figures/shared/config.R")        # sources style.R (palettes, theme, helpers)
source("04_Figures/shared/pathway_utils.R")
library(tidyverse)
library(e1071)
library(ggrepel)
library(patchwork)

RPT_PNG <- "04_Figures/F04_rescue/b_reports/main/png/panels"
DAT     <- "04_Figures/F04_rescue/c_data/panel_traj"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
set.seed(42)

ORDER    <- c("Ctl", "PHE", "PHE_Mito")          # control -> disease -> rescued
AXISLAB  <- c(Ctl = "Control", PHE = "PHE", PHE_Mito = "PHE+Mito")
K        <- 6
M_FUZZ   <- 1.25
PHI_BAND <- 0.25
CLASS_COLORS <- c(`Toward control` = "#1B7837", Unrescued = "#878787",
                  `Away from control` = "#D6604D")

# ── Data: rescue-axis group means + per-protein rescue class ─────────────────
wide <- load_combined_wide()
dal  <- readRDS(P05$imp_rds)
stopifnot(all(ORDER %in% dal$metadata$group))

gmean <- sapply(ORDER, function(g)
  rowMeans(dal$data[, dal$metadata$Col_ID[dal$metadata$group == g], drop = FALSE]))
rownames(gmean) <- rownames(dal$data)

# responsive set = Pi<0.05 in the disease (CTLvPHE) OR rescue (PHEvPHE_MITO) contrast.
# (PHE's disease effect is weak, so anchor on the responsive union, not disease alone.)
sig_ids <- wide |>
  filter(pi_score_CTLvPHE < 0.05 | pi_score_PHEvPHE_MITO < 0.05) |>
  distinct(uniprot_id, gene)

gmean <- gmean[rownames(gmean) %in% sig_ids$uniprot_id, , drop = FALSE]
z <- t(apply(gmean, 1, function(r) (r - mean(r)) / sd(r)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
message(sprintf("  Rescue trajectory clustering: %d responsive proteins (PHE OR rescue)", nrow(z)))

# Robust rescue class from the z-profile (no division): is PHE_Mito closer to
# control than PHE was? -> the transplant moved it back toward Ctl.
MARGIN <- 0.25
zc <- as_tibble(z, rownames = "uniprot_id") |>
  mutate(d_phe = abs(PHE - Ctl), d_mito = abs(PHE_Mito - Ctl),
         rescue_class = case_when(
           abs(PHE_Mito - PHE) < MARGIN ~ "Unrescued",          # transplant barely moved it
           d_mito < d_phe - MARGIN      ~ "Toward control",     # moved back toward Ctl
           d_mito > d_phe + MARGIN      ~ "Away from control",  # pushed further from Ctl
           TRUE                          ~ "Unrescued")) |>
  left_join(sig_ids, by = "uniprot_id")

# ── Fuzzy c-means ────────────────────────────────────────────────────────────
fcm <- cmeans(z, centers = K, m = M_FUZZ, iter.max = 200)
assign <- tibble(uniprot_id = rownames(z), cluster = fcm$cluster,
                 membership = apply(fcm$membership, 1, max)) |>
  left_join(select(zc, uniprot_id, gene, rescue_class), by = "uniprot_id")
write_csv(assign, file.path(DAT, "rescue_trajectory_clusters.csv"))

cent <- as_tibble(fcm$centers, .name_repair = "minimal")
names(cent) <- ORDER
cent <- cent |>
  mutate(cluster = row_number(),
         pattern = case_when(abs(PHE_Mito - PHE) < 0.4         ~ "PHE-locked",
                             abs(PHE_Mito - Ctl) < abs(PHE - Ctl) ~ "toward control",
                             TRUE                               ~ "away from control"))
clab <- cent |>
  mutate(n = as.integer(table(factor(assign$cluster, levels = seq_len(K)))[cluster]),
         lab = sprintf("Cluster %d · %s (n=%d)", cluster, pattern, n))
class_n <- assign |> count(rescue_class)
rescue_dir <- mean(assign$rescue_class == "Toward control", na.rm = TRUE)

# ── Per-cluster ORA (rat) ────────────────────────────────────────────────────
universe <- unique(wide$gene)
pw <- build_pathway_collection(species = "Rattus norvegicus", min_size = 15,
                               max_size = 500, include_goslim = FALSE, exclude_variants = TRUE)
top_path <- map_dfr(seq_len(K), function(k) {
  g <- assign$gene[assign$cluster == k]
  if (length(g) < 5) return(tibble(cluster = k, pathway_label = NA_character_, padj = NA_real_))
  r <- tryCatch(run_ora_deduplicated(genes = g, universe = universe, pathways = pw,
                                     jaccard_cutoff = 0.5, min_size = 15, max_size = 500,
                                     padj_cutoff = 1), error = function(e) tibble())
  if (nrow(r) == 0) return(tibble(cluster = k, pathway_label = NA_character_, padj = NA_real_))
  r |> arrange(padj) |> slice_head(n = 1) |>
    transmute(cluster = k, pathway_label = clean_pathway_name(pathway), padj)
})
clab <- clab |> left_join(top_path, by = "cluster") |>
  mutate(fdr_txt = ifelse(coalesce(padj, 1) < 0.001, "<0.001", sprintf("%.3f", coalesce(padj, 1))),
         path_lab = ifelse(is.na(pathway_label), "",
                           paste0(pathway_label, " (FDR ", fdr_txt, ")")))

# ── Trajectory plot ──────────────────────────────────────────────────────────
long <- as_tibble(z, rownames = "uniprot_id") |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  left_join(select(assign, uniprot_id, cluster, membership, rescue_class), by = "uniprot_id") |>
  mutate(cond = factor(cond, levels = ORDER))
cent_long <- cent |> select(cluster, all_of(ORDER)) |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  mutate(cond = factor(cond, levels = ORDER))
facet_lab <- setNames(clab$lab, clab$cluster)

p_traj <- ggplot(long, aes(cond, zval, group = uniprot_id)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_line(aes(color = rescue_class, alpha = membership), linewidth = 0.25) +
  geom_line(data = cent_long, aes(cond, zval, group = cluster),
            inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_text(data = clab, aes(x = 2, y = Inf, label = path_lab), inherit.aes = FALSE,
            vjust = 1.4, size = 2.4, fontface = "italic", color = "grey25", lineheight = 0.9) +
  facet_wrap(~ cluster, ncol = 3, labeller = labeller(cluster = facet_lab)) +
  scale_color_manual(values = CLASS_COLORS, name = "Rescue class") +
  scale_alpha_continuous(range = c(0.05, 0.5), guide = "none") +
  scale_x_discrete(labels = AXISLAB) +
  labs(x = NULL, y = "standardized abundance (z)") +
  FIG_THEME +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.5, face = "bold"))

composite <- p_traj +
  plot_annotation(
    title = "Mito Rescue: trajectory clustering of the rescue-responsive signature",
    subtitle = sprintf("Fuzzy c-means k=%d | %d responsive proteins | %.1f%% toward control | Toward %d · Unrescued %d · Away %d",
                       K, nrow(z), 100 * rescue_dir,
                       sum(class_n$n[class_n$rescue_class == "Toward control"]),
                       sum(class_n$n[class_n$rescue_class == "Unrescued"]),
                       sum(class_n$n[class_n$rescue_class == "Away from control"])),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "MAIN_panel_traj_composite.png"), composite,
       width = 200, height = 150, units = "mm", dpi = 300)
message(sprintf("Mito rescue trajectory done | %.1f%% rescued", 100 * rescue_dir))
