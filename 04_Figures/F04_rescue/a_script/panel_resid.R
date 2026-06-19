# F04 Rescue — residual volcano: PHE_Mito vs Ctl ("what the transplant didn't fix")
# The rescue residual is not one of the 5 primary DEP contrasts, so it is computed
# here with the SAME non-imputed limma model (~0+group+(1|Replicate)), then shown
# as a volcano + fgsea. Pairs with the rescue trajectory (panel_traj, "what moved
# back toward control").
library(here)
source(here::here("04_Figures/shared/config.R"))
source(here::here("04_Figures/shared/pathway_utils.R"))
suppressPackageStartupMessages({ library(proteoDA)
  library(dplyr); library(ggplot2); library(readr); library(tibble)
  library(fgsea); library(ggrepel); library(patchwork) })

RPT_PNG <- here::here("04_Figures/F04_rescue/b_reports/main/png/panels")
DAT     <- here::here("04_Figures/F04_rescue/c_data/panel_resid")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# ── Residual contrast PHE_Mito − Ctl via the primary (non-imputed) model ─────
dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
dal <- add_design(dal, "~ 0 + group + (1 | Replicate)")
dal <- add_contrasts(dal, contrasts_vector = c("RescueResid = PHE_Mito - Ctl"))
dal <- fit_limma_model(dal)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
res <- as_tibble(dal$results$RescueResid, rownames = "uniprot_id") |>
  left_join(ann, by = "uniprot_id") |>
  mutate(pi = P.Value ^ abs(logFC))

v <- res |>
  transmute(gene, logFC, t, P = pmax(P.Value, 1e-300), pi) |>
  filter(!is.na(logFC), !is.na(P)) |>
  mutate(neglogP = -log10(P), sig = pi < 0.05,
         dir = case_when(!sig ~ "NS", logFC > 0 ~ "Up", TRUE ~ "Down"))
n_sig <- sum(v$sig)
write_csv(v, file.path(DAT, "resid_volcano.csv"))
lab <- v |> filter(sig) |> arrange(pi) |> slice_head(n = 18)

p_volc <- ggplot(v, aes(logFC, neglogP)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(data = filter(v, !sig), color = "grey80", size = 0.5, alpha = 0.4) +
  geom_point(data = filter(v, sig), aes(color = dir), size = 0.9, alpha = 0.8) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.4, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.2, segment.color = "grey50",
                  min.segment.length = 0, box.padding = 0.3, seed = 42) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]), Down = unname(DIR_COLORS["Down"])),
                     name = NULL, labels = c(Down = "lower than control", Up = "higher than control")) +
  labs(x = expression(log[2]*FC~"(PHE+Mito − Control)"), y = expression(-log[10]~italic(P))) +
  FIG_THEME + theme(legend.position = "top", panel.grid.minor = element_blank())

# ── fgsea on the residual t-statistic (rat) ──────────────────────────────────
rr <- res |> select(gene, t) |> filter(is.finite(t)) |> distinct(gene, .keep_all = TRUE)
ranks <- setNames(rr$t, rr$gene)
pw <- build_harmonized_collection(min_size = 15, max_size = 500)
set.seed(42)
gs <- fgsea(pathways = pw, stats = sort(ranks), eps = 0)
top <- gs |> as_tibble() |> filter(padj < 0.10) |>
  mutate(pathway_label = clean_pathway_name(pathway)) |>
  group_by(side = ifelse(NES > 0, "Higher than control", "Lower than control")) |>
  arrange(padj) |> slice_head(n = 8) |> ungroup() |>
  arrange(NES) |> mutate(pathway_label = factor(pathway_label, levels = pathway_label))
write_csv(select(top, pathway, pathway_label, NES, pval, padj, side), file.path(DAT, "resid_fgsea_top.csv"))

p_path <- ggplot(top, aes(NES, pathway_label, color = NES > 0)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_label), linewidth = 0.5) +
  geom_point(aes(size = -log10(padj))) +
  scale_color_manual(values = c(`TRUE` = unname(DIR_COLORS["Up"]), `FALSE` = unname(DIR_COLORS["Down"])),
                     guide = "none") +
  scale_size_continuous(name = expression(-log[10]~FDR), range = c(1.5, 4)) +
  labs(x = "NES (residual)", y = NULL, title = "Pathways still off after transplant") +
  FIG_THEME +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 7), legend.position = "right",
        panel.grid.minor = element_blank())

composite <- (p_volc | p_path) + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = "Mito Rescue: residual contrast (what transplant did NOT return to control)",
    subtitle = sprintf("RescueResid = PHE+Mito − Control | %d proteins still differ (Π<0.05) | fgsea FDR<0.10", n_sig),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "MAIN_panel_resid_composite.png"), composite,
       width = 220, height = 120, units = "mm", dpi = 300)
message(sprintf("Mito residual volcano done | %d residual DEPs", n_sig))
