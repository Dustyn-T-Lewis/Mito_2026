#!/usr/bin/env Rscript
# F06 (4/4) — subunit co-regulation supplement (no fractionation required).
# A ProteomeHD-flavoured view (Kustatscher et al. 2019, Nat Biotechnol,
# doi:10.1038/s41587-019-0298-5): within each complex, derive a module axis (PC1)
# and measure how tightly each subunit tracks it. Subunits that decouple from their
# complex are candidates for independent (stoichiometric) behaviour. Exploratory at
# this n — descriptive correlations only, no inference.

setwd(rprojroot::find_rstudio_root_file())
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr)
  library(limma); library(ggplot2); library(patchwork); library(cowplot)
})
source("04_Figures/shared/style.R")
set.seed(42)

DAT      <- "05_Figures/F06_complex_mito/c_data"
SUPP_PDF <- "05_Figures/F06_complex_mito/b_reports/supp/pdf"
SUPP_PNG <- "05_Figures/F06_complex_mito/b_reports/supp/png"
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

da   <- readRDS(here::here("02_Imputation", "c_data", "01_DAList_imputed.rds"))
meta <- readRDS(file.path(DAT, "analysis_meta.rds"))
gene <- da$annotation$gene[match(rownames(da$data), da$annotation$uniprot_id)]
keep <- !is.na(gene) & nzchar(gene)
mat  <- limma::avereps(da$data[keep, ], ID = gene[keep])

COMPLEX_LV <- c("Complex I", "Complex III", "Complex IV", "Complex V", "Mitoribosome")
members    <- meta$complex_members[COMPLEX_LV]

# Per complex: module axis (PC1 of z-scored subunits), each subunit's correlation to
# it (decoupling), and overall cohesion (mean pairwise Spearman r).
subunit_rows <- list(); complex_rows <- list()
for (cx in COMPLEX_LV) {
  g <- intersect(members[[cx]], rownames(mat))
  if (length(g) < 3) next
  sub <- mat[g, , drop = FALSE]
  z   <- t(scale(t(sub)))                       # gene-wise z across samples
  pc  <- prcomp(t(z), center = TRUE, scale. = FALSE)
  module <- pc$x[, 1]
  if (cor(module, colMeans(z)) < 0) module <- -module   # orient toward mean abundance
  r_to_module <- apply(z, 1, function(x) cor(x, module, method = "spearman"))
  cm  <- cor(t(z), method = "spearman")
  cohesion <- mean(cm[upper.tri(cm)])
  subunit_rows[[cx]] <- tibble(complex = cx, gene = g, r_to_module = unname(r_to_module))
  complex_rows[[cx]] <- tibble(complex = cx, n_subunits = length(g),
                               mean_pairwise_r = cohesion,
                               pc1_var_explained = pc$sdev[1]^2 / sum(pc$sdev^2))
}
subunit_cov <- bind_rows(subunit_rows) |> mutate(complex = factor(complex, levels = COMPLEX_LV))
complex_cov <- bind_rows(complex_rows) |> mutate(complex = factor(complex, levels = COMPLEX_LV))
write_csv(subunit_cov, file.path(DAT, "covariation_subunit.csv"))
write_csv(complex_cov, file.path(DAT, "covariation_complex.csv"))

# Cohesion CI (Fisher z, n = 24 samples) for the bar panel.
complex_cov <- complex_cov |>
  rowwise() |>
  mutate(ci = list(fisher_z_ci(mean_pairwise_r, n = ncol(mat)))) |>
  mutate(lo = ci[1], hi = ci[2]) |> ungroup() |> select(-ci)

# Empirically, within-complex co-regulation is weak in these 24 samples (cohesion near
# zero, PC1 explains <50%), so covariation-based stoichiometry inference is underpowered
# at this n — shown honestly as distributions, not thresholded "decoupled" gene calls.
pA <- ggplot(complex_cov, aes(complex, mean_pairwise_r, fill = complex)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.3, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.3, colour = "grey25") +
  geom_text(aes(label = sprintf("n=%d | PC1 %.0f%%", n_subunits, 100 * pc1_var_explained),
                y = pmax(hi, 0) + 0.03), vjust = 0, size = scale_text(BASE_STAT, 85)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(limits = c(min(0, min(complex_cov$lo)), max(0.4, max(complex_cov$hi) + 0.08)),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Within-complex co-regulation (cohesion)",
       subtitle = "mean pairwise Spearman r across 24 samples (±95% Fisher-z CI) | near-zero ⇒ weak co-regulation at this n",
       x = NULL, y = "mean pairwise r") +
  FIG_THEME + theme(axis.text.x = element_text(face = "bold", size = FIG_AXIS_TEXT))

med <- subunit_cov |> group_by(complex) |> summarise(m = median(r_to_module), .groups = "drop")
pB <- ggplot(subunit_cov, aes(complex, r_to_module, fill = complex)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_violin(alpha = 0.35, colour = "grey45", linewidth = 0.25, scale = "width", width = 0.8) +
  geom_jitter(width = 0.14, height = 0, size = 0.9, alpha = 0.7, colour = "grey25") +
  geom_point(data = med, aes(complex, m), shape = 95, size = 6, colour = "#B2182B", inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(limits = c(min(-0.8, min(subunit_cov$r_to_module)), 1.0)) +
  labs(title = "Subunit correlation to complex module (PC1)",
       subtitle = "distribution per complex (red bar = median) | broad, near-zero spread ⇒ no reliable subunit decoupling signal at n=24",
       x = NULL, y = "r to complex module") +
  FIG_THEME + theme(axis.text.x = element_text(face = "bold", size = FIG_AXIS_TEXT))

COMP_W <- 178; COMP_H <- 130
fig <- (pA / pB) + plot_layout(heights = c(1, 1.25))
txt <- composite_text_sizes(COMP_H)
comp <- ggdraw(fig) +
  draw_label("A", x = 0.005, y = 0.992, size = txt$tag + 4, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.005, y = 0.555, size = txt$tag + 4, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(SUPP_PDF, "SUPP_F06_covariation.pdf"), comp, width = COMP_W, height = COMP_H, units = "mm", device = pdf_device, limitsize = FALSE)
ggsave(file.path(SUPP_PNG, "SUPP_F06_covariation.png"), comp, width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)
message(sprintf("F06 covariation supplement done | %d subunits across %d complexes; cohesion range %.2f–%.2f (weak ⇒ covariation underpowered at n=%d)",
                nrow(subunit_cov), nrow(complex_cov),
                min(complex_cov$mean_pairwise_r), max(complex_cov$mean_pairwise_r), ncol(mat)))
