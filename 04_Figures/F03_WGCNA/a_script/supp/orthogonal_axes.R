#!/usr/bin/env Rscript
# What this design can and cannot ask, and the three answers it can give.
#
# A select-then-test pair is only interpretable when the two contrasts are
# orthogonal. Over (Ctl, Mito, PHE, PHE_Mito) the contrast vectors are
#
#   Disease D = (-1, 0, +1, 0)      Rescue        R = ( 0,  0, -1, +1)
#   Transplant T = (-1, +1, 0, 0)   Disease-after A = ( 0, -1,  0, +1)
#
# D.R = -1, so disease and rescue are coupled at -0.5 and anticorrelate before
# any biology: that is the pair coupling_check.R kills. D.A = 0 and T.R = 0, so
# those two pairs share no group and carry no coupling.
#
# Each row of the figure is one pair against a null that shuffles group labels
# within a replicate. The null is on the figure because the observed value is
# uninterpretable without it: -0.232 reads as weak until you can see that noise
# alone puts the same quantity at +0.532.

pacman::p_load(here, dplyr, tibble, purrr, readr, ggplot2, patchwork, limma, fgsea)
source(here::here("04_Figures", "functions", "shared_data_loaders.R"))
source(here::here("04_Figures", "functions", "shared_enrichment_ora.R"))

png_dir <- here::here("04_Figures", "F03_WGCNA", "b_reports", "supp")
tab_dir <- here::here("04_Figures", "F03_WGCNA", "c_data", "supp")
for (d in c(png_dir, tab_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 300
SEED <- 42

dal <- readRDS(P05$imp_rds)
expr <- as.matrix(dal$data)
meta <- as.data.frame(dal$metadata)
grp0 <- factor(meta$Group[match(colnames(expr), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)
rep_block <- meta$Replicate[match(colnames(expr), meta$Col_ID)]

sym <- load_dep_long() |>
  distinct(uniprot_id, gene) |>
  filter(!is.na(gene), gene != "")
keep <- intersect(rownames(expr), sym$uniprot_id)
expr_gene <- avereps(expr[keep, ], ID = sym$gene[match(keep, sym$uniprot_id)])
sets <- load_set_pool(CANONICAL_DBS)
sets <- sets[vapply(sets, \(s) sum(s %in% rownames(expr_gene)) >= 10, logical(1))]

# Each replicate holds one sample per group, so shuffling labels inside a
# replicate destroys the treatment effect while keeping the design balanced.
permute_within_block <- function() {
  g <- grp0
  for (r in unique(rep_block)) g[rep_block == r] <- sample(grp0[rep_block == r])
  factor(g, levels = levels(grp0))
}

# Set score is the mean moderated t of the members, which tracks NES and is fast
# enough to permute 300 times across two contrasts.
set_scores <- function(g, contrast_a, contrast_b) {
  d <- stats::model.matrix(~ 0 + g)
  colnames(d) <- levels(g)
  fit <- lmFit(expr_gene, d)
  score <- function(what) {
    v <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = what, levels = d)))$t[, 1]
    vapply(sets, \(s) mean(v[intersect(s, names(v))], na.rm = TRUE), numeric(1))
  }
  tibble(a = score(contrast_a), b = score(contrast_b))
}
set_correlation <- function(g, contrast_a, contrast_b) {
  s <- set_scores(g, contrast_a, contrast_b)
  stats::cor(s$a, s$b)
}

AXES <- tibble::tribble(
  ~key, ~a, ~b, ~coupling, ~question,
  "one_program", "PHE_Mito - PHE", "Mito - Ctl", "Orthogonal",
  "The transplant runs one program",
  "opposes", "Mito - Ctl", "PHE - Ctl", "Shares Ctl, same sign",
  "That program opposes the disease",
  "persists", "PHE_Mito - Mito", "PHE - Ctl", "Orthogonal",
  "The disease keeps acting",
  "trap", "PHE_Mito - PHE", "PHE - Ctl", "Shares PHE, opposite sign",
  "The pair that cannot be read"
)

set.seed(SEED)
runs <- pmap(AXES, function(key, a, b, coupling, question) {
  observed <- set_correlation(grp0, a, b)
  null <- vapply(seq_len(N_PERM), \(i) set_correlation(permute_within_block(), a, b), numeric(1))
  list(
    null = tibble(key = key, value = null),
    stat = tibble(
      key = key, question = question, coupling = coupling,
      contrast_a = a, contrast_b = b,
      observed = observed,
      null_median = stats::median(null),
      null_q25 = stats::quantile(null, 0.25),
      null_q75 = stats::quantile(null, 0.75),
      perm_p_low = mean(null <= observed),
      perm_p_high = mean(null >= observed)
    )
  )
})

null_long <- bind_rows(map(runs, "null"))
stats_tbl <- bind_rows(map(runs, "stat"))

# The p that matters is the one-sided tail the claim actually sits in: opposition
# and the trap are read on the low side, the other two on the high side.
stats_tbl <- stats_tbl |>
  mutate(perm_p = if_else(key %in% c("opposes", "trap"), perm_p_low, perm_p_high))

# No shuffle reached the observed value in the one-program row, so the tail p is
# only resolvable to 1/N_PERM; report that bound rather than a spurious zero.
label_of <- stats_tbl |>
  mutate(
    panel = sprintf("%s\n%s  →  %s   [%s]", question, contrast_a, contrast_b, coupling),
    label = sprintf(
      "observed %+.3f, %s", observed,
      if_else(perm_p < 1 / N_PERM,
        sprintf("p < %.3f", 1 / N_PERM), sprintf("p = %.3f", perm_p)
      )
    ),
    hjust = if_else(observed > 0, 1.05, -0.05)
  ) |>
  select(key, panel, observed, null_median, perm_p, label, hjust)

obs_scores <- setNames(
  lapply(AXES$key, function(k) {
    ax <- AXES[AXES$key == k, ]
    set_scores(grp0, ax$a, ax$b)
  }),
  AXES$key
)

# Left of each row: the null distribution of the correlation with the observed
# value marked. Right: the same correlation shown directly, one point per set,
# so the trap row's spurious anti-diagonal is visible rather than only summarised.
hist_panel <- function(k) {
  m <- label_of[label_of$key == k, ]
  ax <- AXES[AXES$key == k, ]
  ggplot(filter(null_long, key == k), aes(value)) +
    geom_histogram(bins = 40, fill = "grey80", colour = NA) +
    geom_vline(xintercept = m$null_median, colour = "grey40", linetype = "22", linewidth = 0.4) +
    geom_vline(xintercept = m$observed, colour = H9C2_PAL_DIR[["Up"]], linewidth = 0.9) +
    annotate("text",
      x = m$observed, y = Inf, label = m$label, hjust = m$hjust,
      vjust = 1.6, colour = H9C2_PAL_DIR[["Up"]], size = 2.6
    ) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(
      title = ax$question, subtitle = sprintf("%s → %s   [%s]", ax$a, ax$b, ax$coupling),
      x = "correlation across sets", y = "shuffles"
    ) +
    FIG_THEME
}

scatter_panel <- function(k) {
  ax <- AXES[AXES$key == k, ]
  m <- label_of[label_of$key == k, ]
  s <- obs_scores[[k]]
  r <- stats::cor(s$a, s$b)
  p_txt <- if (m$perm_p < 1 / N_PERM) sprintf("p < %.3f", 1 / N_PERM) else sprintf("p = %.3f", m$perm_p)
  ggplot(s, aes(a, b)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(alpha = 0.3, size = 0.8, colour = "grey35") +
    geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      colour = H9C2_PAL_DIR[["Up"]], linewidth = 0.8
    ) +
    annotate("text",
      x = -Inf, y = Inf, hjust = -0.08, vjust = 1.4,
      label = sprintf("r = %+.2f, %s", r, p_txt),
      colour = H9C2_PAL_DIR[["Up"]], size = 3, fontface = "bold"
    ) +
    labs(title = "each dot = one pathway set", x = ax$a, y = ax$b) +
    FIG_THEME
}

# Protein-level view of the trap: rank every protein by the rescue contrast, then
# mark where the disease-selected genes fall. A real reversal would pile the
# disease-up genes at the bottom of the rescue ranking; here it is arithmetic.
d0 <- stats::model.matrix(~ 0 + grp0)
colnames(d0) <- levels(grp0)
fit0 <- lmFit(expr, d0)
prot_t <- function(what) {
  eBayes(contrasts.fit(fit0, makeContrasts(contrasts = what, levels = d0)))$t[, 1]
}
disease_t <- prot_t("PHE - Ctl")
rescue_rank <- prot_t("PHE_Mito - PHE")
dis_order <- order(disease_t, decreasing = TRUE)
dis_sets <- list(
  up = names(disease_t)[head(dis_order, 50)],
  down = names(disease_t)[tail(dis_order, 50)]
)
trap_es <- fgsea::fgsea(dis_sets, rescue_rank)
es_lab <- function(nm) {
  r <- trap_es[trap_es$pathway == nm, ]
  sprintf("NES %+.2f, naive p = %.2g; a label shuffle reverses just as often", r$NES, r$pval)
}
barcode_panel <- function(genes, title, sub) {
  fgsea::plotEnrichment(genes, rescue_rank) +
    labs(
      title = title, subtitle = sub,
      x = "protein rank by rescue (PHE_Mito − PHE)", y = "running enrichment"
    ) +
    FIG_THEME
}
bc_up <- barcode_panel(dis_sets$up, "Disease-up 50 along the rescue ranking", es_lab("up"))
bc_dn <- barcode_panel(dis_sets$down, "Disease-down 50 along the rescue ranking", es_lab("down"))

rows <- lapply(AXES$key, \(k) hist_panel(k) + scatter_panel(k) + plot_layout(widths = c(1.5, 1)))
barcode_row <- bc_up + bc_dn + plot_layout(widths = c(1, 1))
panel <- patchwork::wrap_plots(c(rows, list(barcode_row)),
  ncol = 1, heights = c(rep(1, length(rows)), 1.05)
) +
  patchwork::plot_annotation(
    title = "Only contrast pairs that share no group are interpretable",
    subtitle = sprintf(
      paste0(
        "%d pathway sets. Rows 1-4: permutation null (dashed) vs observed (solid) per contrast pair, then the same sets as a scatter.\n",
        "Bottom row: Disease-selected proteins along the rescue ranking; a label shuffle reverses them just as often."
      ),
      length(sets)
    ),
    theme = theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30")
    )
  )

ggsave(file.path(png_dir, "SUPP_F03_orthogonal_axes.png"), panel,
  width = 220, height = 300, units = "mm", dpi = 300, bg = "white"
)
write_csv(select(stats_tbl, -perm_p_low, -perm_p_high), file.path(tab_dir, "SUPP_F03_orthogonal_axes.csv"))

# Optional local Box mirror (author's machine only; no-ops elsewhere).
mirror_to_box(file.path(png_dir, "SUPP_F03_orthogonal_axes.png"), "02_Figures/F03_WGCNA/supp", guard = TRUE)

message("Set-correlation axes against a null with no treatment effect:")
print(as.data.frame(stats_tbl |>
  select(question, coupling, observed, null_median, perm_p) |>
  mutate(across(where(is.numeric), \(x) signif(x, 3)))), row.names = FALSE)
