#!/usr/bin/env Rscript
# Pilot: is the PHE proteome shift reversed by transplant? Four panels build the
# claim from set-level down to per-protein: signature reversal (fry on the disease
# sets and a priori cardiac-stress pathways, tested Disease -> Recovery), a PCA
# trajectory showing each group's distance to control, the per-protein return with
# a mito / direction key, and the formal interaction cap.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "F01-F03_style_palettes_theme.R"))
source(file.path(fns, "F01-F03_data_paths_and_loaders.R"))
source(file.path(fns, "F01-F03_supplementary_workbook.R"))
source(file.path(fns, "F01-F03_composite_layout.R"))
pacman::p_load(limma, fgsea, dplyr, tidyr, tibble, purrr, patchwork, scales)

BASE <- here::here("04_Figures", "test", "pilot_reversal")
for (f in list.files(file.path(BASE, "a_script", "panels"), full.names = TRUE)) source(f)
RPT <- file.path(BASE, "b_reports")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

STRESS_SETS <- c(
  "HALLMARK_MTORC1_SIGNALING", "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
  "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

comb <- load_combined_wide()
mito_genes <- toupper(readRDS(here::here("04_Figures", "shared", "mitocarta3_rat.rds"))$MITOCARTA_ALL)

dal <- readRDS(P05$imp_rds)
expr <- as.matrix(dal$data)
meta <- as.data.frame(dal$metadata)
grp <- factor(meta$Group[match(colnames(expr), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)
rep_block <- meta$Replicate[match(colnames(expr), meta$Col_ID)]
genes <- toupper(dal$annotation$gene[match(rownames(expr), dal$annotation$uniprot_id)])
design <- stats::model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
CTR <- c(Disease = "PHE - Ctl", Recovery = "PHE_Mito - Ctl")
cm <- makeContrasts(contrasts = unname(CTR), levels = design)
colnames(cm) <- names(CTR)

# --- set-level reversal: disease sets + a priori cardiac-stress pathways ---
idx_sig <- list(
  `disease-up` = which(rownames(expr) %in% comb$uniprot_id[comb$sig_pi_CTLvPHE == 1]),
  `disease-down` = which(rownames(expr) %in% comb$uniprot_id[comb$sig_pi_CTLvPHE == -1])
)
hall <- readRDS(here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"))$Hallmark
idx_ext <- lapply(hall[STRESS_SETS], \(g) which(genes %in% toupper(g)))
idx_ext <- idx_ext[vapply(idx_ext, \(i) length(i) >= 5L, logical(1))]

set_rho <- duplicateCorrelation(expr, design, block = rep_block)$consensus

# Recovery (Phe+Mito − Ctl) is derived, not a canonical DE contrast, so compute it
# locally on the imputed matrix and join it in; the reported DE keeps the clean five.
rec_fit <- eBayes(contrasts.fit(
  lmFit(expr, design, block = rep_block, correlation = set_rho),
  makeContrasts("PHE_Mito - Ctl", levels = design)
))
comb <- topTable(rec_fit, number = Inf, sort.by = "none") |>
  as_tibble(rownames = "uniprot_id") |>
  transmute(
    uniprot_id,
    logFC_CTLvPHE_MITO = logFC, t_CTLvPHE_MITO = t, P.Value_CTLvPHE_MITO = P.Value,
    pi_score_CTLvPHE_MITO = P.Value^abs(logFC)
  ) |>
  right_join(comb, by = "uniprot_id")

fdr_col <- function(df, rn) if ("FDR" %in% colnames(df)) df[rn, "FDR"] else df[rn, "PValue"]
run_fry <- function(idx, type) {
  bind_rows(lapply(names(CTR), function(cn) {
    fr <- fry(expr, idx, design, contrast = cm[, cn], block = rep_block, correlation = set_rho)
    tibble(set = rownames(fr), set_type = type, contrast = cn, dir = fr$Direction, fry_fdr = fdr_col(fr, rownames(fr)))
  }))
}
pretty_set <- function(s) ifelse(grepl("^HALLMARK_", s), tolower(gsub("_", " ", sub("^HALLMARK_", "", s))), s)
set_rev <- bind_rows(run_fry(idx_sig, "disease signature"), run_fry(idx_ext, "cardiac-stress (external)")) |>
  mutate(signed = ifelse(dir == "Up", 1, -1) * -log10(fry_fdr), label = pretty_set(set))
ord <- set_rev |>
  filter(contrast == "Disease") |>
  arrange(set_type, signed) |>
  distinct(label) |>
  pull(label)
set_rev$label <- factor(set_rev$label, ord)

# --- proteome geometry: PCA trajectory + full-space distance to control ---
pc <- stats::prcomp(t(expr), scale. = FALSE)
var_exp <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2), 1)
scores <- as_tibble(pc$x[, 1:2], rownames = "sample") |>
  mutate(group = factor(grp[match(sample, colnames(expr))], H9C2_GROUP_LEVELS))
centroids <- scores |>
  group_by(group) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

ccol <- function(cols) rowMeans(expr[, cols, drop = FALSE])
ctl <- ccol(grp == "Ctl")
to_ctl <- function(v) sqrt(sum((v - ctl)^2))
d_phe <- to_ctl(ccol(grp == "PHE"))
d_res <- to_ctl(ccol(grp == "PHE_Mito"))
obs <- d_res / d_phe
stressed <- which(grp %in% c("PHE", "PHE_Mito"))
n_res <- sum(grp == "PHE_Mito")
set.seed(42)
null <- vapply(seq_len(2000), function(i) {
  r <- sample(stressed, n_res)
  to_ctl(ccol(r)) / to_ctl(ccol(setdiff(stressed, r)))
}, numeric(1))
dist_stat <- list(ratio = obs, p = (1 + sum(null <= obs)) / 2001, d_phe = d_phe, d_res = d_res)

# --- per-protein return with the mito / direction key ---
prot <- comb |>
  filter(sig_pi_CTLvPHE != 0) |>
  transmute(
    uniprot_id, gene,
    disease = logFC_CTLvPHE, recovery = logFC_CTLvPHE_MITO,
    direction = ifelse(logFC_CTLvPHE > 0, "disease-up", "disease-down"),
    lens = ifelse(toupper(gene) %in% mito_genes, "MitoCarta", "other"),
    normalized = pi_score_CTLvPHE_MITO >= H9C2_PI_THRESH
  )

# --- fry mountains: a signature set (from set_old) tested with fry over the proteome
# ranked by another contrast (rank_old); leading-edge proteins (broad Disease p<0.05)
# feed top-5 ORA bars, n.s. shaded lighter. ORA labels deduplicated by cleaned name. ---
pw_flat <- unlist(unname(readRDS(here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"))), recursive = FALSE)
universe <- unique(comb$gene[!is.na(comb$gene)])
run_ora <- function(g) {
  g <- unique(g[!is.na(g)])
  if (length(g) < 5) {
    return(tibble())
  }
  fora(pathways = pw_flat, genes = g, universe = universe, minSize = 10, maxSize = 500) |>
    as_tibble() |>
    arrange(padj) |>
    mutate(.key = .pretty_pw(pathway)) |>
    distinct(.key, .keep_all = TRUE) |>
    select(-.key) |>
    slice_head(n = 5)
}
make_mountain <- function(cfg) {
  set_sig <- comb[[paste0("sig_pi_", cfg$set_old)]]
  up_ids <- comb$uniprot_id[set_sig == 1]
  dn_ids <- comb$uniprot_id[set_sig == -1]
  cmx <- makeContrasts(contrasts = cfg$rank_form, levels = design)
  fry1 <- function(ids) {
    idx <- which(rownames(expr) %in% ids)
    fr <- fry(expr, idx, design, contrast = cmx[, 1], block = rep_block, correlation = set_rho)
    tibble(dir = fr$Direction[1], fdr = fdr_col(fr, rownames(fr))[1], n = length(idx))
  }
  fry_up <- fry1(up_ids)
  fry_dn <- fry1(dn_ids)
  tr <- comb |>
    mutate(t_rank_stat = .data[[paste0("t_", cfg$rank_old)]]) |>
    filter(!is.na(t_rank_stat)) |>
    transmute(gene, t_rank_stat, in_up = uniprot_id %in% up_ids, in_down = uniprot_id %in% dn_ids) |>
    arrange(desc(t_rank_stat)) |>
    mutate(rank = row_number(), es_up = running_es(t_rank_stat, in_up), es_down = running_es(t_rank_stat, in_down))
  n_all <- nrow(tr)
  bp <- comb[[paste0("P.Value_", cfg$set_old)]]
  blfc <- comb[[paste0("logFC_", cfg$set_old)]]
  brt <- comb[[paste0("t_", cfg$rank_old)]]
  side <- function(dir) if (dir == "Down") brt < 0 else brt > 0
  ora_up <- run_ora(comb$gene[which(bp < 0.05 & blfc > 0 & side(fry_up$dir))])
  ora_dn <- run_ora(comb$gene[which(bp < 0.05 & blfc < 0 & side(fry_dn$dir))])
  cfg$rank_label <- sprintf(cfg$rank_label, n_all)
  m <- build_mountain(tr, fry_up, fry_dn, ora_up, ora_dn, n_all, cfg)
  ggplot2::ggsave(file.path(RPT, paste0("PILOT_fry_", cfg$id, ".png")), m,
    width = 220, height = 150, units = "mm", dpi = 300, limitsize = FALSE
  )
  tibble(mountain = cfg$id, up_dir = fry_up$dir, up_fdr = fry_up$fdr, dn_dir = fry_dn$dir, dn_fdr = fry_dn$fdr)
}

MOUNTAINS <- list(
  list(
    id = "rescue_v_disease", rank_old = "PHEvPHE_MITO", rank_form = "PHE_Mito - PHE", set_old = "CTLvPHE",
    rank_label = "Rank (Rescue t, n = %d)", title = "Disease signature across the Rescue axis",
    subtitle = "rank = Rescue (Phe+Mito − Phe); set = Π disease; shares PHE, so partly circular",
    up_title = "Disease-up set", dn_title = "Disease-down set",
    ora_up_title = "reversed disease-up · ORA", ora_dn_title = "reversed disease-down · ORA"
  ),
  list(
    id = "recovery_v_disease", rank_old = "CTLvPHE_MITO", rank_form = "PHE_Mito - Ctl", set_old = "CTLvPHE",
    rank_label = "Rank (Recovery t, n = %d)", title = "Disease signature vs control after transplant",
    subtitle = "rank = Recovery (Phe+Mito − Ctl); set = Π disease; shares Ctl, confound-free",
    up_title = "Disease-up set", dn_title = "Disease-down set",
    ora_up_title = "residual disease-up · ORA", ora_dn_title = "residual disease-down · ORA"
  ),
  list(
    id = "disease_v_rescue", rank_old = "CTLvPHE", rank_form = "PHE - Ctl", set_old = "PHEvPHE_MITO",
    rank_label = "Rank (Disease t, n = %d)", title = "Rescue signature across the disease axis (inverse)",
    subtitle = "rank = Disease (Phe − Ctl); set = Π rescue; the reciprocal test",
    up_title = "Rescue-up set", dn_title = "Rescue-down set",
    ora_up_title = "rescue-up · ORA", ora_dn_title = "rescue-down · ORA"
  ),
  list(
    id = "secondary_v_disease", rank_old = "MITOvPHE_MITO", rank_form = "PHE_Mito - Mito", set_old = "CTLvPHE",
    rank_label = "Rank (PHE-on-transplant t, n = %d)", title = "Does PHE still perturb transplanted cells",
    subtitle = "rank = Phe+Mito − Mito (PHE effect on a transplant background); set = Π disease",
    up_title = "Disease-up set", dn_title = "Disease-down set",
    ora_up_title = "disease-up · ORA", ora_dn_title = "disease-down · ORA"
  )
)
mountain_stats <- purrr::map_dfr(MOUNTAINS, make_mountain)

pa <- build_signature_reversal(set_rev)
pb <- build_pca_trajectory(scores, centroids, var_exp, dist_stat)
pc_panel <- build_return_key(prot)
pd <- build_interaction_null(comb)

fig <- (add_tag(pa$plot, "A") | add_tag(pb$plot, "B")) /
  (add_tag(pc_panel$plot, "C") | add_tag(pd$plot, "D"))
ggplot2::ggsave(file.path(RPT, "PILOT_reversal.png"), fig,
  width = 250, height = 185, units = "mm", dpi = 300, limitsize = FALSE
)

build_workbook(
  file.path(DAT, "pilot_reversal.xlsx"),
  figure_title = "Pilot: is the PHE proteome shift reversed by transplant?",
  sheet_specs = list(
    list(name = "signature_reversal", df = pa$table, role = "Panel A", contents = "fry direction and FDR for the disease sets and cardiac-stress pathways, Disease vs Recovery"),
    list(name = "pca_distance", df = pb$table, role = "Panel B", contents = "group-centroid distances to Ctl and the label-permutation p"),
    list(name = "protein_return", df = pc_panel$table, role = "Panel C", contents = "disease-signature proteins: disease vs recovery log2FC, direction, MitoCarta lens, normalized flag"),
    list(name = "interaction", df = pd$table, role = "Panel D", contents = "interaction contrast per protein: log2FC, p, FDR, pi")
  )
)

message("pilot_reversal built")
