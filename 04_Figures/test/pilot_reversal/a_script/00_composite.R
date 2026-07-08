#!/usr/bin/env Rscript
# Pilot: is the PHE proteome shift reversed by transplant? Four panels build the
# claim from set-level down to per-protein: signature reversal (fry on the disease
# sets and a priori cardiac-stress pathways, tested Disease -> Recovery), a PCA
# trajectory showing each group's distance to control, the per-protein return with
# a mito / direction key, and the formal interaction cap.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "01_style_palettes_theme.R"))
source(file.path(fns, "02_data_paths_and_loaders.R"))
source(file.path(fns, "06_supplementary_workbook.R"))
source(file.path(fns, "08_composite_layout.R"))
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

# --- fry mountain: disease sets ranked by the Rescue t (transplant under stress) ---
cm_r <- makeContrasts("PHE_Mito - PHE", levels = design)
fry_rescue <- function(idx, expected) {
  fr <- fry(expr, idx, design, contrast = cm_r[, 1], block = rep_block, correlation = set_rho)
  tibble(dir = fr$Direction[1], fdr = fdr_col(fr, rownames(fr))[1], n = length(idx), consistent = fr$Direction[1] == expected)
}
fry_up <- fry_rescue(idx_sig$`disease-up`, "Down")
fry_dn <- fry_rescue(idx_sig$`disease-down`, "Up")

up_ids <- comb$uniprot_id[comb$sig_pi_CTLvPHE == 1]
dn_ids <- comb$uniprot_id[comb$sig_pi_CTLvPHE == -1]
t_rank <- comb |>
  filter(!is.na(t_PHEvPHE_MITO)) |>
  transmute(uniprot_id, gene, t_rescue = t_PHEvPHE_MITO, in_up = uniprot_id %in% up_ids, in_down = uniprot_id %in% dn_ids) |>
  arrange(desc(t_rescue)) |>
  mutate(rank = dplyr::row_number(), es_up = running_es(t_rescue, in_up), es_down = running_es(t_rescue, in_down))
n_all <- nrow(t_rank)
circ_r <- cor(comb$t_CTLvPHE, comb$t_PHEvPHE_MITO, use = "complete.obs")

pw_flat <- unlist(unname(readRDS(here::here("04_Figures", "shared", "c_data", "rat_gene_sets.rds"))), recursive = FALSE)
universe <- unique(comb$gene[!is.na(comb$gene)])
run_ora <- function(g) {
  if (length(g) < 5) {
    return(tibble())
  }
  fora(pathways = pw_flat, genes = unique(g), universe = universe, minSize = 10, maxSize = 500) |>
    as_tibble() |>
    filter(padj < 0.1) |>
    arrange(padj) |>
    mutate(.key = .pretty_pw(pathway)) |>
    distinct(.key, .keep_all = TRUE) |>
    select(-.key)
}
# ORA characterises the broader disease-responsive reversal (Disease p<0.05), since
# the strict Pi signature is too small to over-represent; the fry test above keeps Pi.
broad <- comb |> filter(P.Value_CTLvPHE < 0.05, !is.na(t_PHEvPHE_MITO))
ora_up <- run_ora(broad$gene[broad$logFC_CTLvPHE > 0 & broad$t_PHEvPHE_MITO < 0])
ora_dn <- run_ora(broad$gene[broad$logFC_CTLvPHE < 0 & broad$t_PHEvPHE_MITO > 0])

pa <- build_signature_reversal(set_rev)
pb <- build_pca_trajectory(scores, centroids, var_exp, dist_stat)
pc_panel <- build_return_key(prot)
pd <- build_interaction_null(comb)
mountain <- build_fry_mountain(t_rank, fry_up, fry_dn, ora_up, ora_dn, circ_r, n_all)
ggplot2::ggsave(file.path(RPT, "PILOT_reversal_fry.png"), mountain,
  width = 220, height = 150, units = "mm", dpi = 300, limitsize = FALSE
)

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
