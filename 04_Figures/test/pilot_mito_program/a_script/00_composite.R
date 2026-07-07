#!/usr/bin/env Rscript
# Pilot: did transplant just add mitochondrial proteins, or reprogram the cell?
# Programs are MitoCarta3.0 MitoPathways (Rath 2021). Tested with the same engine
# as the WGCNA modules: fry (self-contained, valid at n=6; Wu 2010) as the gate,
# camera (competitive, correlation-adjusted; Wu & Smyth 2012) for selectivity vs
# the rest of the proteome, on the imputed matrix with the paired Replicate block.
# Median log2FC is descriptive only and never the inferential claim.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "01_style_palettes_theme.R"))
source(file.path(fns, "02_data_paths_and_loaders.R"))
source(file.path(fns, "06_supplementary_workbook.R"))
source(file.path(fns, "08_composite_layout.R"))
pacman::p_load(limma, dplyr, tidyr, tibble, purrr, patchwork)

BASE <- here::here("04_Figures", "test", "pilot_mito_program")
for (f in list.files(file.path(BASE, "a_script", "panels"), full.names = TRUE)) source(f)
RPT <- file.path(BASE, "b_reports")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

CAT_LABELS <- c(
  OXPHOS = "OXPHOS",
  MITOCHONDRIAL_CENTRAL_DOGMA = "mtDNA & translation",
  METABOLISM = "Metabolism",
  MITOCHONDRIAL_DYNAMICS_AND_SURVEILLANCE = "Dynamics & mitophagy",
  PROTEIN_IMPORT_SORTING_AND_HOMEOSTASIS = "Import & proteostasis",
  SIGNALING = "Signaling",
  SMALL_MOLECULE_TRANSPORT = "Small-molecule transport"
)
STRUCTURAL_CATS <- c("OXPHOS", "MITOCHONDRIAL_CENTRAL_DOGMA")

comb <- load_combined_wide()

mito <- readRDS(here::here("04_Figures", "shared", "mitocarta3_rat.rds"))
mito_genes <- toupper(mito$MITOCARTA_ALL)
cat_long <- mito[grepl("__", names(mito))] |>
  imap(~ tibble(gene = toupper(.x), category = sub("__.*$", "", sub("^MITOCARTA_", "", .y)))) |>
  bind_rows() |>
  distinct(gene, category)
cat_genes <- split(cat_long$gene, cat_long$category)

# self-contained fry + competitive camera per program, on the gap-free imputed
# matrix (fry/camera need a complete matrix), same design + block as the DE.
dal <- readRDS(P05$imp_rds)
expr <- as.matrix(dal$data)
meta <- as.data.frame(dal$metadata)
grp <- factor(meta$Group[match(colnames(expr), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)
rep_block <- meta$Replicate[match(colnames(expr), meta$Col_ID)]
design <- stats::model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
CTR <- c(Transplant = "Mito - Ctl", Rescue = "PHE_Mito - PHE", Interaction = "(PHE_Mito - PHE) - (Mito - Ctl)")
cm <- makeContrasts(contrasts = unname(CTR), levels = design)
colnames(cm) <- names(CTR)

genes <- toupper(dal$annotation$gene[match(rownames(expr), dal$annotation$uniprot_id)])
cat_genes <- cat_genes[vapply(cat_genes, \(g) sum(genes %in% g) >= 5L, logical(1))]
idx_cat <- lapply(cat_genes, \(g) which(genes %in% g))
idx_all <- list(MitoCarta = which(genes %in% mito_genes))

set_rho <- duplicateCorrelation(expr, design, block = rep_block)$consensus
# camera/fry omit the FDR column for a single set; fall back to the p-value there.
fdr_col <- function(df, rn) if ("FDR" %in% colnames(df)) df[rn, "FDR"] else df[rn, "PValue"]
run_sets <- function(idx) {
  bind_rows(lapply(names(CTR), function(cn) {
    fr <- fry(expr, idx, design, contrast = cm[, cn], block = rep_block, correlation = set_rho)
    ca <- camera(expr, idx, design, contrast = cm[, cn])
    tibble(
      set = rownames(fr), contrast = cn,
      fry_dir = fr$Direction, fry_fdr = fdr_col(fr, rownames(fr)),
      camera_dir = ca[rownames(fr), "Direction"], camera_fdr = fdr_col(ca, rownames(fr))
    )
  }))
}
cat_tests <- run_sets(idx_cat) |>
  mutate(
    arm = ifelse(set %in% STRUCTURAL_CATS, "structural", "regulatory"),
    label = dplyr::coalesce(CAT_LABELS[set], set)
  )
all_tests <- run_sets(idx_all)

pa <- build_mito_vs_nonmito(comb, mito_genes, all_tests)
pb <- build_program_breakdown(cat_tests)
pc <- build_context_dependence(cat_tests)

fig <- add_tag(pa$plot, "A") | add_tag(pb$plot, "B") | add_tag(pc$plot, "C")
ggplot2::ggsave(file.path(RPT, "PILOT_mito_program.png"), fig,
  width = 260, height = 105, units = "mm", dpi = 300, limitsize = FALSE
)

build_workbook(
  file.path(DAT, "pilot_mito_program.xlsx"),
  figure_title = "Pilot: added mitochondrial cargo or proteome reprogramming?",
  sheet_specs = list(
    list(
      name = "mito_vs_nonmito", df = pa$table, role = "Panel A",
      contents = "MitoCarta vs other median log2FC per contrast, with the whole-set camera competitive test"
    ),
    list(
      name = "program_fry", df = pb$table, role = "Panel B",
      contents = "per MitoPathway: fry (self-contained) direction and FDR, Transplant and Rescue"
    ),
    list(
      name = "program_interaction", df = pc$table, role = "Panel C",
      contents = "per MitoPathway: fry FDR for the MT x PHE interaction contrast"
    )
  )
)

message("pilot_mito_program built")
