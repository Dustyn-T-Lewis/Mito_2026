#!/usr/bin/env Rscript
# Pilot: the transplant main effect. Pool both mito-receiving groups against both
# non-mito groups, (Mito + PHE_Mito)/2 - (Ctl + PHE)/2, the sixth factorial
# contrast the reported DE omits. It is interpretable because the MT x PHE
# interaction is ~null (pilot_mito_program panel C), which also doubles the
# samples behind the estimate. Fit on the imputed missForest matrix with the
# paired Replicate block, same design as the DE; refit on non-imputed if this
# graduates from a pilot to a reported figure.

fns <- here::here("04_Figures", "functions")
source(file.path(fns, "F01-F03_style_palettes_theme.R"))
source(file.path(fns, "F01-F03_data_paths_and_loaders.R"))
source(file.path(fns, "F01-F03_supplementary_workbook.R"))
source(file.path(fns, "F01-F03_composite_layout.R"))
pacman::p_load(limma, dplyr, tibble, patchwork)

BASE <- here::here("04_Figures", "test", "pilot_mito_main")
for (f in list.files(file.path(BASE, "a_script", "panels"), full.names = TRUE)) source(f)
RPT <- file.path(BASE, "b_reports")
DAT <- file.path(BASE, "c_data")
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(P05$imp_rds)
expr <- as.matrix(dal$data)
meta <- as.data.frame(dal$metadata)
grp <- factor(meta$Group[match(colnames(expr), meta$Col_ID)], levels = H9C2_GROUP_LEVELS)
rep_block <- meta$Replicate[match(colnames(expr), meta$Col_ID)]
design <- stats::model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
cm <- makeContrasts(MitoMain = (Mito + PHE_Mito) / 2 - (Ctl + PHE) / 2, levels = design)

dc <- duplicateCorrelation(expr, design, block = rep_block)
fit <- lmFit(expr, design, block = rep_block, correlation = dc$consensus)
fit <- eBayes(contrasts.fit(fit, cm))

mito <- readRDS(here::here("04_Figures", "shared", "mitocarta3_rat.rds"))
mito_genes <- toupper(mito$MITOCARTA_ALL)

main <- topTable(fit, coef = "MitoMain", number = Inf, sort.by = "none") |>
  as_tibble(rownames = "uniprot_id") |>
  mutate(
    gene = dal$annotation$gene[match(uniprot_id, dal$annotation$uniprot_id)],
    lens = ifelse(toupper(gene) %in% mito_genes, "MitoCarta", "other"),
    pi_score = P.Value^abs(logFC),
    hit = pi_score < H9C2_PI_THRESH
  )

pa <- build_mito_main_volcano(main)
pb <- build_driving_proteins(main)

fig <- add_tag(pa$plot, "A") | add_tag(pb$plot, "B")
ggplot2::ggsave(file.path(RPT, "PILOT_mito_main.png"), fig,
  width = 220, height = 100, units = "mm", dpi = 300, limitsize = FALSE
)

build_workbook(
  file.path(DAT, "pilot_mito_main.xlsx"),
  figure_title = "Pilot: the transplant main effect (mito vs non-mito, pooled over stress)",
  sheet_specs = list(
    list(
      name = "mito_main", df = pa$table, role = "Panel A",
      contents = "per protein: pooled main-effect logFC, p, FDR, pi-score, hit, MitoCarta flag"
    ),
    list(
      name = "drivers", df = pb$table, role = "Panel B",
      contents = "pi-score movers ranked by |moderated t|, signed logFC and MitoCarta flag"
    )
  )
)

message("pilot_mito_main built")
