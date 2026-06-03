# GO Slim pathway assignment — shared by F2/panel_F.R and F3/panel_F.R
# Maps genes to 15 consolidated pathways via GO Slim Generic BP + GOBPANCESTOR traversal.

requireNamespace("GO.db",       quietly = TRUE)
requireNamespace("org.Rn.eg.db", quietly = TRUE)
requireNamespace("AnnotationDbi", quietly = TRUE)
library(dplyr)
library(tidyr)

# Mito (H9c2 rat cardiomyoblast) uses org.Rn.eg.db (Rattus norvegicus).
# Gene symbols in 03_combined_results.csv are TitleCase rat symbols (e.g. "Abcf2").
.SLIM_ORG_DB <- org.Rn.eg.db::org.Rn.eg.db

bp_slim <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135"
)

# signaling + nervous system excluded — too broad / irrelevant for muscle

SLIM_CONSOLIDATED <- c(
  "GO:0003012" = "Muscle & Contractile",
  "GO:0003013" = "Circulatory System",
  "GO:0030198" = "ECM & Adhesion",      "GO:0007155" = "ECM & Adhesion",
  "GO:0034330" = "ECM & Adhesion",      "GO:0042060" = "ECM & Adhesion",
  "GO:0007010" = "Cytoskeleton & Motility", "GO:0048870" = "Cytoskeleton & Motility",
  "GO:0007018" = "Cytoskeleton & Motility", "GO:0007163" = "Cytoskeleton & Motility",
  "GO:0044782" = "Cytoskeleton & Motility",
  "GO:0002376" = "Immune & Inflammation", "GO:0006954" = "Immune & Inflammation",
  "GO:0098542" = "Immune & Inflammation",
  "GO:0006629" = "Lipid Metabolism",     "GO:0042180" = "Lipid Metabolism",
  "GO:0005975" = "Carbohydrate & Energy Metabolism",
  "GO:0006091" = "Carbohydrate & Energy Metabolism",
  "GO:1901135" = "Carbohydrate & Energy Metabolism",
  "GO:0006520" = "Amino Acid & Cofactor Metabolism",
  "GO:0055086" = "Amino Acid & Cofactor Metabolism",
  "GO:0006766" = "Amino Acid & Cofactor Metabolism",
  "GO:0071941" = "Amino Acid & Cofactor Metabolism",
  "GO:0098754" = "Amino Acid & Cofactor Metabolism",
  "GO:0007586" = "Amino Acid & Cofactor Metabolism",
  "GO:0007005" = "Mitochondria & Energy", "GO:0007031" = "Mitochondria & Energy",
  "GO:0006457" = "Protein Homeostasis",  "GO:0030163" = "Protein Homeostasis",
  "GO:0006914" = "Protein Homeostasis",  "GO:0051604" = "Protein Homeostasis",
  "GO:0065003" = "Protein Homeostasis",  "GO:0009100" = "Protein Homeostasis",
  "GO:0055085" = "Transport",            "GO:0016192" = "Transport",
  "GO:0006886" = "Transport",            "GO:0006913" = "Transport",
  "GO:0072659" = "Transport",            "GO:0061024" = "Transport",
  "GO:0002181" = "Translation & Ribosome", "GO:0042254" = "Translation & Ribosome",
  "GO:0006399" = "Translation & Ribosome",
  "GO:0006351" = "Transcription & Chromatin", "GO:0006355" = "Transcription & Chromatin",
  "GO:0016071" = "Transcription & Chromatin", "GO:0006325" = "Transcription & Chromatin",
  "GO:0006281" = "DNA & Cell Cycle",     "GO:0006260" = "DNA & Cell Cycle",
  "GO:0006310" = "DNA & Cell Cycle",     "GO:0032200" = "DNA & Cell Cycle",
  "GO:0000278" = "DNA & Cell Cycle",     "GO:0140014" = "DNA & Cell Cycle",
  "GO:0007059" = "DNA & Cell Cycle",     "GO:0000910" = "DNA & Cell Cycle",
  "GO:0007126" = "DNA & Cell Cycle",
  "GO:0048856" = "Development",          "GO:0030154" = "Development",
  "GO:0012501" = "Development",          "GO:0003014" = "Development",
  "GO:0003016" = "Development"
)


if (!exists("CONSOLIDATED_PATHWAY_ORDER")) {
  CONSOLIDATED_PATHWAY_ORDER <- c(
    "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
    "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
    "Amino Acid & Cofactor Metabolism",
    "Mitochondria & Energy", "Protein Homeostasis",
    "Transport", "Translation & Ribosome", "Transcription & Chromatin",
    "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
    "Development", "Other"
  )
}

if (!exists("CONSOLIDATED_COLORS")) {
  CONSOLIDATED_COLORS <- c(
    "Muscle & Contractile"              = "#E57373",
    "Cytoskeleton & Motility"           = "#FFB74D",
    "ECM & Adhesion"                    = "#FFF176",
    "Lipid Metabolism"                  = "#AED581",
    "Carbohydrate & Energy Metabolism"  = "#81C784",
    "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
    "Mitochondria & Energy"             = "#4DB6AC",
    "Protein Homeostasis"               = "#4FC3F7",
    "Transport"                         = "#7986CB",
    "Translation & Ribosome"            = "#BA68C8",
    "Transcription & Chromatin"         = "#AB47BC",
    "Immune & Inflammation"             = "#A1887F",
    "DNA & Cell Cycle"                  = "#90A4AE",
    "Circulatory System"                = "#CE93D8",
    "Development"                       = "#B0BEC5",
    "Other"                             = "#D0D0D0"
  )
}

assign_go_slim_consolidated <- function(fg_genes, all_genes, min_cat_size = 2) {

  suppressMessages({
    all_entrez <- AnnotationDbi::mapIds(.SLIM_ORG_DB, keys = all_genes,
                        keytype = "SYMBOL", column = "ENTREZID",
                        multiVals = "first")
    all_go <- AnnotationDbi::select(.SLIM_ORG_DB,
                 keys = as.character(na.omit(all_entrez)),
                 keytype = "ENTREZID",
                 columns = c("SYMBOL", "GO", "ONTOLOGY"))
  })
  all_bp <- all_go |>
    filter(ONTOLOGY == "BP", !is.na(GO)) |>
    distinct(SYMBOL, GO)

  ancestors  <- as.list(GO.db::GOBPANCESTOR)
  all_go_ids <- unique(all_bp$GO)

  go_to_slim <- setNames(
    lapply(all_go_ids, function(go_id) {
      hits <- character(0)
      if (go_id %in% bp_slim) hits <- go_id
      anc <- ancestors[[go_id]]
      if (!is.null(anc)) hits <- c(hits, intersect(anc, bp_slim))
      unique(hits)
    }),
    all_go_ids
  )

  all_gene_slim <- all_bp |>
    mutate(slim_list = go_to_slim[GO]) |>
    unnest(slim_list) |>
    dplyr::select(SYMBOL, slim = slim_list) |>   # qualify: AnnotationDbi may mask select
    distinct()

  all_gene_consolidated <- all_gene_slim |>
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) |>
    filter(!is.na(consolidated)) |>
    distinct(SYMBOL, consolidated)

  fg_gene_slim <- all_gene_slim |> filter(SYMBOL %in% fg_genes)
  fg_gene_consolidated <- fg_gene_slim |>
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) |>
    filter(!is.na(consolidated))

  fg_term_counts <- fg_gene_slim |> count(slim, name = "n_fg")

  best_consolidated <- fg_gene_consolidated |>
    left_join(fg_term_counts, by = "slim") |>
    mutate(priority = ifelse(consolidated == "Development", 2, 1)) |>
    arrange(priority, n_fg) |>
    group_by(SYMBOL) |>
    slice_head(n = 1) |>
    ungroup()

  unmapped <- setdiff(fg_genes, best_consolidated$SYMBOL)
  if (length(unmapped) > 0) {
    best_consolidated <- bind_rows(best_consolidated,
      tibble(SYMBOL = unmapped, slim = "OTHER", consolidated = "Other",
             n_fg = NA_integer_, priority = 3L))
  }

  small_cats <- best_consolidated |> count(consolidated) |>
    filter(n < min_cat_size, consolidated != "Other") |> pull(consolidated)
  if (length(small_cats) > 0) {
    best_consolidated <- best_consolidated |>
      mutate(consolidated = ifelse(consolidated %in% small_cats, "Other", consolidated))
  }

  best_consolidated |>
    transmute(gene = SYMBOL, slim, consolidated) |>
    mutate(consolidated = factor(consolidated, levels = CONSOLIDATED_PATHWAY_ORDER))
}

