# MitoCarta-restricted enrichment helpers for F0X_mito variants.
# Builds on pathway_utils.R; uses only the MitoCarta subset of rat_gene_sets.rds.

# Mito-keyword regex applied to pathway names from Hallmark/Reactome/KEGG/GO:BP.
# MitoCarta pathways pass through unfiltered (prefix MITOCARTA_*). Single
# source of truth for all _mito figure variants (F03_mito / F04_mito / F05_mito
# / F06_mito).
MITO_PATHWAY_REGEX <- paste(
  "MITOCHOND",
  "OXPHOS|OXIDATIVE_PHOSPH",
  "AEROBIC_RESPIRATION|ELECTRON_TRANSPORT_CHAIN|RESPIRATORY_(CHAIN|ELECTRON)",
  "CITRATE_CYCLE|TRICARBOXYLIC|KREBS_CYCLE",
  "MITOPHAG",
  "MITO.+FISSION|FISSION.+MITO|MITO.+FUSION|FUSION.+MITO",
  "CRISTAE",
  "PROTEIN_IMPORT_INTO_MITO|MITOCH.+PROTEIN_IMPORT",
  "MICOS",
  "BETA_OXIDATION",
  "CARDIOLIPIN",
  "CARNITINE",
  "HEME_BIO|HEME_METAB",
  "FE_S_CLUSTER|IRON_SULFUR_CLUSTER",
  "PYRUVATE_(METABOLISM|DEHYDROGENASE)",
  "KETONE_(BODY|BIOSYNTH|METAB)",
  "ONE_CARBON_(POOL|METABOL)|GLYCINE_CLEAVAGE",
  "MITORIBOSOM",
  "NADH_DEHYDROGENASE|SUCCINATE_DEHYDROGENASE|CYTOCHROME_C_OXIDASE|ATP_SYNTHASE",
  "INTRINSIC_(APOPTOTIC|PATHWAY_FOR_APOPTOSIS)",
  "RELEASE_OF_CYTOCHROME_C|CYTOCHROME_COMPLEX_ASSEMBLY",
  sep = "|")

# Subset a multi-DB pathway list (output of build_pathway_collection) to the
# mito lens: keep MitoCarta entries unchanged + retain other-DB entries whose
# name matches MITO_PATHWAY_REGEX. Use for ORA panels in F0X_mito variants.
mitofilter_pathway_collection <- function(pw_collection,
                                          regex = MITO_PATHWAY_REGEX) {
  keep <- grepl("^MITOCARTA_", names(pw_collection)) |
            grepl(regex, names(pw_collection), perl = TRUE)
  pw_collection[keep]
}

# Sub-mitochondrial compartment palette (consistent across all _mito figures).
mitocarta_compartment_palette <- c(
  Matrix        = "#3D5A80",  # deep blue
  IMM           = "#98C1D9",  # mid blue (inner membrane)
  OMM           = "#EE6C4D",  # warm coral (outer membrane)
  IMS           = "#E0FBFC",  # pale cyan (intermembrane space)
  OXPHOS        = "#293241",  # near-black
  Metabolism    = "#81B29A",  # muted green
  `Mitochondrial dynamics` = "#F2CC8F",  # warm tan
  Other         = "#BFC0C0"
)

# Load only the MitoCarta block from the shared rat_gene_sets cache.
load_mitocarta_collection <- function(
  cache_path = here::here("04_Figures", "shared", "rat_gene_sets.rds"),
  drop_all   = TRUE
) {
  if (!file.exists(cache_path)) {
    stop("Cache not found: ", cache_path,
         "\nRun 04_Figures/shared/fetch_rat_gene_sets.R first.")
  }
  gs <- readRDS(cache_path)
  if (!"MitoCarta" %in% names(gs)) {
    stop("rat_gene_sets.rds has no 'MitoCarta' block. Re-run fetch_rat_gene_sets.R.")
  }
  mc <- gs[["MitoCarta"]]
  if (drop_all) mc[["MITOCARTA_ALL"]] <- NULL
  mc
}

# Split MitoCarta sets into named "tiers" for figure logic:
#   compartment: MATRIX / IMM / OMM / IMS
#   oxphos:      OXPHOS > Complex I..V (parsed from MITOCARTA_OXPHOS__...)
#   pathway:     all other MitoPathway sets
classify_mitocarta_tier <- function(set_names) {
  dplyr::case_when(
    set_names %in% c("MITOCARTA_MATRIX", "MITOCARTA_IMM",
                     "MITOCARTA_OMM",    "MITOCARTA_IMS") ~ "compartment",
    grepl("^MITOCARTA_OXPHOS",  set_names) ~ "oxphos",
    set_names == "MITOCARTA_ALL"           ~ "all",
    TRUE                                   ~ "pathway"
  )
}

# Run fGSEA restricted to MitoCarta sets. Returns dedup'd long-form data.frame
# matching pathway_utils.R::run_enrichment_pipeline() conventions.
run_fgsea_mitocarta <- function(stats_list,
                                jaccard_cutoff = 0.5,
                                nperm = 10000,
                                min_size = 5, max_size = 500,
                                padj_cutoff = 0.05) {
  if (!exists("run_enrichment_pipeline", mode = "function")) {
    source(here::here("04_Figures", "shared", "pathway_utils.R"))
  }
  mc <- load_mitocarta_collection()
  run_enrichment_pipeline(
    stats_list    = stats_list,
    pw_list       = mc,
    jaccard_cutoff = jaccard_cutoff,
    nperm         = nperm,
    min_size      = min_size,
    max_size      = max_size,
    padj_cutoff   = padj_cutoff
  )
}

# ORA restricted to MitoCarta sets (mirrors run_ora_deduplicated).
run_ora_mitocarta <- function(genes, universe,
                              jaccard_cutoff = 0.5,
                              min_size = 5, max_size = 500,
                              padj_cutoff = 0.05) {
  if (!exists("run_ora_deduplicated", mode = "function")) {
    source(here::here("04_Figures", "shared", "pathway_utils.R"))
  }
  mc <- load_mitocarta_collection()
  run_ora_deduplicated(
    genes        = genes,
    universe     = universe,
    pathways     = mc,
    jaccard_cutoff = jaccard_cutoff,
    min_size     = min_size,
    max_size     = max_size,
    padj_cutoff  = padj_cutoff
  )
}

# Human-readable label for a MITOCARTA_* set name, used in plot axis text.
prettify_mitocarta_label <- function(set_names) {
  s <- sub("^MITOCARTA_", "", set_names)
  s <- gsub("__", " > ", s)
  s <- gsub("_", " ", s)
  s <- tools::toTitleCase(tolower(s))
  s <- gsub("Oxphos", "OXPHOS", s)
  s <- gsub("Imm", "IMM", s, fixed = TRUE)
  s <- gsub("Omm", "OMM", s, fixed = TRUE)
  s <- gsub("Ims", "IMS", s, fixed = TRUE)
  s <- gsub("Tca", "TCA", s, fixed = TRUE)
  s
}
