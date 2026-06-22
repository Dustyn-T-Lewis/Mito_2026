# MitoCarta helpers: the mito-keyword regex (single source of truth for the mito
# lens) and a loader for the MitoCarta block of rat_gene_sets.rds.

# Mito-keyword regex applied to Hallmark/Reactome/KEGG/GO-Slim pathway names so
# mito-relevant non-MitoCarta sets are flagged alongside the MITOCARTA_* sets.
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

# The MitoCarta block of the shared cache. drop_all = TRUE drops the MITOCARTA_ALL
# aggregate, leaving the sub-localization/pathway sets.
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
