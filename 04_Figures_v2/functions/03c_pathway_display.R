# Pathway display labels + functional category dictionaries. Depends on
# clean_pathway_name() from 01_style_palettes_theme.R.

# Single shared display-name dictionary (was duplicated as DISPLAY_LABELS in each
# _panel_A_quadrant.R and display_overrides in each 01_main_panels.R). Keys are
# clean_pathway_name() output (Title Case, "_"->" ", DB prefix stripped,
# MitoCarta reduced to its most-specific last segment).
PATHWAY_DISPLAY_OVERRIDES <- c(
  "Oxidative Phosphorylation"                              = "OXPHOS",
  "Oxphos Subunits"                                        = "OXPHOS Subunits",
  "Cv Subunits"                                            = "Complex V Subunits",
  "Ci Subunits"                                            = "Complex I Subunits",
  "Civ Subunits"                                           = "Complex IV Subunits",
  "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic Resp. + ETC",
  "Respiratory Electron Transport"                         = "Respiratory ETC",
  "Complex I Biogenesis"                                   = "Complex I Biogenesis",
  "Citric Acid Cycle Tca Cycle"                            = "TCA Cycle",
  "Tca Cycle"                                              = "TCA Cycle",
  "Mitochondrial Translation"                              = "Mito Translation",
  "Mitochondrial Ribosome"                                 = "Mitoribosome",
  "Mitochondrial Protein Degradation"                      = "Mito Protein Deg.",
  "Mitochondrial Calcium Ion Transport"                    = "Mito Ca²⁺ Transport",
  "Pink1 Prkn Mediated Mitophagy"                          = "PINK1/PRKN Mitophagy",
  "Fatty Acid Beta Oxidation"                              = "FA β-Oxidation",
  "Fatty Acid Oxidation"                                   = "FA Oxidation",
  "Fatty Acid Metabolism"                                  = "FA Metabolism",
  "Branched Chain Amino Acid Metabolism"                   = "BCAA Metabolism",
  "Cardiolipin Biosynthesis"                               = "Cardiolipin Biosynth.",
  "Mitochondrion Organization"                             = "Mito Organization",
  "Mitochondrial Organization"                             = "Mito Org.",
  "Mitochondrial Transport"                                = "Mito Transport",
  "Mitochondrial Protein Import"                           = "Mito Protein Import",
  "Small Molecule Transport"                               = "Sm. Molecule Transport",
  "Slc25a Family"                                          = "SLC25A Carriers",
  "Proteases"                                              = "Mito Proteases",
  "Ros And Glutathione Metabolism"                         = "ROS/Glutathione Metab.",
  "Xenobiotic Metabolism"                                  = "Xenobiotic Metab.",
  "Lysine Metabolism"                                      = "Lysine Metab.",
  "Amino Acid Metabolic Process"                           = "Amino Acid Metab.",
  "Amino Acid Metabolism"                                  = "Amino Acid Metab.",
  "Extracellular Matrix Organization"                      = "ECM Organization",
  "Epithelial Mesenchymal Transition"                      = "EMT",
  "Cytoplasmic Translation"                                = "Cytoplasmic Transl.",
  "Glycolysis Gluconeogenesis"                             = "Glycolysis/Gluconeo.",
  "Generation Of Precursor Metabolites And Energy"         = "Precursor Metab. & Energy",
  "Organophosphate Metabolic Process"                      = "Organophosphate Metab.",
  "Organic Acid Metabolic Process"                         = "Organic Acid Metab.",
  "Carboxylic Acid Metabolic Process"                      = "Carboxylic Acid Metab.",
  "Monocarboxylic Acid Metabolic Process"                  = "Monocarbox. Acid Metab.",
  "Small Molecule Catabolic Process"                       = "Sm. Molecule Catab.",
  "Small Molecule Metabolic Process"                       = "Sm. Molecule Metab.",
  "Ketone Metabolic Process"                               = "Ketone Metabolism",
  "Protein Localization To Plasma Membrane"                = "Plasma Membr. Protein Loc.",
  "Striated Muscle Contraction"                            = "Striated Muscle Contr."
)

# Human-readable pathway label for any DB. MitoCarta names are a
# "MITOCARTA_A__B__TERM" hierarchy -> keep the most specific last segment, drop
# the prefix; everything else goes through clean_pathway_name(). Then apply the
# shared override dictionary (extra = figure-specific extension, merged on top).
clean_display_label <- function(pathway, extra = NULL) {
  is_mito <- grepl("^MITOCARTA_", pathway)
  lab <- ifelse(
    is_mito,
    clean_pathway_name(gsub("_", " ", sub(".*__", "", sub("^MITOCARTA_", "", pathway)))),
    clean_pathway_name(pathway))
  dict <- PATHWAY_DISPLAY_OVERRIDES
  if (!is.null(extra)) dict[names(extra)] <- extra
  hit <- match(lab, names(dict))
  ifelse(!is.na(hit), unname(dict[hit]), lab)
}

# MSigDB pathway ID -> 15 consolidated categories (keyword rules)
CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

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

classify_pathway_func <- function(ids) {
  rules <- list(
    "Muscle & Contractile"              = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND",
    "Cytoskeleton & Motility"           = "CYTOSKELET|ACTIN_BIND|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL",
    "ECM & Adhesion"                    = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING",
    "Lipid Metabolism"                  = "FATTY_ACID|LIPID|ADIPOGEN|STEROID|SPHINGOLIPID|PHOSPHOLIPID|KETONE",
    "Carbohydrate & Energy Metabolism"  = "GLYCOLY|GLUCONEO|CARBOHYDRATE|PENTOSE|PRECURSOR_METABOL",
    "Amino Acid & Cofactor Metabolism"  = "AMINO_ACID|VITAMIN|COFACTOR|NITROGEN|DETOXIF|DIGEST|XENOBIOT",
    "Mitochondria & Energy"             = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT",
    "Protein Homeostasis"               = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD",
    "Transport"                         = "TRANSPORT(?!.*ELECTRON)|VESICLE|ENDOCYT|SECRETI",
    "Translation & Ribosome"            = "TRANSLAT|RIBOSOM|TRNA|MYC_TARGET",
    "Transcription & Chromatin"         = "TRANSCRIPT|SPLICEOSOM|E2F_TARGET|CHROMATIN|MRNA_PROC",
    "Immune & Inflammation"             = "IMMUN|INFLAMMA|INTERFERON|IL2|IL6|TNFA|NF.KB|COMPLEMENT",
    "DNA & Cell Cycle"                  = "DNA_REPAIR|CELL_CYCLE|MITOTIC|P53_PATHWAY",
    "Circulatory System"                = "ANGIOGEN|BLOOD_VESSEL|HYPOXIA",
    "Development"                       = "UV_RESPONSE|GROWTH_FACTOR|WNT|HEDGEHOG|NOTCH|TGF_BETA|KRAS"
  )
  vapply(toupper(ids), function(id) {
    matches <- character(0)
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) matches <- c(matches, cat)
    }
    if (length(matches) > 1) {
      warning("classify_pathway_func: '", id, "' matches multiple categories [",
              paste(matches, collapse = ", "), "]; using first match: ", matches[1])
    }
    if (length(matches) >= 1) return(matches[1])
    if (grepl("METABOL", id, perl = TRUE)) return("Amino Acid & Cofactor Metabolism")
    "Other"
  }, character(1), USE.NAMES = FALSE)
}
