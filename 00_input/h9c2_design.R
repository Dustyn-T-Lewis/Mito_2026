# H9c2 mito-transplantation proteomics: design constants, contrasts, palettes,
# thresholds, and metadata helpers. Sourced by each stage's _setup.R.
#
# Design: 4 groups (Ctl, Mito, PHE, PHE_Mito), n = 6 per group, paired by
# passage / plating day / plate (block: Replicate). PHE = phenylephrine
# (right-HF stress model); Mito = mitochondrial transplantation.

# Factor order — every script must use these levels.
H9C2_GROUP_LEVELS <- c("Ctl", "Mito", "PHE", "PHE_Mito")

# Raw label -> standardized Group. Originals kept in metadata$Original_Label.
H9C2_RAW_GROUP_MAP <- c(
  "CTL"      = "Ctl",
  "CTL-Mito" = "Mito",
  "PE"       = "PHE",
  "PE-Mito"  = "PHE_Mito"
)

H9C2_DESIGN_FORMULA <- "~ 0 + group + (1|Replicate)"

# Contrasts

H9C2_CONTRASTS <- c(
  "CTLvPHE       = PHE - Ctl",
  "CTLvMITO      = Mito - Ctl",
  "PHEvPHE_MITO  = PHE_Mito - PHE",
  "Interaction   = (PHE_Mito - PHE) - (Mito - Ctl)",
  "MITOvPHE_MITO = PHE_Mito - Mito"
)

H9C2_CORE_CONTRASTS <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction")

H9C2_CONTRAST_ROLES <- c(
  CTLvPHE       = "Disease",
  CTLvMITO      = "Intervention",
  PHEvPHE_MITO  = "Rescue",
  Interaction   = "Interaction",
  MITOvPHE_MITO = "Secondary"
)

# Subtitle labels — printed verbatim in panels. Unicode minus (U+2212).
H9C2_CONTRAST_LABELS <- c(
  CTLvPHE       = "Disease  =  PHE − Ctl",
  CTLvMITO      = "Intervention  =  Mito − Ctl",
  PHEvPHE_MITO  = "Rescue  =  PHE_Mito − PHE",
  Interaction   = "Interaction  =  Rescue − Intervention",
  MITOvPHE_MITO = "Secondary  =  PHE_Mito − Mito"
)

H9C2_CONTRAST_DESC <- c(
  CTLvPHE       = "PHE-induced proteomic remodeling (disease/stress effect)",
  CTLvMITO      = "Mitochondrial transplantation-alone effect",
  PHEvPHE_MITO  = "Mitochondrial transplantation effect under PHE stress (rescue)",
  Interaction   = "PHE-dependent mitochondrial transplantation response",
  MITOvPHE_MITO = "PHE effect in mito-transplanted cells"
)

# Palettes

H9C2_PAL_GROUP <- c(
  Ctl      = "#4393C3",  # blue   — control
  Mito     = "#009E73",  # green  — intervention (Okabe-Ito, colourblind-safe)
  PHE      = "#D6604D",  # red    — disease
  PHE_Mito = "#984EA3"   # purple — disease + intervention
)

H9C2_PAL_DIR <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

H9C2_PAL_CONTRAST <- c(
  CTLvPHE       = "#D6604D",  # Disease
  CTLvMITO      = "#4393C3",  # Intervention
  PHEvPHE_MITO  = "#4DAF4A",  # Rescue
  Interaction   = "#7B5EA7",  # Interaction
  MITOvPHE_MITO = "#FF8C00"   # Secondary
)

# Thresholds
# Pi-score = P.Value ^ |logFC|  (Xiao et al. 2014, PMID 22321699).
# Collaborator's DAP list used p/FDR/pi = 0.05 and |lfc| >= 0.6.
H9C2_PVAL_THRESH <- 0.10
H9C2_PI_THRESH   <- 0.05
H9C2_FDR_EXPLOR  <- 0.10

# Contaminants
# Keratins (handling): rat symbols are Title-case. Match gene OR description
# so LOC-named keratin entries without a Krt* symbol are still caught.
H9C2_KERATIN_GENE_PATTERN <- "^Krt[0-9]|^Krtap"
H9C2_KERATIN_DESC_PATTERN <- "^Keratin"

# FBS carry-over from culture medium: Geyer (2016) plasma panel,
# case-insensitive on gene symbol.
H9C2_SERUM_CONTAMINANTS <- c(
  "Hba1", "Hba2", "Hbb", "Hba-a1", "Hbb-b1",            # hemoglobins
  "Alb", "Trf", "Tf", "Hp", "Hpx", "Gc",                # carriers
  "Apoa1", "Apoa2", "Apob", "Apoc1", "Apoc2", "Apoc3", "Apoe",
  "Fga", "Fgb", "Fgg", "F2", "Plg",                     # coagulation
  "C3", "C4a", "C4b", "C5", "C6", "C7", "C8a", "C8b", "C8g", "C9",
  "Cfb", "Cfh", "Cfi", "C1qa", "C1qb", "C1qc", "C1r", "C1s", "C2",
  "Serpina1", "Serpina3", "A2m", "Orm1", "Orm2",        # acute-phase
  "Ahsg", "Fetub", "Itih4", "Agt", "Ambp", "Kng1", "Hrg", "Vtn")

# Filter / outlier constants
H9C2_MIN_REPS  <- 4L      # detected in >= 4 of 6 in >= 1 group
H9C2_OUTLIER_K <- 3L      # outlier: flagged by >= K of 4 methods
H9C2_MAHAL_P   <- 0.01    # PCA-Mahalanobis tail probability
H9C2_MAD_K     <- 3       # MAD multiplier for intensity / correlation flags

# Imputation constant
H9C2_MISS_UNRELIABLE <- 50  # % missing above which an imputed protein is flagged

# Helpers

h9c2_strip_raw <- function(x) sub("\\.raw$", "", x)

load_h9c2_metadata <- function(path) {
  meta <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("Col_ID", "Sample_ID", "Group", "Original_Label",
                "Replicate", "Reinjected")
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L)
    stop("H9c2 metadata missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(meta$Col_ID))
    stop("H9c2 metadata has duplicate Col_ID values")
  if (!setequal(unique(meta$Group), H9C2_GROUP_LEVELS))
    stop("H9c2 metadata Group must be exactly {",
         paste(H9C2_GROUP_LEVELS, collapse = ", "), "}; found {",
         paste(sort(unique(meta$Group)), collapse = ", "), "}")
  meta$Group <- factor(meta$Group, levels = H9C2_GROUP_LEVELS)
  rownames(meta) <- meta$Col_ID
  meta
}

# Parse "name = expression" strings into a named character vector.
h9c2_parse_contrasts <- function(contrasts = H9C2_CONTRASTS) {
  parts <- strsplit(contrasts, "=", fixed = TRUE)
  nm    <- trimws(vapply(parts, `[`, character(1), 1L))
  rhs   <- trimws(vapply(parts, `[`, character(1), 2L))
  stats::setNames(rhs, nm)
}

assert_h9c2_group_sizes <- function(meta, min_n = 3L) {
  tab <- table(meta$Group)
  if (any(tab < min_n))
    stop("H9c2 design needs >= ", min_n, " samples per group; found ",
         paste(sprintf("%s=%d", names(tab), tab), collapse = ", "))
  invisible(tab)
}
