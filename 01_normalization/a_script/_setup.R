# Stage 01 setup — sourced by 01_normalize.R and 02_generate_reports.R.

suppressPackageStartupMessages(library(here))

source(here::here("00_input", "h9c2_design.R"))

# Paths
INPUT <- here::here("00_input")
RPT   <- here::here("01_normalization", "b_reports")
DAT   <- here::here("01_normalization", "c_data")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Inputs
RAW_XLSX <- file.path(INPUT, "H9c2_raw.xlsx")
META_CSV <- file.path(INPUT, "H9c2_meta.csv")

# DIA-NN annotation columns (everything else is a sample column).
ANNOT_COLS_RAW <- c("Protein.Group", "Protein.Names", "Genes",
                    "First.Protein.Description", "N.Sequences",
                    "N.Proteotypic.Sequences")
