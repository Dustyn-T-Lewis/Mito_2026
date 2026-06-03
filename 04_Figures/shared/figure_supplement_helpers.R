# Shared helpers for 90_stitch_figure.R scripts: safe CSV reads, workbook assembly, cleanup.

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
  library(readxl)
})

add_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data)
  hs <- createStyle(textDecoration = "bold")
  addStyle(wb, name, hs, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
  cat(sprintf("    + %s: %d x %d\n", name, nrow(data), ncol(data)))
}

safe_read <- function(path) {
  if (file.exists(path)) {
    as.data.frame(read_csv(path))
  } else {
    cat(sprintf("    SKIP (not found): %s\n", path))
    NULL
  }
}

build_workbook <- function(out_file, title = NULL, description = NULL,
                            overview_df = NULL, sheet_specs) {
  # No Overview sheet: matches the no-Overview style of S01–S03 stage workbooks.
  # title/description/overview_df accepted for backward compatibility, ignored.
  wb <- createWorkbook()
  for (spec in sheet_specs) {
    df <- if (!is.null(spec$df)) spec$df else safe_read(spec$path)
    if (!is.null(df)) add_sheet(wb, spec$name, df)
  }
  saveWorkbook(wb, out_file, overwrite = TRUE)
  cat(sprintf("  Saved: %s (%.0f KB)\n\n", out_file, file.size(out_file) / 1e3))
}

# Cross-figure readers: prefer these over reading raw CSVs from another figure's c_data/.
read_sheet_df <- function(xlsx, sheet) {
  stopifnot("supplementary workbook missing" = file.exists(xlsx))
  as.data.frame(readxl::read_excel(xlsx, sheet = sheet))
}

read_matrix_sheet <- function(xlsx, sheet, row_col = "sample_id") {
  df <- read_sheet_df(xlsx, sheet)
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df[[row_col]]
  mat
}

read_vector_sheet <- function(xlsx, sheet) read_sheet_df(xlsx, sheet)[[1]]

matrix_to_df <- function(mat, row_col = "sample_id") {
  df <- data.frame(rn = rownames(mat),
                   as.data.frame(mat, check.names = FALSE),
                   check.names = FALSE, stringsAsFactors = FALSE)
  names(df)[1] <- row_col
  df
}

# Remove per-panel CSV intermediates after the workbook is written.
# preserve_patterns: regexes for paths that must not be deleted (upstream stages, shared cache).
cleanup_after_workbook <- function(sheet_specs,
                                    extra_subdirs = character(),
                                    extra_files = character(),
                                    preserve_patterns = c(
                                      "^00_input/",
                                      "^01_normalization/",
                                      "^02_Imputation/",
                                      "^03_DEP/",
                                      "^04_Figures/shared/")) {
  is_preserved <- function(path) {
    rel <- sub(paste0("^", here::here(), "/?"), "", path)
    any(vapply(preserve_patterns, function(p) grepl(p, path) || grepl(p, rel), logical(1)))
  }
  removed <- 0L; preserved <- 0L
  for (sp in sheet_specs) {
    if (!is.null(sp$path) && file.exists(sp$path)) {
      if (is_preserved(sp$path)) {
        preserved <- preserved + 1L
      } else {
        unlink(sp$path); removed <- removed + 1L
      }
    }
  }
  for (sd in extra_subdirs) {
    if (dir.exists(sd)) { unlink(sd, recursive = TRUE); removed <- removed + 1L }
  }
  for (f in extra_files) {
    if (file.exists(f) && !is_preserved(f)) { unlink(f); removed <- removed + 1L }
  }
  cat(sprintf("  cleanup: removed %d intermediate(s); preserved %d upstream/shared path(s)\n",
              removed, preserved))
}
