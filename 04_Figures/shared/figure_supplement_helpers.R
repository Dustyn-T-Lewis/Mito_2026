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
