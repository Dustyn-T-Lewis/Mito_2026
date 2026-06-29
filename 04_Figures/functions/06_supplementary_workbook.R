# 04_Figures/functions/06_supplementary_workbook.R
# Pipeline step 06: assemble each figure's single supplementary .xlsx.
# Every figure writes ONE workbook to its c_data/. The first sheet is always an
# Overview that documents, per data sheet: the sheet name, how it was used in the
# figure, and what its rows/columns contain. Result tables live only in the
# workbook (no redundant loose CSVs). This builder never deletes anything.

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Excel sheet name: <=31 chars, none of []:*?/\, unique within the workbook.
.sheet_name <- function(base, taken) {
  nm <- substr(gsub("[\\[\\]:*?/\\\\]", "_", base), 1, 31)
  i <- 2L
  while (nm %in% taken) {
    nm <- paste0(substr(nm, 1, 29), "_", i)
    i <- i + 1L
  }
  nm
}

.add_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data)
  addStyle(wb, name, createStyle(textDecoration = "bold"),
    rows = 1, cols = seq_len(max(1, ncol(data))), gridExpand = TRUE
  )
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(max(1, ncol(data))), widths = "auto")
  cat(sprintf("    + %s: %d x %d\n", name, nrow(data), ncol(data)))
}

# out_file     path to the .xlsx (its dirname is the figure's c_data/).
# figure_title one-line title shown atop the Overview sheet.
# sheet_specs  list of list(name=, df=, role=, contents=):
#                name     -> sheet label (sanitized + de-duplicated)
#                df       -> data.frame / tibble of results
#                role     -> short phrase: how this table was used in the figure
#                contents -> what the rows/columns mean
build_workbook <- function(out_file, figure_title, sheet_specs) {
  wb <- createWorkbook()

  # Reserve the Overview sheet in position 1; fill it after sheet names settle.
  addWorksheet(wb, "Overview")

  taken <- "Overview"
  ov <- data.frame(
    Sheet = character(0), `Role in figure` = character(0),
    Contents = character(0), check.names = FALSE
  )
  for (spec in sheet_specs) {
    df <- spec$df
    if (is.null(df) && !is.null(spec$path) && file.exists(spec$path)) {
      df <- as.data.frame(read_csv(spec$path, show_col_types = FALSE))
    }
    if (is.null(df)) {
      cat(sprintf("    SKIP (no data): %s\n", spec$name))
      next
    }
    nm <- .sheet_name(spec$name, taken)
    .add_sheet(wb, nm, df)
    taken <- c(taken, nm)
    ov <- rbind(ov, data.frame(
      Sheet = nm,
      `Role in figure` = spec$role %||% "",
      Contents = spec$contents %||% sprintf("%d rows x %d cols", nrow(df), ncol(df)),
      check.names = FALSE
    ))
  }

  # Overview: title row, then the sheet directory.
  writeData(wb, "Overview", figure_title, startRow = 1, startCol = 1)
  addStyle(wb, "Overview", createStyle(textDecoration = "bold", fontSize = 12),
    rows = 1, cols = 1
  )
  writeData(wb, "Overview", ov, startRow = 3, startCol = 1)
  addStyle(wb, "Overview", createStyle(textDecoration = "bold"),
    rows = 3, cols = seq_len(ncol(ov)), gridExpand = TRUE
  )
  setColWidths(wb, "Overview", cols = 1:3, widths = c(28, 46, 70))
  freezePane(wb, "Overview", firstActiveRow = 4)

  saveWorkbook(wb, out_file, overwrite = TRUE)
  cat(sprintf(
    "  Saved: %s (%.0f KB) — %d data sheets + Overview\n\n",
    out_file, file.size(out_file) / 1e3, length(taken) - 1L
  ))
}
