# R conventions

One page. The point is that any figure script reads the same way: load upstream
data with a canonical reader, transform with clearly-named objects, hand panels
to a composite. Keep it parsimonious.

## Script layout

- One panel builder per file under a figure's `a_script/panels/`.
- A `00_composite.R` sources the panels, runs the shared fit or load, stitches
  the panels, and writes the figure + workbook.
- Stats live with the panel they support, not in a separate stats script.

## Paths — `here::here`, always

Every path is built from the project root with `here::here("04_Figures", ...)`.

A bare relative path (`read_csv("c_data/x.csv")`) resolves against R's working
directory, which changes with how R was launched, so the same string works from
one entry point and fails from another. `here::here()` anchors to the project
root (found by the `.git` / `.Rproj` marker) and returns a stable absolute path.
No `setwd()`, no machine paths.

## Data in and out — canonical only

| Need | Use |
|---|---|
| Read/write a table | `readr::read_csv()` / `readr::write_csv()` |
| Read/write an R object (DAList, fit) | `readRDS()` / `saveRDS()` |
| Locate the file | `here::here(...)` |

Name the upstream source on the way in: `norm_long <- read_csv(here::here(...))`.

## Object names — specific nouns, shared suffixes

Name for what the object holds, then carry a consistent suffix so a reader tracks
shape and processing at a glance.

| Suffix | Means | Example |
|---|---|---|
| `_long` / `_wide` | tidy shape | `de_long`, `combined_wide` |
| `_mat` | matrix | `expr_mat` |
| `_norm` | normalized | `intensity_norm` |
| `_imp` | imputed | `expr_imp` |
| `_scores` / `_loadings` | PCA output | `pca_scores` |
| `meta` | sample metadata | `sample_meta` |

Never `df`, `data`, `input`, `result`, `temp`, `temp_df`, `output_df`. They tell
a reader nothing and are the first thing a reviewer flags.

## Comments — annotate the non-obvious, stay silent on the canonical

- One header line per script: what it produces.
- Annotate a **non-canonical or non-obvious** step with a concise, plain note on
  what it does and why. Hand-rolled permutations, unusual contrasts, a statistic
  the reader can't infer from the call — these earn a line.
- Do **not** restate a canonical or self-evident line. No `# read the csv`, no
  `# Step 1:`, no banners, no emoji.

```r
# good — the contrast is not obvious from the call
cm <- makeContrasts(MitoMain = (Mito + PHE_Mito) / 2 - (Ctl + PHE) / 2, levels = design)

# noise — the function already says this
norm_long <- read_csv(path)  # read the normalized data
```

## Helpers — tagged by the figures that use them

Every helper carries the figure span it serves, so the filename says where it is
used. Helpers live in `04_Figures/functions/`; a helper used by a single figure
lives in that figure's own `a_script/helpers/`.

- **One figure** → `F01_<topic>.R`, in that figure's `a_script/helpers/`.
- **A span of figures** → `F01-F03_<topic>.R`, in `functions/`.
- **All figures** → `F01-F06_<topic>.R`, in `functions/`.

Tag by the span of figures that source the helper. Three reported figures exist
(F01–F03), and each sources all five shared engines directly or transitively, so
they carry the full `F01-F03` span.

Current layout (applied):

| File | Location |
|---|---|
| `F01-F03_style_palettes_theme.R` | `functions/` |
| `F01-F03_data_paths_and_loaders.R` | `functions/` |
| `F01-F03_pathway_enrichment_dedup_ora.R` | `functions/` |
| `F01-F03_supplementary_workbook.R` | `functions/` |
| `F01-F03_composite_layout.R` | `functions/` |
| `F01_mitocarta_lens_lookup.R` | `F01_Proteome_Overview/a_script/helpers/` |

## Definition of done

- [ ] Paths through `here::here`; no `setwd`, no machine paths.
- [ ] `read_csv` / `write_csv` / `readRDS` / `saveRDS` for I/O.
- [ ] Objects named specifically with the suffix vocabulary; no `df` / `result`.
- [ ] Header line present; non-canonical steps annotated; nothing obvious restated.
- [ ] Stochastic steps seeded.
- [ ] Helpers tagged `shared_` or `F0x_`; outputs unchanged after the edit.
