# Figures v2 — Pathway bars + Cluster pilots Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

## ⚠ Post-implementation status (2026-06-20)

The implementation diverged from this plan in seven places during build-out;
the spec doc next to this plan
(`docs/superpowers/specs/2026-06-20-figures-v2-pathway-cluster-redesign.md`)
is now authoritative. Each plan task that diverged is flagged inline with an
`> [UPDATED]` callout. Summary of divergences (with commit hash that
introduced the change):

| # | Plan said | Code does | Commit |
|---|---|---|---|
| 1 | Pilots 1-3 use fixed `c = 6` | Per-pilot c via Dmin elbow + Mardia √(N/2) cap; observed c = 3 / 4 / 3 | `8392a4e` |
| 2 | WGCNA classifier `r_floor = 0.4` | `r_floor = 0.25` (n = 24 collapses every module to "Other" at 0.40); spec hedges with a "ranking aid, not inferential claim" note | `d7b2f3b`, `6e113d4` |
| 3 | Pilot 5 k-means `centers = 6` | k from gap statistic firstSEmax (Tibshirani 2001, B = 50); happens to return 6 for this dataset but is now data-driven | `8392a4e` |
| 4 | Cluster rows = "Hallmark only on main panels" | Hallmark **+ MitoCarta** side-by-side on every cluster row; workbook gains per-pilot `*_ora_mito` sheets (Reimand 2019 two-source pairing) | `d0eb731` |
| 5 | RRHO2 ranking = signed `-log10(P)*sign(logFC)` | Signed limma moderated-t (Smyth 2004); harmonizes with the fgsea cache | `ed28a1e` |
| 6 | No RRHO2 fallback for small quadrants | Top-20% percentile fallback when peak-overlap < 5 genes (UU + DD fire on this dataset; UD + DU use rrho2_peak); tagged in workbook (`source` column) + figure subtitle | `b47e2b5` |
| 7 | F06 workbook has a single `Overview` sheet | Auto-generated `Overview` (data dictionary) + driver-built `Pilot_summary` (per-pilot run summary); original collision between the two was the root of the "duplicate Overview" workbook bug | `b47e2b5` |
| 8 | WGCNA Interaction r computed from `c(Ctl=-1, Mito=+1, PHE=+1, PHE_Mito=-1)` | Sign-flipped to `c(Ctl=+1, Mito=-1, PHE=-1, PHE_Mito=+1)` to match the limma `(PHE_Mito − Mito) − (PHE − Ctl)` algebra | `b47e2b5` |

Beyond methodology divergences, the **file layout was reorganized post-audit**:

| File | Status |
|---|---|
| `functions/03_pathway_enrichment_dedup_ora.R` | Split into stub (sources 3 siblings) + `03a_dedup_engine.R` + `03b_enrichment_runners.R` + `03c_pathway_display.R`. Odds-ratio Fisher 2×2 calc deduplicated into `fora_odds_ratio()` in 03a (was 3 duplicate sites). Commit `0157237`. |
| `functions/05_volcano_ring_plot_builder.R` | Split into stub + 5 siblings (05a labels / 05b ring data / 05c volcano layers / 05d ring layers / 05e composite). |
| `06_Cluster/a_script/01_main_panels.R` | Pilots 4-6 extracted into `02_pilot_wgcna.R`, `03_pilot_logfc.R`, `04_pilot_rrho2.R`. Driver shrinks 738 → 248 lines. Commit `b289417`. |
| `06_Cluster/a_script/01_main_panels.R` | AI-tell pass: dropped 16 banner-divider lines + 4 `===` message banners + the redundant `Overview` collision. Commit `b47e2b5`. |

All five verifiers (`F04`, `F06 c-means`, `F06 WGCNA`, `F06 logfc`, `F06 RRHO2`) PASS on the post-cleanup state.

The task bodies below are the **original** step-by-step plan, preserved as
historical reference for understanding why specific decisions were made. They
are **no longer the source of truth** — read the spec next door for current
behavior, and the `> [UPDATED]` callouts inline for the specific lines that
were superseded.

---

**Goal:** Rebuild `04_Figures_v2/04_Pathway_bars` to recapitulate the Panel-D count summary across 5 DBs, and replace `04_Figures_v2/06_Cluster` with a 6-pilot framework (3 c-means significance gates + WGCNA + logFC-vector k-means + RRHO2) that uses a uniform per-cluster row layout (trajectory left | Hallmark ORA right).

**Architecture:** Two figure-script rewrites consuming existing `functions/` engines, plus one new shared helper file (`functions/07_cluster_row_layout.R`) carrying the per-cluster row builder, the significance-gate utility, and the WGCNA module loader. All upstream stages (00–03) and figure caches stay untouched. Each pilot saves its own standalone figure; a single supplementary workbook bundles all six pilots.

**Tech Stack:** R 4.x, `e1071::cmeans`, `WGCNA`, `RRHO2`, `fgsea`, `ggplot2`, `patchwork`, `scales`, `dplyr`/`tidyr`/`readr`, `here`, `openxlsx`.

## Global Constraints

- Reuse `04_Figures_v2/functions/` engines via `here::here()` source paths. Never re-derive style, gene sets, or fgsea cache.
- Data caches stay single-sourced from `04_Figures/shared/` (`fgsea_tstat_all_h9c2.csv`, `rat_gene_sets.rds`, `wgcna_network.rds` under `F05_modules/c_data/`).
- Inputs: `P05$imp_rds` = `DAList_imputed_missforest.rds`; `P05$comb` = `combined_results_pi.csv`.
- Contrasts: `CORE <- c("CTLvMITO","CTLvPHE","PHEvPHE_MITO","Interaction")`. Disease=`CTLvPHE`, Rescue=`PHEvPHE_MITO`.
- Thresholds: `H9C2_PI_THRESH = 0.05`, `H9C2_FDR_EXPLOR = 0.10`, p-value threshold = `0.05`.
- Style invariants: `FIG_THEME`, `PANEL_MD = 178`, `get_pdf_device()`, png `dpi = 300`, `units = "mm"`. Palettes: `H9C2_PAL_GROUP`, `H9C2_PAL_DIR`, `CONTRAST_COLORS`.
- Redundancy collapse: `0.5*overlap + 0.5*jaccard ≥ 0.375` (existing `deduplicate_enrichment()`).
- Seeds: `set.seed(42)` before any stochastic call (cmeans, kmeans, RRHO2, WGCNA-derived computations).
- No new `library()` calls in style helpers (already loaded); scripts use `suppressPackageStartupMessages()` around their library block.
- Naming: outputs are `MAIN_F04_*` or `MAIN_F06_<pilot_key>*`; supplementary diagnostics to `b_reports/supp/{pdf,png}/`.
- Commit message style (project rule): one short sentence, no AI/Claude/Anthropic mentions, no `Co-Authored-By` trailers.

---

## File Structure

| File | Role | Status |
|---|---|---|
| `04_Figures_v2/README.md` | Suite overview | **modify** (imp4p → missforest) |
| `04_Figures_v2/BUILD_PROMPT.md` | Regeneration spec | **modify** (imp4p → missforest) |
| `04_Figures_v2/functions/07_cluster_row_layout.R` | Shared row builder + gate + WGCNA loader + RRHO2 helper | **create** |
| `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R` | Panel-D count-summary script | **rewrite** |
| `04_Figures_v2/06_Cluster/a_script/01_main_panels.R` | 6-pilot driver | **rewrite** |

Outputs:
- `04_Pathway_bars/b_reports/main/{pdf,png}/MAIN_F04_pathway_bars.{pdf,png}` (~120×70 mm)
- `04_Pathway_bars/c_data/F04_supplementary.xlsx`
- `06_Cluster/b_reports/main/{pdf,png}/MAIN_F06_<pilot_key>.{pdf,png}` × 6 pilots
- `06_Cluster/b_reports/supp/{pdf,png}/MAIN_F06_<pilot_key>_selection.{pdf,png}` × 3 c-means pilots
- `06_Cluster/c_data/F06_supplementary.xlsx`

---

## Task 0: Documentation rot fix

**Files:**
- Modify: `04_Figures_v2/README.md`
- Modify: `04_Figures_v2/BUILD_PROMPT.md`

**Interfaces:**
- Consumes: nothing
- Produces: corrected docs

- [ ] **Step 1: Patch README.md**

Run:
```bash
grep -n imp4p /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026/04_Figures_v2/README.md
```
Expected: line(s) mentioning `DAList_imputed_imp4p.rds`.

Replace `DAList_imputed_imp4p.rds` with `DAList_imputed_missforest.rds` in README.md.

- [ ] **Step 2: Patch BUILD_PROMPT.md**

Run:
```bash
grep -n imp4p /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026/04_Figures_v2/BUILD_PROMPT.md
```
Expected: line(s) mentioning `DAList_imputed_imp4p.rds`.

Replace each `DAList_imputed_imp4p.rds` with `DAList_imputed_missforest.rds` in BUILD_PROMPT.md.

- [ ] **Step 3: Verify no `imp4p` references remain in v2 docs**

Run:
```bash
grep -rn imp4p /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026/04_Figures_v2/
```
Expected: no output.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/README.md 04_Figures_v2/BUILD_PROMPT.md
git commit -m "fix the imp4p reference in the v2 docs; figures actually load missForest"
```

---

## Task 1: Rebuild F04 as the Panel-D count summary

**Files:**
- Rewrite: `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R`
- Delete: `04_Figures_v2/04_Pathway_bars/c_data/shown_pathways.csv` (if present after rewrite)

**Interfaces:**
- Consumes: `04_Figures_v2/functions/02_data_paths_and_loaders.R`, `functions/03_pathway_enrichment_dedup_ora.R`, `functions/04_mitocarta_lens_lookup.R`, `functions/06_supplementary_workbook.R`; `04_Figures/shared/fgsea_tstat_all_h9c2.csv`; `04_Figures/shared/rat_gene_sets.rds`.
- Produces: `MAIN_F04_pathway_bars.{pdf,png}`, `F04_supplementary.xlsx` with sheets: `Overview`, `dep_pathway_counts`, `Transplant_Mito_sig_pathways`, `Disease_Phe_sig_pathways`, `Rescue_Mito+Phe_sig_pathways`, `Interaction_Mito_sig_pathways`.

- [ ] **Step 1: Write the verification script (tinytest-style)**

Create file `/tmp/verify_F04.R`:

```r
# Verification: F04 outputs exist with the expected shape.
suppressPackageStartupMessages({ library(openxlsx) })
root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
base <- file.path(root, "04_Figures_v2", "04_Pathway_bars")
png  <- file.path(base, "b_reports", "main", "png", "MAIN_F04_pathway_bars.png")
pdf  <- file.path(base, "b_reports", "main", "pdf", "MAIN_F04_pathway_bars.pdf")
xlsx <- file.path(base, "c_data",   "F04_supplementary.xlsx")
stopifnot("PNG missing" = file.exists(png), file.info(png)$size > 10000)
stopifnot("PDF missing" = file.exists(pdf), file.info(pdf)$size > 1000)
stopifnot("XLSX missing" = file.exists(xlsx))
sheets <- getSheetNames(xlsx)
expected <- c("Overview", "dep_pathway_counts",
              "Transplant_Mito_sig_pathways", "Disease_Phe_sig_pathways",
              "Rescue_Mito+Phe_sig_pathways", "Interaction_Mito_sig_pathways")
missing <- setdiff(expected, sheets)
if (length(missing)) stop("Missing sheets: ", paste(missing, collapse=", "))
cat("F04 verification PASS\n")
```

- [ ] **Step 2: Run verification — confirm FAIL (artifacts predate the rewrite)**

```bash
Rscript /tmp/verify_F04.R
```
Expected: FAIL on the sheet check (current xlsx has `pathway_bars` + per-contrast names but not the new schema).

- [ ] **Step 3: Rewrite the F04 script**

Overwrite `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R` with:

```r
#!/usr/bin/env Rscript
# F04 MAIN — Panel-D-style pathway count summary across 5 DBs.
# Per contrast, Up bar and Down bar drawn from 0; mito subset overdrawn darker.
# Pool = CANONICAL_DBS (Hallmark + Reactome + KEGG + MitoCarta + GO Slim);
# padj<0.05, size>=10, MITO_DROP_SETS excluded; cross-DB Jaccard+Overlap dedup
# at 0.375 per contrast (EnrichmentMap, Merico 2010 / Reimand 2019).
# Mito flag: database == MitoCarta OR pathway matches MITO_PATHWAY_REGEX.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(scales)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "04_mitocarta_lens_lookup.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE    <- here::here("04_Figures_v2", "04_Pathway_bars")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

CORE      <- c("CTLvMITO", "CTLvPHE", "PHEvPHE_MITO", "Interaction")
PADJ_CUT  <- 0.05
MIN_SIZE  <- 10
SIM_CUT   <- 0.375

fgsea_all     <- read_csv(here::here("04_Figures", "shared", "fgsea_tstat_all_h9c2.csv"),
                          show_col_types = FALSE)
rat_gene_sets <- readRDS(here::here("04_Figures", "shared", "rat_gene_sets.rds"))
set_pool      <- do.call(c, unname(rat_gene_sets[CANONICAL_DBS]))

# ---- Per-contrast: pool 5 DBs, dedup, classify mito ---------------------------
per_contrast <- function(ctr) {
  sig <- fgsea_all |>
    filter(contrast == ctr, database %in% CANONICAL_DBS,
           !is.na(padj), padj < PADJ_CUT,
           size >= MIN_SIZE, !pathway %in% MITO_DROP_SETS) |>
    arrange(padj)
  if (nrow(sig) > 1)
    sig <- deduplicate_enrichment(as.data.frame(sig), pathways = set_pool,
                                  cutoff = SIM_CUT, cross_db = TRUE) |> as_tibble()
  sig |>
    mutate(direction = if_else(NES > 0, "Up", "Down"),
           is_mito = database == "MitoCarta" |
                     grepl(MITO_PATHWAY_REGEX, pathway, perl = TRUE),
           contrast = ctr)
}
sig_pw <- bind_rows(lapply(CORE, per_contrast))
message(sprintf("F04: %d total post-dedup sig pathways across %d contrasts",
                nrow(sig_pw), length(CORE)))

# ---- Bar data: contrast × direction → total + mito subset --------------------
bar_df <- sig_pw |>
  summarise(total = n(), mito = sum(is_mito), .by = c(contrast, direction)) |>
  complete(contrast = CORE, direction = c("Up", "Down"),
           fill = list(total = 0L, mito = 0L)) |>
  mutate(contrast = factor(contrast, levels = CORE),
         direction = factor(direction, levels = c("Up", "Down")))

# Fills: light = total, dark = mito subset.
FILL_TOTAL <- c(Up = "#F4A582", Down = "#92C5DE")
FILL_MITO  <- c(Up = "#B2182B", Down = "#2166AC")

BAR_W   <- 0.34; GAP <- 0.06
centers <- setNames(seq_along(CORE), CORE)
bar_df <- bar_df |>
  mutate(x_center  = centers[as.character(contrast)] +
                     ifelse(direction == "Up", -(BAR_W/2 + GAP/2), (BAR_W/2 + GAP/2)),
         fill_tot  = FILL_TOTAL[as.character(direction)],
         fill_mito = FILL_MITO[as.character(direction)])

bg_df <- tibble(xmin = seq_along(CORE) - 0.5, xmax = seq_along(CORE) + 0.5,
                fill = unname(CONTRAST_COLORS[CORE]))

p <- ggplot() +
  geom_rect(data = bg_df, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = bg_df$fill, alpha = 0.16, color = "grey75", linewidth = 0.2) +
  geom_rect(data = bar_df,
            aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2,
                ymin = 0, ymax = total),
            fill = bar_df$fill_tot, color = "black", linewidth = 0.3) +
  geom_rect(data = filter(bar_df, mito > 0),
            aes(xmin = x_center - BAR_W/2, xmax = x_center + BAR_W/2,
                ymin = 0, ymax = mito),
            fill = filter(bar_df, mito > 0)$fill_mito,
            color = "black", linewidth = 0.3) +
  geom_text(data = filter(bar_df, total > 0),
            aes(x = x_center, y = total, label = total),
            vjust = -0.4, size = 2.4, fontface = "bold") +
  scale_x_continuous(breaks = seq_along(CORE),
                     labels = setNames(gsub("_", "\n", contrast_brief(CORE)), CORE),
                     expand = expansion(mult = 0)) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.08)),
               breaks = c(5, 25, 50, 100, 200, 300, 500)) +
  coord_cartesian(clip = "off") +
  labs(title    = "Pathway enrichment (Up / Down)",
       subtitle = "Hallmark + Reactome + KEGG + MitoCarta + GO Slim | dark = mito | sqrt y, cross-DB Jaccard-deduped",
       x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(axis.text.x = element_text(size = FIG_AXIS_TEXT, lineheight = 0.85, face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(2, 4, 2, 2))

# Inline 4-key legend (matches Panel D).
key_df <- tibble(y = 4:1,
                 lab  = c("Up total", "Up mito", "Down total", "Down mito"),
                 fill = c(FILL_TOTAL["Up"], FILL_MITO["Up"],
                          FILL_TOTAL["Down"], FILL_MITO["Down"]))
p_key <- ggplot(key_df) +
  geom_point(aes(0, y), shape = 22, size = 2, fill = key_df$fill,
             color = "grey30", stroke = 0.3) +
  geom_text(aes(0.22, y, label = lab), hjust = 0, size = 1.6, color = "grey20") +
  scale_x_continuous(limits = c(-0.1, 2.6)) +
  scale_y_continuous(limits = c(0.5, 4.5)) +
  theme_void()
p <- p + patchwork::inset_element(p_key, left = 0.02, right = 0.40,
                                  top = 0.99, bottom = 0.62)

FIG_W <- 120; FIG_H <- 70
ggsave(file.path(RPT_PDF, "MAIN_F04_pathway_bars.pdf"), p,
       width = FIG_W, height = FIG_H, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F04_pathway_bars.png"), p,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300, limitsize = FALSE)

# ---- Workbook ----------------------------------------------------------------
counts_tab <- bar_df |>
  transmute(contrast = contrast_brief(as.character(contrast)),
            direction, total, mito,
            mito_fraction = ifelse(total > 0, mito / total, NA_real_)) |>
  arrange(contrast, direction)

per_ctr_tabs <- lapply(CORE, function(ctr) {
  sig_pw |> filter(contrast == ctr) |>
    transmute(pathway, database, clean_label = clean_display_label(pathway),
              NES = round(NES, 3), padj = signif(padj, 3),
              size, direction, is_mito) |>
    arrange(padj)
})

sheet_specs <- c(
  list(list(name = "dep_pathway_counts", df = counts_tab,
            role     = "F04 bar heights — total significant pathways and mito subset per contrast x direction",
            contents = "contrast, direction (Up/Down), total significant pathways, mito subset, mito_fraction")),
  Map(function(ctr, tab) list(
        name = paste0(contrast_brief(ctr), "_sig_pathways"),
        df = tab,
        role     = sprintf("Post-dedup significant pathways behind the %s bars", contrast_brief(ctr)),
        contents = "pathway, database, clean_label, NES, padj, size, direction (Up/Down), is_mito"),
      CORE, per_ctr_tabs)
)

build_workbook(file.path(DAT, "F04_supplementary.xlsx"),
               figure_title = "F04 — Panel-D pathway count summary (5-DB pool, mito subset overlaid)",
               sheet_specs = sheet_specs)

# Remove the legacy shown_pathways.csv if a prior run left it behind — no consumer now.
old_csv <- file.path(DAT, "shown_pathways.csv")
if (file.exists(old_csv)) file.remove(old_csv)

message(sprintf("F04 done | %d sig pathways shown | output %s", nrow(sig_pw), RPT_PNG))
```

- [ ] **Step 4: Run the rewritten script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R
```
Expected: messages reporting per-contrast post-dedup sig counts, ending with `F04 done | N sig pathways shown`. Exit code 0.

- [ ] **Step 5: Run verification — confirm PASS**

```bash
Rscript /tmp/verify_F04.R
```
Expected: `F04 verification PASS`.

- [ ] **Step 6: Visual check**

Read `04_Figures_v2/04_Pathway_bars/b_reports/main/png/MAIN_F04_pathway_bars.png` and compare against `04_Figures/F01_QC_overview/b_reports/main/png/panels/MAIN_panel_F_enrichment.png`. Same structure: 4 contrasts × Up/Down dodged bars, mito subset overdrawn darker, sqrt y, contrast-tinted bands. The v2 totals are larger because the pool is 5 DBs vs the original 3.

- [ ] **Step 7: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/04_Pathway_bars/
git commit -m "rebuild the F04 pathway bars as the panel D count summary with mito subset overlaid"
```

---

## Task 2: Shared helpers (significance gate, WGCNA loader, row layout)

> **[UPDATED — see Post-implementation status #2, #4]** Helper signature
> `classify_module_sign_pattern(disease_r, rescue_r, r_floor = 0.4)` ships
> with default 0.4 but is **called from `02_pilot_wgcna.R` with `r_floor = 0.25`**.
> `build_cluster_row()` was extended to accept an `ora_plot2` argument (the
> MitoCarta column), and a sibling `run_mitocarta_ora()` helper was added next
> to `run_hallmark_ora()`. The shared odds-ratio block is now
> `fora_odds_ratio()` in `functions/03a_dedup_engine.R` (called by both
> `run_*_ora()` helpers in 07).

**Files:**
- Create: `04_Figures_v2/functions/07_cluster_row_layout.R`

**Interfaces:**
- Consumes: `04_Figures_v2/functions/01_style_palettes_theme.R` (already sourced by `02_data_paths_and_loaders.R`), `04_Figures_v2/functions/03_pathway_enrichment_dedup_ora.R`.
- Produces (all callable from pilot scripts):
  - `filter_sig_in_any_contrast(comb_long, col, threshold, contrasts, op = "lt")` → character vector of gene symbols.
  - `load_wgcna_modules(rds_path)` → `list(modules = tibble(gene, module), MEs = matrix samples × modules, color_lookup = named character)`.
  - `compute_me_contrast_correlations(MEs, meta, contrasts)` → `tibble(module, contrast, r, p)`; `meta` carries `Group` keyed to `MEs` rownames.
  - `classify_module_sign_pattern(disease_r, rescue_r, r_floor = 0.4)` → factor with levels `c("Reversal","Concordant up","Concordant down","Other")`.
  - `cluster_palette(n)` → named character vector length `n` of `viridis::turbo` hex colors.
  - `build_trajectory_panel(z_mat, cluster, x_levels, x_lab, color, kind = c("line","barlogfc"))` → ggplot.
  - `build_ora_bar_panel(ora_df, color, max_n = 6)` → ggplot.
  - `build_cluster_row(traj_plot, ora_plot, header_text, color, widths = c(1, 1.4))` → patchwork object.
  - `stack_cluster_rows(rows_list, title, subtitle)` → patchwork object (vertical stack).
  - `run_hallmark_ora(genes, universe)` → tibble (`pathway`, `padj`, `overlap`, `size`, `odds_ratio`) restricted to Hallmark only.

- [ ] **Step 1: Write smoke tests for the helpers (tinytest-style)**

Create `/tmp/test_helpers.R`:

```r
# Sanity checks on functions/07_cluster_row_layout.R helpers.
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(patchwork)
})

root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
source(file.path(root, "04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(file.path(root, "04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(file.path(root, "04_Figures_v2", "functions", "07_cluster_row_layout.R"))

# filter_sig_in_any_contrast: a single OR across contrasts.
fake <- tibble(
  gene = rep(c("A","B","C","D"), each = 2),
  contrast = rep(c("X","Y"), 4),
  pi_score = c(0.01, 0.5, 0.5, 0.5, 0.04, 0.6, 0.5, 0.5))
out <- filter_sig_in_any_contrast(fake, col = "pi_score",
                                  threshold = 0.05,
                                  contrasts = c("X","Y"))
stopifnot(setequal(out, c("A","C")))

# classify_module_sign_pattern: sign pairs.
lv <- levels(classify_module_sign_pattern(0, 0))
stopifnot(identical(lv,
  c("Reversal","Concordant up","Concordant down","Other")))
stopifnot(as.character(classify_module_sign_pattern( 0.6,-0.6)) == "Reversal")
stopifnot(as.character(classify_module_sign_pattern( 0.6, 0.6)) == "Concordant up")
stopifnot(as.character(classify_module_sign_pattern(-0.6,-0.6)) == "Concordant down")
stopifnot(as.character(classify_module_sign_pattern( 0.1, 0.1)) == "Other")

# cluster_palette
pal <- cluster_palette(4); stopifnot(length(pal) == 4, all(grepl("^#", pal)))

# WGCNA loader returns expected slots on the real artifact (smoke).
wpath <- file.path(root, "04_Figures", "F05_modules", "c_data", "wgcna_network.rds")
if (file.exists(wpath)) {
  w <- load_wgcna_modules(wpath)
  stopifnot(all(c("modules","MEs","color_lookup") %in% names(w)))
  stopifnot(nrow(w$modules) > 100)
  stopifnot(ncol(w$MEs) >= 2)
}

cat("Helper smoke tests PASS\n")
```

- [ ] **Step 2: Run smoke tests — confirm FAIL (helpers don't exist yet)**

```bash
Rscript /tmp/test_helpers.R
```
Expected: error sourcing `functions/07_cluster_row_layout.R` (file not found).

- [ ] **Step 3: Create the helper file**

Write `04_Figures_v2/functions/07_cluster_row_layout.R`:

```r
# 04_Figures_v2/functions/07_cluster_row_layout.R
# Shared helpers for the F06 pilot framework:
#   * filter_sig_in_any_contrast  — per-protein significance gate
#   * load_wgcna_modules          — slim wrapper around the F05 WGCNA artifact
#   * compute_me_contrast_correlations / classify_module_sign_pattern
#   * cluster_palette             — viridis::turbo-derived cluster colors
#   * build_trajectory_panel / build_ora_bar_panel / build_cluster_row
#   * stack_cluster_rows          — patchwork vertical stack
#   * run_hallmark_ora            — Hallmark-only ORA (fora) with no dedup

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(ggplot2); library(patchwork); library(scales); library(viridis)
})

filter_sig_in_any_contrast <- function(comb_long, col, threshold,
                                       contrasts, op = c("lt", "le")) {
  op <- match.arg(op)
  stopifnot(col %in% names(comb_long), "gene" %in% names(comb_long))
  d <- comb_long[comb_long$contrast %in% contrasts, , drop = FALSE]
  v <- d[[col]]
  hit <- if (op == "lt") !is.na(v) & v < threshold
         else            !is.na(v) & v <= threshold
  unique(d$gene[hit & !is.na(d$gene) & nzchar(d$gene)])
}

load_wgcna_modules <- function(rds_path) {
  obj <- readRDS(rds_path)
  # The build_wgcna_network.R artifact carries either a flat list of
  # (gene, module, MEs, ...) or a nested 'net' element; handle both.
  mods <- if (!is.null(obj$module_assignments)) obj$module_assignments
          else if (!is.null(obj$net) && !is.null(obj$net$colors))
            tibble(gene = obj$gene_ids %||% names(obj$net$colors),
                   module = unname(obj$net$colors))
          else stop("Unrecognised WGCNA artifact shape: cannot find module assignments")
  MEs <- if (!is.null(obj$MEs)) as.matrix(obj$MEs)
         else if (!is.null(obj$net) && !is.null(obj$net$MEs)) as.matrix(obj$net$MEs)
         else stop("Unrecognised WGCNA artifact shape: cannot find module eigengenes")
  # Module names are colors; the color lookup is the identity unless the
  # artifact already provides a translation.
  color_lookup <- setNames(unique(mods$module), unique(mods$module))
  list(modules = tibble::as_tibble(mods), MEs = MEs, color_lookup = color_lookup)
}

# Pearson r + Student p for each ME column against each binary contrast
# indicator. `meta` must carry a Group column matching MEs rownames; contrasts
# are interpreted as named lists list("Disease" = c("Ctl","PHE"), ...) where
# the second level is the +1 group and the first is the -1 group.
compute_me_contrast_correlations <- function(MEs, meta, contrasts) {
  stopifnot(nrow(MEs) == nrow(meta))
  out <- list()
  for (cn in names(contrasts)) {
    pair <- contrasts[[cn]]
    indic <- ifelse(meta$Group == pair[2],  1,
              ifelse(meta$Group == pair[1], -1, NA_real_))
    keep <- !is.na(indic)
    if (sum(keep) < 4) next
    for (mod in colnames(MEs)) {
      ct <- suppressWarnings(cor.test(MEs[keep, mod], indic[keep], method = "pearson"))
      out[[length(out) + 1]] <- tibble(module = mod, contrast = cn,
                                       r = unname(ct$estimate),
                                       p = ct$p.value)
    }
  }
  bind_rows(out)
}

classify_module_sign_pattern <- function(disease_r, rescue_r, r_floor = 0.4) {
  cls <- dplyr::case_when(
    abs(disease_r) >= r_floor & abs(rescue_r) >= r_floor &
      sign(disease_r) != sign(rescue_r)              ~ "Reversal",
    disease_r >=  r_floor & rescue_r >=  r_floor     ~ "Concordant up",
    disease_r <= -r_floor & rescue_r <= -r_floor     ~ "Concordant down",
    TRUE                                             ~ "Other"
  )
  factor(cls, levels = c("Reversal", "Concordant up", "Concordant down", "Other"))
}

cluster_palette <- function(n) {
  stopifnot(n >= 1)
  pal <- viridis::turbo(n + 2)[seq_len(n + 1)[-1]]   # trim extremes
  setNames(pal, as.character(seq_len(n)))
}

build_trajectory_panel <- function(z_mat, cluster, x_levels, x_lab,
                                   color, kind = c("line", "barlogfc")) {
  kind <- match.arg(kind)
  z_df <- as.data.frame(z_mat) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(cols = colnames(z_mat), names_to = "x", values_to = "expr") |>
    dplyr::mutate(x = factor(.data$x, levels = x_levels))
  if (kind == "line") {
    p <- ggplot(z_df, aes(x, expr, group = gene)) +
      geom_line(color = color, alpha = 0.30, linewidth = 0.25) +
      stat_summary(aes(group = 1), fun = mean, geom = "line",
                   color = color, linewidth = 0.9)
  } else {
    means <- z_df |> summarise(mean_expr = mean(.data$expr), .by = x)
    p <- ggplot(means, aes(x, mean_expr)) +
      geom_col(fill = color, color = "grey20", linewidth = 0.25, width = 0.7) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3)
  }
  p +
    labs(x = x_lab, y = if (kind == "line") "z" else "mean logFC") +
    FIG_THEME +
    theme(plot.margin = margin(2, 2, 2, 2),
          axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT))
}

build_ora_bar_panel <- function(ora_df, color, max_n = 6) {
  if (is.null(ora_df) || nrow(ora_df) == 0)
    return(ggplot() + annotate("text", 0, 0, label = "no Hallmark hits",
                               size = 2.4, color = "grey40") +
             theme_void())
  d <- ora_df |>
    arrange(.data$padj) |>
    head(max_n) |>
    mutate(label = clean_display_label(.data$pathway),
           neglog10 = -log10(.data$padj))
  ggplot(d, aes(reorder(.data$label, .data$neglog10), .data$neglog10)) +
    geom_col(fill = color, color = "grey20", linewidth = 0.2, width = 0.78) +
    coord_flip() +
    labs(x = NULL, y = "-log10 padj") +
    FIG_THEME +
    theme(plot.margin = margin(2, 4, 2, 2),
          axis.text.y = element_text(size = FIG_AXIS_TEXT - 0.5))
}

build_cluster_row <- function(traj_plot, ora_plot, header_text, color,
                              widths = c(1, 1.4)) {
  header <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = header_text,
             hjust = 0, vjust = 0.5, size = 2.3, color = color, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()
  body <- traj_plot + ora_plot + patchwork::plot_layout(widths = widths)
  header / body + patchwork::plot_layout(heights = c(0.10, 1))
}

stack_cluster_rows <- function(rows_list, title, subtitle) {
  if (length(rows_list) == 0)
    return(ggplot() + annotate("text", 0, 0, label = "no clusters") + theme_void())
  patchwork::wrap_plots(rows_list, ncol = 1) +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
      theme = theme(plot.title    = element_text(face = "bold", size = 8),
                    plot.subtitle = element_text(face = "italic", size = 5,
                                                 color = "grey30")))
}

# Hallmark-only ORA: fora over the Hallmark sublist of rat_gene_sets.rds. No
# dedup needed at the Hallmark level (50 sets, low redundancy).
run_hallmark_ora <- function(genes, universe,
                             rat_gene_sets_path = here::here("04_Figures", "shared",
                                                             "rat_gene_sets.rds")) {
  gs   <- readRDS(rat_gene_sets_path)
  hall <- gs$Hallmark
  if (length(hall) == 0) return(NULL)
  res <- fgsea::fora(pathways = hall, genes = genes, universe = universe,
                     minSize = 5, maxSize = 500)
  res <- as.data.frame(res)
  if (nrow(res) == 0) return(NULL)
  K <- length(intersect(genes, universe))
  N <- length(universe)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]; b <- K - a
    c <- res$size[i] - a; d <- N - K - c
    if (b <= 0 || c <= 0) Inf else (a * d) / (b * c)
  }, numeric(1))
  tibble::as_tibble(res)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
```

- [ ] **Step 4: Run smoke tests — confirm PASS**

```bash
Rscript /tmp/test_helpers.R
```
Expected: `Helper smoke tests PASS`.

- [ ] **Step 5: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/functions/07_cluster_row_layout.R
git commit -m "add shared helpers for the F06 pilot row layout"
```

---

## Task 3: F06 driver scaffold + Pilots 1–3 (c-means × p / Π / FDR gates)

> **[UPDATED — see Post-implementation status #1, #4, #7, code-organization rows]**
> - **c is per-pilot (Dmin elbow)**, not fixed 6. Calls `pick_c_dmin()` from
>   `07_cluster_row_layout.R` with `drop_frac = 0.10` and a `DEFAULT_C = 6L`
>   fallback. Observed: pilot_p c = 3, pilot_pi c = 4, pilot_fdr c = 3.
> - Cluster rows render **Hallmark + MitoCarta side-by-side** (extra ORA panel
>   per row); workbook gains `pilot_<key>_ora_mito` sheets per pilot.
> - **The first workbook spec is `name = "Pilot_summary"`, NOT `"Overview"`**
>   (the auto-generated Overview from `build_workbook()` was colliding with
>   the driver's hand-built Overview, producing the duplicate-`Overview_2`
>   workbook bug).
> - Driver file shrunk from 738 lines to **248 lines after pilots 4-6 were
>   extracted into siblings** (`02_pilot_wgcna.R`, `03_pilot_logfc.R`,
>   `04_pilot_rrho2.R` in the same `a_script/` directory, sourced after the
>   c-means block).
> - The script preamble comment claiming "Fixed c = 6" and `===` message
>   banners between pilots were dropped (commit `b47e2b5`).

**Files:**
- Rewrite: `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: helpers from Task 2; existing `e1071::cmeans`-based engine inlined from prior 06 script; `02_data_paths_and_loaders.R`; `03_pathway_enrichment_dedup_ora.R`; `06_supplementary_workbook.R`.
- Produces: `MAIN_F06_pilot_p.{pdf,png}`, `MAIN_F06_pilot_pi.{pdf,png}`, `MAIN_F06_pilot_fdr.{pdf,png}`; three corresponding `b_reports/supp/.../*_selection.{pdf,png}`; first three sheet groups in `F06_supplementary.xlsx`.

- [ ] **Step 1: Write verification script**

Create `/tmp/verify_F06_cmeans.R`:

```r
suppressPackageStartupMessages({ library(openxlsx) })
root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
base <- file.path(root, "04_Figures_v2", "06_Cluster")
for (k in c("pilot_p","pilot_pi","pilot_fdr")) {
  png <- file.path(base, "b_reports", "main", "png", paste0("MAIN_F06_", k, ".png"))
  stopifnot(paste("missing", png) = file.exists(png), file.info(png)$size > 10000)
  supp <- file.path(base, "b_reports", "supp", "png",
                    paste0("MAIN_F06_", k, "_selection.png"))
  stopifnot(paste("missing", supp) = file.exists(supp))
}
xlsx <- file.path(base, "c_data", "F06_supplementary.xlsx")
stopifnot(file.exists(xlsx))
sh <- getSheetNames(xlsx)
needed <- c("Overview",
            "pilot_p_membership","pilot_p_ora",
            "pilot_pi_membership","pilot_pi_ora",
            "pilot_fdr_membership","pilot_fdr_ora")
miss <- setdiff(needed, sh)
if (length(miss)) stop("missing sheets: ", paste(miss, collapse=", "))
cat("F06 c-means pilots verification PASS\n")
```

- [ ] **Step 2: Run verification — confirm FAIL**

```bash
Rscript /tmp/verify_F06_cmeans.R
```
Expected: FAIL (the new pilots don't exist yet).

- [ ] **Step 3: Rewrite the F06 driver with Pilots 1–3**

Overwrite `04_Figures_v2/06_Cluster/a_script/01_main_panels.R` with the driver below. It builds the c-means pilots only; later tasks add Pilots 4–6 by extending the `PILOTS` list and re-running.

```r
#!/usr/bin/env Rscript
# F06 MAIN — multi-pilot cluster framework. Each pilot is one figure with a
# per-cluster row layout (trajectory left | Hallmark ORA right, color-coded by
# cluster). Pilots 1-3 (this task): fuzzy c-means on group means, gated on
# p<0.05 / Pi<0.05 / FDR<0.10 in >=1 core contrast. Fixed c = 6, seed 42.
# Cluster-selection Dmin sweeps still emitted to b_reports/supp for record.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(readr); library(ggplot2); library(patchwork); library(scales)
  library(e1071); library(limma)
})

source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))
source(here::here("04_Figures_v2", "functions", "03_pathway_enrichment_dedup_ora.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))
source(here::here("04_Figures_v2", "functions", "07_cluster_row_layout.R"))

BASE     <- here::here("04_Figures_v2", "06_Cluster")
MAIN_PDF <- file.path(BASE, "b_reports", "main", "pdf")
MAIN_PNG <- file.path(BASE, "b_reports", "main", "png")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
DAT      <- file.path(BASE, "c_data")
for (d in c(MAIN_PDF, MAIN_PNG, SUPP_PDF, SUPP_PNG, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

FIG_W      <- PANEL_MD       # 178 mm
FIXED_C    <- 6L
SEED       <- 42L
CORE_MEMB  <- 0.5
C_RANGE    <- 2:12
SEEDS_DMIN <- 1:5

# Inputs: matrices + tables
dal     <- readRDS(P05$imp_rds)
prot    <- dal$data
ann     <- dal$annotation
meta    <- dal$metadata
gene_v  <- ann$gene[match(rownames(prot), ann$uniprot_id)]
keep_g  <- !is.na(gene_v) & nzchar(gene_v)
prot    <- prot[keep_g, , drop = FALSE]
gene_v  <- gene_v[keep_g]
gene_mat <- limma::avereps(prot, ID = gene_v)
ALL_GENES <- rownames(gene_mat)

grp <- meta$Group[match(colnames(gene_mat), meta$Col_ID)]
grp <- factor(grp, levels = H9C2_GROUP_LEVELS)
group_mat <- vapply(H9C2_GROUP_LEVELS, function(g)
  rowMeans(gene_mat[, grp == g, drop = FALSE]), numeric(nrow(gene_mat)))
rownames(group_mat) <- ALL_GENES

comb_long <- read_csv(P05$comb, show_col_types = FALSE)
CORE      <- H9C2_CORE_CONTRASTS

# c-means engine (carried over, simplified)
standardise_genes <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.finite(rowSums(z)), , drop = FALSE]
}
mestimate_fuzzifier <- function(z) {
  N <- nrow(z); D <- ncol(z)
  1 + (1418 / N + 22.05) * D^(-2) +
    (12.33 / N + 0.243) * D^(-0.0406 * log(N) - 0.1134)
}
mfuzz_cmeans <- function(z, c, m, seed = SEED) {
  set.seed(seed)
  e1071::cmeans(z, centers = c, m = m, method = "cmeans", iter.max = 200)
}
mean_dmin <- function(z, c, m, seeds = SEEDS_DMIN) {
  mean(vapply(seeds, function(s) min(dist(mfuzz_cmeans(z, c, m, seed = s)$centers)),
              numeric(1)))
}

# Per-pilot runner (c-means flavor)
run_cmeans_pilot <- function(key, gene_set, gate_label) {
  message(sprintf("\n=== %s (n_genes = %d) ===", key, length(gene_set)))
  mat <- group_mat[intersect(gene_set, rownames(group_mat)), , drop = FALSE]
  if (nrow(mat) < FIXED_C * 2)
    stop(sprintf("pilot %s: too few genes (%d) for c = %d", key, nrow(mat), FIXED_C))
  z <- standardise_genes(mat)
  m <- mestimate_fuzzifier(z)

  # selection sweep (diagnostic only; not used to pick c)
  dmin_tbl <- tibble(c = C_RANGE,
                     mean_Dmin = vapply(C_RANGE, function(c) mean_dmin(z, c, m),
                                        numeric(1)))
  fit  <- mfuzz_cmeans(z, FIXED_C, m, seed = SEED)
  hard <- fit$cluster
  max_mem <- apply(fit$membership, 1, max)
  memb <- tibble(gene = rownames(z), cluster = as.integer(hard),
                 membership = as.numeric(max_mem),
                 core = max_mem > CORE_MEMB)

  pal <- cluster_palette(FIXED_C)
  rows <- lapply(sort(unique(memb$cluster)), function(cl) {
    g_in_cl <- memb$gene[memb$cluster == cl]
    z_cl   <- z[g_in_cl, , drop = FALSE]
    ora    <- run_hallmark_ora(genes = memb$gene[memb$cluster == cl & memb$core],
                               universe = ALL_GENES)
    color  <- pal[as.character(cl)]
    hdr    <- sprintf("Cluster %d  |  n = %d (core %d)  |  Hallmark ORA",
                      cl, length(g_in_cl), sum(memb$cluster == cl & memb$core))
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl, cluster = cl,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(ora, color = color, max_n = 6),
      header_text = hdr, color = color)
  })

  fig <- stack_cluster_rows(rows,
    title    = sprintf("F06 %s — c = %d fuzzy c-means", key, FIXED_C),
    subtitle = sprintf("gate: %s in >=1 of {%s}; m = %.3f; rows = clusters; right = Hallmark top-6",
                       gate_label, paste(CORE, collapse = " / "), m))
  h_mm <- 32 + 32 * FIXED_C
  ggsave(file.path(MAIN_PDF, sprintf("MAIN_F06_%s.pdf", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, sprintf("MAIN_F06_%s.png", key)), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp: Dmin-vs-c diagnostic
  sup <- ggplot(dmin_tbl, aes(c, mean_Dmin)) +
    geom_line(color = "grey50", linewidth = 0.4) +
    geom_point(size = 1.4, color = "grey30") +
    geom_vline(xintercept = FIXED_C, color = "#D6604D",
               linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = C_RANGE) +
    labs(title = sprintf("F06 %s cluster-selection diagnostic", key),
         subtitle = sprintf("fixed c = %d (no auto-pick); m = %.3f", FIXED_C, m),
         x = "number of clusters (c)",
         y = "mean min centroid distance (Dmin)") + FIG_THEME
  ggsave(file.path(SUPP_PDF, sprintf("MAIN_F06_%s_selection.pdf", key)), sup,
         width = 120, height = 80, units = "mm", device = pdf_dev)
  ggsave(file.path(SUPP_PNG, sprintf("MAIN_F06_%s_selection.png", key)), sup,
         width = 120, height = 80, units = "mm", dpi = 300)

  # ORA across all clusters (for the workbook)
  ora_all <- bind_rows(lapply(sort(unique(memb$cluster)), function(cl) {
    g <- memb$gene[memb$cluster == cl & memb$core]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, cluster = cl)
  }))

  list(key = key, m = m, fixed_c = FIXED_C, n_genes = nrow(z),
       memb = memb, ora = ora_all, dmin_tbl = dmin_tbl)
}

# PILOTS 1-3 (c-means × p / Π / FDR)
PILOTS_CMEANS <- list(
  list(key = "pilot_p",   col = "P.Value",   threshold = 0.05,
       gate_label = "p < 0.05"),
  list(key = "pilot_pi",  col = "pi_score",  threshold = 0.05,
       gate_label = "Π < 0.05"),
  list(key = "pilot_fdr", col = "adj.P.Val", threshold = 0.10,
       gate_label = "FDR < 0.10")
)

results <- list()
for (p in PILOTS_CMEANS) {
  genes <- filter_sig_in_any_contrast(comb_long, col = p$col,
                                      threshold = p$threshold, contrasts = CORE)
  results[[p$key]] <- run_cmeans_pilot(key = p$key, gene_set = genes,
                                       gate_label = p$gate_label)
}

# Single supplementary workbook (Pilots 1-3 only at this stage)
overview <- tibble(
  Pilot = vapply(results, `[[`, character(1), "key"),
  Method = "fuzzy c-means on group means",
  Gate = vapply(PILOTS_CMEANS, `[[`, character(1), "gate_label"),
  N_genes = vapply(results, `[[`, integer(1), "n_genes"),
  Fuzzifier_m = vapply(results, function(r) round(r$m, 3), numeric(1)),
  Cluster_c = vapply(results, `[[`, integer(1), "fixed_c"))

sheet_specs <- list(list(
  name = "Overview", df = overview,
  role = "F06 pilots — run-level summary",
  contents = "pilot key, method, gate, gene count, fuzzifier m, cluster count c"))

for (r in results) {
  sheet_specs <- c(sheet_specs, list(
    list(name = paste0(r$key, "_membership"), df = r$memb,
         role = sprintf("Soft-cluster assignment for %s", r$key),
         contents = "gene, hard cluster, max membership, core flag (>0.5)"),
    list(name = paste0(r$key, "_ora"),
         df = if (is.null(r$ora)) tibble(cluster = integer(), pathway = character(),
                                         padj = numeric(), overlap = integer(),
                                         size = integer(), odds_ratio = numeric())
              else r$ora,
         role = sprintf("Hallmark ORA per cluster for %s", r$key),
         contents = "cluster, pathway, padj, overlap, size, odds_ratio")))
}

build_workbook(file.path(DAT, "F06_supplementary.xlsx"),
               figure_title = "F06 — multi-pilot cluster framework",
               sheet_specs = sheet_specs)

message(sprintf("F06 c-means pilots done: %s", paste(names(results), collapse = ", ")))
```

- [ ] **Step 4: Run the script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R
```
Expected: per-pilot run messages, ending with `F06 c-means pilots done: pilot_p, pilot_pi, pilot_fdr`. Exit code 0. Warnings about pdf device fallback (cairo → quartz) are acceptable.

- [ ] **Step 5: Run verification**

```bash
Rscript /tmp/verify_F06_cmeans.R
```
Expected: `F06 c-means pilots verification PASS`.

- [ ] **Step 6: Visual check (one pilot)**

Read `04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_pilot_pi.png`. Confirm 6 cluster rows, each with a trajectory line on the left and Hallmark bars on the right, both in the same cluster color.

- [ ] **Step 7: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/06_Cluster/
git commit -m "replace the F06 cluster figure with the three c-means significance-gate pilots"
```

---

## Task 4: Pilot 4 — WGCNA modules

> **[UPDATED — see Post-implementation status #2, #4, #8]**
> - Pilot lives in `06_Cluster/a_script/02_pilot_wgcna.R` (sourced by the F06
>   driver after pilots 1-3); it is NOT inline in the driver.
> - **`r_floor = 0.25` not 0.40** when calling `classify_module_sign_pattern()`.
>   At n = 24, the 0.40 floor collapses every module to "Other". The pilot
>   carries an in-code comment hedging that the column is a "ranking aid, not
>   an inferential claim" — propagate that hedge into any figure caption.
> - **Interaction r vector** is `c(Ctl = +1, Mito = -1, PHE = -1, PHE_Mito = +1)`
>   (the limma algebra `(PHE_Mito − Mito) − (PHE − Ctl)` expands to those signs;
>   the original plan had the signs reversed).
> - Cluster rows render **Hallmark + MitoCarta**; workbook gains
>   `pilot_wgcna_ora_mito` sheet.

**Files:**
- Modify: `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: `load_wgcna_modules()`, `compute_me_contrast_correlations()`, `classify_module_sign_pattern()` from Task 2; `wgcna_network.rds` from `04_Figures/F05_modules/c_data/`.
- Produces: `MAIN_F06_pilot_wgcna.{pdf,png}`; supp ME-trait heatmap `MAIN_F06_pilot_wgcna_me_traits.{pdf,png}`; workbook sheets `pilot_wgcna_membership`, `pilot_wgcna_me_traits`, `pilot_wgcna_ora`.

- [ ] **Step 1: Inspect WGCNA artifact shape (read-only)**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript -e 'w <- readRDS("04_Figures/F05_modules/c_data/wgcna_network.rds"); str(w, max.level=2)'
```
Expected: a list with elements likely including `module_assignments` (tibble of gene→module) and `MEs` (samples × modules matrix). If the names differ, adjust `load_wgcna_modules()` in `functions/07_cluster_row_layout.R` to match before continuing — keep its return contract identical.

- [ ] **Step 2: Extend verification script**

Append to `/tmp/verify_F06_cmeans.R` (save as `/tmp/verify_F06_wgcna.R`):

```r
# Inherits the c-means checks; adds WGCNA pilot.
source("/tmp/verify_F06_cmeans.R")
root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
base <- file.path(root, "04_Figures_v2", "06_Cluster")
suppressPackageStartupMessages({ library(openxlsx) })
png <- file.path(base, "b_reports", "main", "png", "MAIN_F06_pilot_wgcna.png")
stopifnot("pilot_wgcna png missing" = file.exists(png), file.info(png)$size > 10000)
supp <- file.path(base, "b_reports", "supp", "png", "MAIN_F06_pilot_wgcna_me_traits.png")
stopifnot("pilot_wgcna me_traits supp missing" = file.exists(supp))
sh <- getSheetNames(file.path(base, "c_data", "F06_supplementary.xlsx"))
need <- c("pilot_wgcna_membership","pilot_wgcna_me_traits","pilot_wgcna_ora")
miss <- setdiff(need, sh)
if (length(miss)) stop("missing wgcna sheets: ", paste(miss, collapse=", "))
cat("F06 WGCNA pilot verification PASS\n")
```

- [ ] **Step 3: Run — confirm FAIL**

```bash
Rscript /tmp/verify_F06_wgcna.R
```
Expected: FAIL on the wgcna png check.

- [ ] **Step 4: Add the WGCNA pilot to the F06 driver**

In `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`, append the following block **before** the `build_workbook(...)` call (so the workbook picks the new sheets up):

```r
# PILOT 4 — WGCNA modules
WGCNA_RDS <- here::here("04_Figures", "F05_modules", "c_data", "wgcna_network.rds")
if (file.exists(WGCNA_RDS)) {
  message("\n=== pilot_wgcna ===")
  w <- load_wgcna_modules(WGCNA_RDS)
  mods <- w$modules |> filter(.data$module != "grey")
  MEs  <- w$MEs

  # Build the per-sample (Group) indicator vectors over the 4 contrasts the
  # eigengenes are correlated against. Group ordering must match MEs rownames.
  me_meta <- tibble(Col_ID = rownames(MEs)) |>
    left_join(as_tibble(meta) |> select(Col_ID, Group), by = "Col_ID")
  contrast_pairs <- list(
    Disease    = c("Ctl", "PHE"),
    Transplant = c("Ctl", "Mito"),
    Rescue     = c("PHE", "PHE_Mito"))
  me_corr <- compute_me_contrast_correlations(MEs, me_meta, contrast_pairs)

  signs <- me_corr |>
    pivot_wider(id_cols = module, names_from = contrast, values_from = r,
                names_prefix = "r_") |>
    mutate(sign_pattern = classify_module_sign_pattern(r_Disease, r_Rescue))
  # Row order: Reversal first, then Concordant up/down, then Other, then by |r_Rescue|
  mod_order <- signs |>
    arrange(.data$sign_pattern, desc(abs(.data$r_Rescue))) |>
    pull(.data$module)

  pal <- setNames(viridis::turbo(length(mod_order) + 2)[seq_along(mod_order) + 1],
                  mod_order)

  rows <- lapply(mod_order, function(mod) {
    g_in <- mods$gene[mods$module == mod]
    g_in <- intersect(g_in, rownames(gene_mat))
    if (length(g_in) < 5) return(NULL)
    # trajectory on group means (z-scored per gene over conditions)
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal[mod]
    sp   <- signs$sign_pattern[match(mod, signs$module)]
    rD   <- signs$r_Disease[match(mod, signs$module)]
    rR   <- signs$r_Rescue [match(mod, signs$module)]
    rT   <- signs$r_Transplant[match(mod, signs$module)]
    hdr  <- sprintf("Module %s  |  n = %d  |  %s  |  r(D)=%.2f r(R)=%.2f r(T)=%.2f",
                    mod, length(g_in), sp, rD, rR, rT)
    ora  <- run_hallmark_ora(g_in, universe = ALL_GENES)
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl, cluster = mod,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(ora, color = color, max_n = 6),
      header_text = hdr, color = color)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fig <- stack_cluster_rows(rows,
    title    = "F06 pilot_wgcna — modules from F05 WGCNA artifact",
    subtitle = "rows = modules; ordered by Disease<->Rescue sign pattern (reversal first); right = Hallmark top-6")
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_wgcna.pdf"), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_wgcna.png"), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # supp ME-trait heatmap
  me_hm <- ggplot(me_corr, aes(.data$contrast, .data$module, fill = .data$r)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.6, color = "grey15") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(title = "F06 pilot_wgcna — module eigengene × contrast indicator",
         x = NULL, y = NULL) + FIG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = FIG_AXIS_TEXT))
  ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_wgcna_me_traits.pdf"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_wgcna_me_traits.png"), me_hm,
         width = 100, height = 6 + 4 * length(unique(me_corr$module)),
         units = "mm", dpi = 300, limitsize = FALSE)

  # ORA across all modules for the workbook
  ora_w <- bind_rows(lapply(mod_order, function(mod) {
    g <- mods$gene[mods$module == mod]
    o <- run_hallmark_ora(g, universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, module = mod)
  }))

  results$pilot_wgcna <- list(
    key = "pilot_wgcna",
    sheets = list(
      pilot_wgcna_membership = mods,
      pilot_wgcna_me_traits  = signs,
      pilot_wgcna_ora        = if (is.null(ora_w)) tibble() else ora_w))
} else {
  warning("WGCNA artifact not found at ", WGCNA_RDS, " — skipping pilot_wgcna")
}
```

And **extend the `sheet_specs` loop** in the same script (just before `build_workbook(...)`) to pick up the wgcna sheets:

```r
if (!is.null(results$pilot_wgcna)) {
  wg <- results$pilot_wgcna$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_wgcna_membership", df = wg$pilot_wgcna_membership,
              role = "Gene -> WGCNA module assignment (grey excluded)",
              contents = "gene, module (color label from F05 build)")),
    list(list(name = "pilot_wgcna_me_traits", df = wg$pilot_wgcna_me_traits,
              role = "Module eigengene Pearson r vs contrast indicator vectors",
              contents = "module, r_Disease, r_Transplant, r_Rescue, sign_pattern")),
    list(list(name = "pilot_wgcna_ora", df = wg$pilot_wgcna_ora,
              role = "Hallmark ORA per module",
              contents = "module, pathway, padj, overlap, size, odds_ratio")))
  # Add a row to Overview
  overview2 <- tibble(
    Pilot = "pilot_wgcna", Method = "WGCNA modules (existing artifact)",
    Gate = "all genes (network gate)",
    N_genes = nrow(wg$pilot_wgcna_membership),
    Fuzzifier_m = NA_real_,
    Cluster_c = length(unique(wg$pilot_wgcna_membership$module)))
  sheet_specs[[1]]$df <- bind_rows(sheet_specs[[1]]$df, overview2)
}
```

- [ ] **Step 5: Run the script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R
```
Expected: c-means messages followed by `=== pilot_wgcna ===` and a successful render; exit 0.

- [ ] **Step 6: Run verification**

```bash
Rscript /tmp/verify_F06_wgcna.R
```
Expected: `F06 WGCNA pilot verification PASS`.

- [ ] **Step 7: Visual check**

Read `04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_pilot_wgcna.png`. Confirm rows are labeled with `Module <color>  |  n = …  |  Reversal/Concordant up/...` and the first rows are reversal modules (opposite-sign r_Disease and r_Rescue).

- [ ] **Step 8: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/06_Cluster/
git commit -m "add the WGCNA module pilot for F06, ordered by disease-rescue reversal pattern"
```

---

## Task 5: Pilot 5 — k-means on per-protein 4-D logFC vector

> **[UPDATED — see Post-implementation status #3, #4]**
> - Pilot lives in `06_Cluster/a_script/03_pilot_logfc.R` (sourced by the F06
>   driver after pilot 4); it is NOT inline in the driver.
> - **k chosen by gap statistic firstSEmax** (Tibshirani, Walther & Hastie 2001,
>   `cluster::clusGap(B = 50)`) via `pick_c_gap()` from
>   `07_cluster_row_layout.R`. Returns k = 6 for this dataset (same as the
>   plan's hard-coded value, but now data-driven).
> - **Contrast vector uses canonical `H9C2_CORE_CONTRASTS` order**:
>   `(logFC_CTLvMITO, logFC_CTLvPHE, logFC_PHEvPHE_MITO, logFC_Interaction)`.
>   k-means is order-invariant so this doesn't change cluster geometry.
> - Quadrant labeler still keys off `logFC_CTLvPHE` (Disease) and
>   `logFC_PHEvPHE_MITO` (Rescue) regardless of vector position.
> - Cluster rows render **Hallmark + MitoCarta**; workbook gains
>   `pilot_logfc_ora_mito` sheet.
> - Supp diagnostic is now the **gap(k) ± 1 SE bar chart** (not silhouette or
>   within-cluster SS as the plan body suggests).

**Files:**
- Modify: `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: `load_combined_wide()`; `cluster_palette()`, `build_*_panel()`, `build_cluster_row()`, `stack_cluster_rows()`, `run_hallmark_ora()` from Task 2.
- Produces: `MAIN_F06_pilot_logfc.{pdf,png}`; supp `MAIN_F06_pilot_logfc_selection.{pdf,png}` (within-SS elbow); workbook sheets `pilot_logfc_membership`, `pilot_logfc_centroids`, `pilot_logfc_ora`.

- [ ] **Step 1: Extend verification script**

Save as `/tmp/verify_F06_logfc.R`:

```r
source("/tmp/verify_F06_wgcna.R")
root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
base <- file.path(root, "04_Figures_v2", "06_Cluster")
suppressPackageStartupMessages({ library(openxlsx) })
png <- file.path(base, "b_reports", "main", "png", "MAIN_F06_pilot_logfc.png")
stopifnot("pilot_logfc png missing" = file.exists(png), file.info(png)$size > 10000)
supp <- file.path(base, "b_reports", "supp", "png", "MAIN_F06_pilot_logfc_selection.png")
stopifnot("pilot_logfc selection supp missing" = file.exists(supp))
sh <- getSheetNames(file.path(base, "c_data", "F06_supplementary.xlsx"))
need <- c("pilot_logfc_membership","pilot_logfc_centroids","pilot_logfc_ora")
miss <- setdiff(need, sh)
if (length(miss)) stop("missing logfc sheets: ", paste(miss, collapse=", "))
cat("F06 logfc pilot verification PASS\n")
```

- [ ] **Step 2: Run — confirm FAIL**

```bash
Rscript /tmp/verify_F06_logfc.R
```
Expected: FAIL on the logfc png check.

- [ ] **Step 3: Append Pilot 5 to the F06 driver**

Append before the WGCNA block in `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`:

```r
# PILOT 5 — k-means on per-protein 4-D logFC vector
message("\n=== pilot_logfc ===")
comb_wide <- load_combined_wide()
lf_cols  <- paste0("logFC_", CORE)
stopifnot(all(lf_cols %in% names(comb_wide)))
lf <- comb_wide |>
  select(gene, all_of(lf_cols)) |>
  filter(!is.na(.data$gene), nzchar(.data$gene),
         if_all(all_of(lf_cols), ~ !is.na(.x)))
lf_mat <- as.matrix(lf[, lf_cols])
rownames(lf_mat) <- lf$gene
set.seed(SEED)
km <- kmeans(lf_mat, centers = FIXED_C, nstart = 50, iter.max = 100)

# Quadrant-style label from (logFC_Disease, logFC_Rescue) centroid signs.
disease_idx <- which(lf_cols == "logFC_CTLvPHE")
rescue_idx  <- which(lf_cols == "logFC_PHEvPHE_MITO")
quadrant_label <- function(dis, res) {
  if (abs(dis) < 0.1 & abs(res) < 0.1) return("Neutral")
  if (dis > 0 & res < 0) return("Reversed Down")     # disease up, rescue brings down
  if (dis < 0 & res > 0) return("Reversed Up")
  if (dis > 0 & res > 0) return("Concordant Up")
  if (dis < 0 & res < 0) return("Concordant Down")
  "Other"
}
centroids <- as_tibble(km$centers, rownames = "cluster") |>
  mutate(label = mapply(quadrant_label, km$centers[, disease_idx],
                                         km$centers[, rescue_idx]),
         n = unname(table(km$cluster))[as.integer(cluster)])

pal <- cluster_palette(FIXED_C)
rows <- lapply(seq_len(FIXED_C), function(cl) {
  g_in <- rownames(lf_mat)[km$cluster == cl]
  color <- pal[as.character(cl)]
  # trajectory = mean-logFC bar per contrast (kind = "barlogfc")
  z_cl  <- lf_mat[g_in, , drop = FALSE]   # already in logFC space
  colnames(z_cl) <- gsub("^logFC_", "", colnames(z_cl))
  ora   <- run_hallmark_ora(g_in, universe = unique(rownames(lf_mat)))
  cent  <- centroids[centroids$cluster == as.character(cl), ]
  hdr   <- sprintf("Cluster %d  |  n = %d  |  %s",
                   cl, cent$n, cent$label)
  build_cluster_row(
    traj_plot = build_trajectory_panel(z_cl, cluster = cl,
                                       x_levels = colnames(z_cl),
                                       x_lab = "contrast (centroid logFC)",
                                       color = color, kind = "barlogfc"),
    ora_plot  = build_ora_bar_panel(ora, color = color, max_n = 6),
    header_text = hdr, color = color)
})

fig <- stack_cluster_rows(rows,
  title    = sprintf("F06 pilot_logfc — k-means on per-protein 4-D logFC vector (c = %d)", FIXED_C),
  subtitle = "rows = clusters; cluster labels derived from (Disease, Rescue) sign quadrant; right = Hallmark top-6")
h_mm <- 32 + 32 * FIXED_C
ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_logfc.pdf"), fig,
       width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_logfc.png"), fig,
       width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

# supp: WSS elbow
elbow_tbl <- tibble(c = 2:10,
                    tot_wss = vapply(2:10, function(k) {
                      set.seed(SEED)
                      kmeans(lf_mat, centers = k, nstart = 25, iter.max = 100)$tot.withinss
                    }, numeric(1)))
sup <- ggplot(elbow_tbl, aes(c, tot_wss)) +
  geom_line(color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.4, color = "grey30") +
  geom_vline(xintercept = FIXED_C, color = "#D6604D",
             linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(breaks = 2:10) +
  labs(title = "F06 pilot_logfc cluster-selection diagnostic",
       subtitle = sprintf("fixed c = %d (no auto-pick); within-cluster SS by k", FIXED_C),
       x = "number of clusters (c)", y = "total within-cluster SS") + FIG_THEME
ggsave(file.path(SUPP_PDF, "MAIN_F06_pilot_logfc_selection.pdf"), sup,
       width = 120, height = 80, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "MAIN_F06_pilot_logfc_selection.png"), sup,
       width = 120, height = 80, units = "mm", dpi = 300)

memb_lf <- tibble(gene = rownames(lf_mat), cluster = as.integer(km$cluster))
ora_lf <- bind_rows(lapply(seq_len(FIXED_C), function(cl) {
  g <- memb_lf$gene[memb_lf$cluster == cl]
  o <- run_hallmark_ora(g, universe = unique(rownames(lf_mat)))
  if (is.null(o) || nrow(o) == 0) return(NULL)
  mutate(o, cluster = cl)
}))

results$pilot_logfc <- list(
  key = "pilot_logfc",
  sheets = list(
    pilot_logfc_membership = memb_lf,
    pilot_logfc_centroids  = centroids,
    pilot_logfc_ora        = if (is.null(ora_lf)) tibble() else ora_lf))
```

And **extend the sheet-spec assembly** (just before `build_workbook(...)`) to also handle `pilot_logfc`:

```r
if (!is.null(results$pilot_logfc)) {
  lg <- results$pilot_logfc$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_logfc_membership", df = lg$pilot_logfc_membership,
              role = "Gene -> k-means cluster on per-protein 4-D logFC vector",
              contents = "gene, cluster (1..c)")),
    list(list(name = "pilot_logfc_centroids", df = lg$pilot_logfc_centroids,
              role = "k-means centroids in logFC space + quadrant label",
              contents = "cluster, logFC_<contrast>... mean centroid, quadrant label, n")),
    list(list(name = "pilot_logfc_ora", df = lg$pilot_logfc_ora,
              role = "Hallmark ORA per logFC cluster",
              contents = "cluster, pathway, padj, overlap, size, odds_ratio")))
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(Pilot = "pilot_logfc", Method = "k-means on per-protein 4-D logFC vector",
           Gate = "all genes with non-NA logFC across core contrasts",
           N_genes = nrow(lg$pilot_logfc_membership),
           Fuzzifier_m = NA_real_, Cluster_c = FIXED_C))
}
```

- [ ] **Step 4: Run + verify**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R
Rscript /tmp/verify_F06_logfc.R
```
Expected: `F06 logfc pilot verification PASS`.

- [ ] **Step 5: Visual check**

Read `04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_pilot_logfc.png`. Each row's left panel is a bar of mean logFC per contrast (not a line), and the header says e.g. `Cluster 3  |  n = 612  |  Reversed Down`.

- [ ] **Step 6: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/06_Cluster/
git commit -m "add the logFC-vector k-means pilot for F06 with quadrant-labeled clusters"
```

---

## Task 6: Pilot 6 — RRHO2 threshold-free Disease ↔ Rescue concordance map

> **[UPDATED — see Post-implementation status #4, #5, #6]**
> - Pilot lives in `06_Cluster/a_script/04_pilot_rrho2.R` (sourced by the F06
>   driver after pilot 5); it is NOT inline in the driver.
> - **Ranking statistic = signed limma moderated-t** (`cw$t_CTLvPHE`,
>   `cw$t_PHEvPHE_MITO`), NOT the plan's original `-log10(P)*sign(logFC)`.
>   Variance-stabilized at n = 24 (Smyth 2004) and consistent with fgsea
>   upstream. RRHO2's hypergeometric tail depends on rank order, not metric
>   value, so the two are essentially equivalent on signal quadrants but the
>   mod-t form is better-behaved on the small-overlap tail.
> - **Top-20% percentile fallback** for sparse quadrants: when the RRHO2
>   peak-overlap set is < 5 genes, the quadrant gene set is filled from the
>   intersection of the top-20% rank fraction of both lists in the quadrant's
>   implied direction. The workbook tags each gene with `source ∈
>   {rrho2_peak, pct_fallback}` and the figure subtitle annotates which
>   quadrants used the fallback. **The fallback yields a biologically
>   meaningful set but is NOT RRHO2-significant** — manuscript captions must
>   inherit this hedge. For the current dataset: UU + DD fire on
>   pct_fallback (35 + 36 genes), UD + DU on rrho2_peak (1524 + 1160).
> - Per-quadrant rows render **Hallmark + MitoCarta**; workbook gains
>   `pilot_rrho2_ora_mito` sheet.
> - The heatmap saves to a **separate file** `MAIN_F06_pilot_rrho2_heatmap.{pdf,png}`
>   (base-graphics output that cannot share a patchwork canvas with the
>   ggplot per-quadrant rows). The PILOT_KEYS allowlist in the driver
>   includes `pilot_rrho2_heatmap` to preserve it during orphan cleanup.
> - Quadrant gene lists ship in a single long-format sheet
>   `pilot_rrho2_genelists` (columns: `quadrant, role, gene, source`), not
>   one sheet per quadrant.

**Files:**
- Modify: `04_Figures_v2/06_Cluster/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: `RRHO2` package; `load_combined_wide()`; `run_hallmark_ora()`, `build_*_panel()`, `build_cluster_row()`, `stack_cluster_rows()` from Task 2.
- Produces: `MAIN_F06_pilot_rrho2.{pdf,png}`; workbook sheets `pilot_rrho2_genelists` (long), `pilot_rrho2_ora`.

- [ ] **Step 1: Inspect RRHO2 exposed slots (read-only)**

```bash
Rscript -e 'suppressPackageStartupMessages(library(RRHO2)); cat(paste(ls("package:RRHO2"), collapse="\n"), "\n")'
```
Expected output includes `RRHO2_initialize` and `RRHO2_heatmap`; gene lists are accessible via `obj$genelist_uu`, `obj$genelist_dd`, `obj$genelist_ud`, `obj$genelist_du` on the returned object. (If your installed RRHO2 names them differently, adapt the field names below — keep the four quadrant-set extractions consistent.)

- [ ] **Step 2: Extend verification script**

Save as `/tmp/verify_F06_rrho2.R`:

```r
source("/tmp/verify_F06_logfc.R")
root <- "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026"
base <- file.path(root, "04_Figures_v2", "06_Cluster")
suppressPackageStartupMessages({ library(openxlsx) })
png <- file.path(base, "b_reports", "main", "png", "MAIN_F06_pilot_rrho2.png")
stopifnot("pilot_rrho2 png missing" = file.exists(png), file.info(png)$size > 10000)
sh <- getSheetNames(file.path(base, "c_data", "F06_supplementary.xlsx"))
need <- c("pilot_rrho2_genelists","pilot_rrho2_ora")
miss <- setdiff(need, sh)
if (length(miss)) stop("missing rrho2 sheets: ", paste(miss, collapse=", "))
cat("F06 RRHO2 pilot verification PASS\n")
```

- [ ] **Step 3: Run — confirm FAIL**

```bash
Rscript /tmp/verify_F06_rrho2.R
```
Expected: FAIL on the rrho2 png check.

- [ ] **Step 4: Append Pilot 6 to the F06 driver**

Append just before the `build_workbook(...)` call:

```r
# PILOT 6 — RRHO2 Disease <-> Rescue threshold-free concordance map
if (requireNamespace("RRHO2", quietly = TRUE)) {
  message("\n=== pilot_rrho2 ===")
  cw <- load_combined_wide()
  rank_vec <- function(p, lfc) {
    s <- -log10(p) * sign(lfc)
    s[!is.finite(s)] <- 0
    s
  }
  d_rank <- rank_vec(cw$P.Value_CTLvPHE,      cw$logFC_CTLvPHE)
  r_rank <- rank_vec(cw$P.Value_PHEvPHE_MITO, cw$logFC_PHEvPHE_MITO)
  keep <- !is.na(d_rank) & !is.na(r_rank) & !is.na(cw$gene) & nzchar(cw$gene)
  list1 <- data.frame(gene = cw$gene[keep], score = d_rank[keep])
  list2 <- data.frame(gene = cw$gene[keep], score = r_rank[keep])
  set.seed(SEED)
  rr <- RRHO2::RRHO2_initialize(list1, list2,
                                labels = c("Disease (CTLvPHE)", "Rescue (PHEvPHE_MITO)"),
                                log10.ind = TRUE,
                                stepsize = ceiling(sqrt(nrow(list1))),
                                boundary = 0.025)

  # Heatmap (top half of figure)
  hm_pdf <- file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2_heatmap.pdf")
  hm_png <- file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2_heatmap.png")
  pdf_dev(hm_pdf, width = 120 / 25.4, height = 120 / 25.4); RRHO2::RRHO2_heatmap(rr); dev.off()
  png(hm_png, width = 120, height = 120, units = "mm", res = 300); RRHO2::RRHO2_heatmap(rr); dev.off()

  # Quadrant gene lists
  quad_lists <- list(
    UU = rr$genelist_uu,
    DD = rr$genelist_dd,
    UD = rr$genelist_ud,    # Disease up, Rescue down  -> reversed
    DU = rr$genelist_du)    # Disease down, Rescue up  -> reversed
  quad_lists <- lapply(quad_lists, function(x)
    if (is.null(x)) character(0) else as.character(x))
  pal_q <- c(UU = "#2E7D32", DD = "#1565C0", UD = "#B2182B", DU = "#D6604D")
  quad_role <- c(UU = "Concordant Up", DD = "Concordant Down",
                 UD = "Reversed (Disease Up / Rescue Down)",
                 DU = "Reversed (Disease Down / Rescue Up)")
  rows <- lapply(names(quad_lists), function(q) {
    g_in <- quad_lists[[q]]
    if (length(g_in) < 5) return(NULL)
    g_in <- intersect(g_in, rownames(group_mat))
    if (length(g_in) < 5) return(NULL)
    z_cl <- standardise_genes(group_mat[g_in, , drop = FALSE])
    color <- pal_q[q]
    ora <- run_hallmark_ora(g_in, universe = ALL_GENES)
    hdr <- sprintf("Quadrant %s  |  n = %d  |  %s", q, length(g_in), quad_role[q])
    build_cluster_row(
      traj_plot = build_trajectory_panel(z_cl, cluster = q,
                                         x_levels = H9C2_GROUP_LEVELS,
                                         x_lab = "condition (group mean)",
                                         color = color, kind = "line"),
      ora_plot  = build_ora_bar_panel(ora, color = color, max_n = 6),
      header_text = hdr, color = color)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]

  fig <- stack_cluster_rows(rows,
    title    = "F06 pilot_rrho2 — Disease<->Rescue quadrants (per-quadrant ORA)",
    subtitle = "heatmap saved separately; rows = significant RRHO2 quadrants")
  h_mm <- 32 + 32 * length(rows)
  ggsave(file.path(MAIN_PDF, "MAIN_F06_pilot_rrho2.pdf"), fig,
         width = FIG_W, height = h_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(MAIN_PNG, "MAIN_F06_pilot_rrho2.png"), fig,
         width = FIG_W, height = h_mm, units = "mm", dpi = 300, limitsize = FALSE)

  # Workbook payloads
  genelist_long <- bind_rows(lapply(names(quad_lists), function(q)
    tibble(quadrant = q, role = quad_role[q], gene = quad_lists[[q]])))
  ora_rr <- bind_rows(lapply(names(quad_lists), function(q) {
    o <- run_hallmark_ora(quad_lists[[q]], universe = ALL_GENES)
    if (is.null(o) || nrow(o) == 0) return(NULL)
    mutate(o, quadrant = q, role = quad_role[q])
  }))

  results$pilot_rrho2 <- list(
    key = "pilot_rrho2",
    sheets = list(
      pilot_rrho2_genelists = genelist_long,
      pilot_rrho2_ora       = if (is.null(ora_rr)) tibble() else ora_rr))
} else {
  warning("RRHO2 package missing — skipping pilot_rrho2")
}
```

And **extend the sheet-spec assembly** for `pilot_rrho2`:

```r
if (!is.null(results$pilot_rrho2)) {
  rr <- results$pilot_rrho2$sheets
  sheet_specs <- c(sheet_specs,
    list(list(name = "pilot_rrho2_genelists", df = rr$pilot_rrho2_genelists,
              role = "Per-quadrant gene lists from RRHO2 (UU/DD/UD/DU)",
              contents = "quadrant (UU=concordant up, DD=concordant down, UD/DU=reversed), role, gene")),
    list(list(name = "pilot_rrho2_ora", df = rr$pilot_rrho2_ora,
              role = "Hallmark ORA per RRHO2 quadrant",
              contents = "quadrant, role, pathway, padj, overlap, size, odds_ratio")))
  sheet_specs[[1]]$df <- bind_rows(
    sheet_specs[[1]]$df,
    tibble(Pilot = "pilot_rrho2",
           Method = "RRHO2 threshold-free Disease<->Rescue map",
           Gate = "all genes with non-NA P.Value+logFC in Disease and Rescue",
           N_genes = nrow(rr$pilot_rrho2_genelists),
           Fuzzifier_m = NA_real_, Cluster_c = 4L))
}
```

- [ ] **Step 5: Run + verify**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R
Rscript /tmp/verify_F06_rrho2.R
```
Expected: `F06 RRHO2 pilot verification PASS`.

- [ ] **Step 6: Visual check**

Read `04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_pilot_rrho2.png` and `MAIN_F06_pilot_rrho2_heatmap.png`. The heatmap shows the rank-rank concordance hotspots; the row figure shows up to 4 quadrant rows each with trajectory + Hallmark ORA, headers labeled `Quadrant UD  |  n = …  |  Reversed (Disease Up / Rescue Down)`.

- [ ] **Step 7: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git add 04_Figures_v2/06_Cluster/
git commit -m "add the RRHO2 pilot for F06 with quadrant trajectories and per-quadrant ORA"
```

---

## Task 7: End-to-end verification + drop the obsolete `shown_pathways.csv`

> **[UPDATED — see Post-implementation status]** Current verifier set lives at
> `/tmp/verify_F04.R` and `/tmp/verify_F06_{cmeans,wgcna,logfc,rrho2}.R`. The
> RRHO2 verifier transitively sources the other three F06 verifiers, so a
> single `Rscript /tmp/verify_F06_rrho2.R` exercises all four sub-pilots in
> one call. Last confirmed PASS state: commit `0157237` (2026-06-20).

**Files:**
- (None modified; runtime verification only.)

**Interfaces:**
- Consumes: every output produced by Tasks 0–6.

- [ ] **Step 1: Full re-run from a clean shell**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R
```
Expected: both scripts exit 0. F06 messages enumerate all 6 pilots in order.

- [ ] **Step 2: Run the full verification chain**

```bash
Rscript /tmp/verify_F04.R
Rscript /tmp/verify_F06_rrho2.R   # transitively runs the others via source()
```
Expected: both report PASS.

- [ ] **Step 3: Confirm no stale orphans**

```bash
ls 04_Figures_v2/06_Cluster/b_reports/main/png/
```
Expected: exactly six `MAIN_F06_pilot_*.png` (plus the RRHO2 heatmap), no `MAIN_F06_group_all*.png`, `MAIN_F06_sample_all*.png`, or `MAIN_F06_group_pi*.png` from the old design.

If old orphans remain, delete them:
```bash
rm -f 04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_group_all*.png \
      04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_sample_all*.png \
      04_Figures_v2/06_Cluster/b_reports/main/png/MAIN_F06_group_pi*.png \
      04_Figures_v2/06_Cluster/b_reports/main/pdf/MAIN_F06_group_all*.pdf \
      04_Figures_v2/06_Cluster/b_reports/main/pdf/MAIN_F06_sample_all*.pdf \
      04_Figures_v2/06_Cluster/b_reports/main/pdf/MAIN_F06_group_pi*.pdf \
      04_Figures_v2/06_Cluster/b_reports/supp/png/MAIN_F06_group_all*.png \
      04_Figures_v2/06_Cluster/b_reports/supp/png/MAIN_F06_sample_all*.png \
      04_Figures_v2/06_Cluster/b_reports/supp/png/MAIN_F06_group_pi*.png \
      04_Figures_v2/06_Cluster/b_reports/supp/pdf/MAIN_F06_group_all*.pdf \
      04_Figures_v2/06_Cluster/b_reports/supp/pdf/MAIN_F06_sample_all*.pdf \
      04_Figures_v2/06_Cluster/b_reports/supp/pdf/MAIN_F06_group_pi*.pdf
```

Verify `04_Pathway_bars/c_data/shown_pathways.csv` is gone (the rewritten F04 deletes it on run):
```bash
ls 04_Figures_v2/04_Pathway_bars/c_data/
```
Expected: no `shown_pathways.csv`.

- [ ] **Step 4: Final commit (if any orphans were swept)**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git status
# if anything is staged after the orphan sweep:
git add -u
git commit -m "sweep the orphaned F06 outputs from the prior design"
```

- [ ] **Step 5: Print a short summary**

```bash
echo "F04: $(ls 04_Figures_v2/04_Pathway_bars/b_reports/main/png/)"
echo "F06: $(ls 04_Figures_v2/06_Cluster/b_reports/main/png/)"
```
Expected: F04 = `MAIN_F04_pathway_bars.png`; F06 = the 6 pilot PNGs + RRHO2 heatmap PNG.

---

## Self-review checklist

- Spec coverage: every section of `docs/superpowers/specs/2026-06-20-figures-v2-pathway-cluster-redesign.md` maps to a task above (doc-rot → Task 0; F04 → Task 1; helpers → Task 2; Pilots 1–3 → Task 3; Pilot 4 → Task 4; Pilot 5 → Task 5; Pilot 6 → Task 6; end-to-end + cleanup → Task 7).
- Placeholder scan: no TBDs; every code step shows actual code; every shell step shows the actual command and expected output.
- Type consistency: helpers introduced in Task 2 are called by name in Tasks 3–6 with the signatures defined here (`filter_sig_in_any_contrast`, `cluster_palette`, `build_trajectory_panel`, `build_ora_bar_panel`, `build_cluster_row`, `stack_cluster_rows`, `run_hallmark_ora`, `load_wgcna_modules`, `compute_me_contrast_correlations`, `classify_module_sign_pattern`). The pilot driver appends to a single `results` list and a single `sheet_specs` accumulator throughout.
- Scope: one figure suite, two scripts; no architecture changes upstream.
