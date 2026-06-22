# 04_Figures_v2 Coherence + Composites Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the five C1/C2 panels of the Mito (H9c2) figure suite visually coherent, fix their method/code issues, then stitch them into two R-built composites.

**Architecture:** Centralize palettes/naming in `functions/01`; repoint the existing `CONTRAST_DISPLAY_MAP` to role names so every `contrast_brief()` caller becomes role-named for free; refactor each panel script so a `build_*()` function returns the ggplot object(s) that both the standalone driver and a new `07_Composites/` assembler consume (build-once). C3 clustering is out of scope.

**Tech Stack:** R, ggplot2, patchwork, eulerr, vegan, fgsea cache (precomputed), `here`, `openxlsx` (via `build_workbook`).

## Global Constraints

- Role names everywhere: `Disease`, `Transplant`, `Rescue`, `Interaction`, `Secondary`. Algebra stated once per composite caption.
- Direction palette owns red/blue: `Up = #D6604D`, `Down = #4393C3`, `NS = grey70`; mito subset = same hue, darker shade.
- Group palette (PCA only, recolored): `Ctl = #3B7DB5`, `Mito = #009E73`, `PHE = #E08214`, `PHE_Mito = #8073AC`.
- Narrative contrast order everywhere: Disease → Transplant → Rescue → Interaction.
- C2 main rings: Disease + Transplant + Rescue. Interaction + Secondary supplementary.
- White panel backgrounds; no contrast tints.
- No edits to stages 00–03 or to `04_Figures/`. No new gene-set derivation. Read-only inputs.
- Width invariant 178 mm (`PANEL_MD`); pdf via `get_pdf_device()`, png 300 dpi, `units = "mm"`.
- Commit style: one short lowercase sentence, no AI/attribution trailers.
- Verification loop for figure scripts = run via `Rscript`, assert outputs exist + dims, eyeball the PNG. No testthat harness exists for figures; pure helpers get inline `stopifnot` checks.

---

### Task 1: Recolor + rename in the style engine, trim dead code

**Files:**
- Modify: `04_Figures_v2/functions/01_style_palettes_theme.R`
- Modify: `04_Figures_v2/functions/02_data_paths_and_loaders.R:14-20` (repoint map)

**Interfaces:**
- Produces: `DIR_COLORS` (named `Up`/`Down`/`NS`), `DIR_COLORS_MITO` (named `Up`/`Down`, darker), `GROUP_COLORS` (recolored), `role_label(old_name)` returning role strings, `H9C2_CONTRAST_ORDER` (`c("CTLvPHE","CTLvMITO","PHEvPHE_MITO","Interaction")`). `contrast_brief()` now returns role names.

- [ ] **Step 1: Repoint the display map to role names** in `02_data_paths_and_loaders.R`, replacing the `CONTRAST_DISPLAY_MAP` block (lines 14-20):

```r
CONTRAST_DISPLAY_MAP <- c(
  CTLvMITO      = "Transplant",   # Mito - Ctl
  CTLvPHE       = "Disease",      # PHE - Ctl
  PHEvPHE_MITO  = "Rescue",       # PHE_Mito - PHE
  Interaction   = "Interaction",  # orthogonal 2x2
  MITOvPHE_MITO = "Secondary"     # PHE_Mito - Mito
)
```

- [ ] **Step 2: Add `role_label` alias** at the end of the naming block in `02_data_paths_and_loaders.R` (after `contrast_brief`, ~line 36):

```r
# Role name for a contrast (Disease / Transplant / Rescue / Interaction / Secondary).
role_label <- contrast_brief
```

- [ ] **Step 3: Recolor groups + add direction constants** in `01_style_palettes_theme.R`. Replace the `H9C2_PAL_GROUP` and `H9C2_PAL_DIR` lines (15-16) with:

```r
H9C2_PAL_GROUP <- c(Ctl = "#3B7DB5", Mito = "#009E73", PHE = "#E08214", PHE_Mito = "#8073AC")
H9C2_PAL_DIR   <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")
H9C2_PAL_DIR_MITO <- c(Up = "#B2182B", Down = "#2166AC")   # darker mito subset
H9C2_CONTRAST_ORDER <- c("CTLvPHE", "CTLvMITO", "PHEvPHE_MITO", "Interaction")
```

- [ ] **Step 4: Export the direction aliases** near the existing alias block (after line 27 `DIR_COLORS <- H9C2_PAL_DIR`):

```r
DIR_COLORS_MITO <- H9C2_PAL_DIR_MITO
```

- [ ] **Step 5: Delete the dead YvO-ported block** — remove lines 95-145 of `01_style_palettes_theme.R` (`SIG_COLORS_F2` through `ORA_QUAD_COLORS_F3`) and lines 147-173 (`classify_proteins_f2`, `classify_proteins_f3`), plus `composite_text_sizes` (68-74), `strip_for_composite` (50-53), and `make_sigmoid_ribbon` (228-237). Leave `is_light_color` (it has 1 consumer).

- [ ] **Step 6: Verify the engine loads and constants are right**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript -e 'source(here::here("04_Figures_v2","functions","02_data_paths_and_loaders.R"));
  stopifnot(DIR_COLORS[["Up"]]=="#D6604D", GROUP_COLORS[["PHE"]]=="#E08214",
            role_label("CTLvPHE")=="Disease", role_label("PHEvPHE_MITO")=="Rescue",
            !exists("classify_proteins_f2"), !exists("make_sigmoid_ribbon"));
  cat("OK\n")'
```
Expected: `OK`

- [ ] **Step 7: Confirm no live consumer of the deleted symbols**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026/04_Figures_v2
grep -rEl "classify_proteins_f[23]|SIG_COLORS_F[23]|ORA_QUAD_COLORS_F[23]|make_sigmoid_ribbon|composite_text_sizes|strip_for_composite" --include="*.R" . | grep -v 01_style_palettes_theme.R
```
Expected: no output.

- [ ] **Step 8: Commit**

```bash
git add 04_Figures_v2/functions/01_style_palettes_theme.R 04_Figures_v2/functions/02_data_paths_and_loaders.R
git commit -m "recolor groups, give direction its own red/blue, trim the dead yvo helpers"
```

---

### Task 2: Composite layout helper

**Files:**
- Create: `04_Figures_v2/functions/08_composite_layout.R`

**Interfaces:**
- Consumes: `get_pdf_device()`, `FIG_THEME`, `PANEL_MD` from `01`.
- Produces: `add_tag(p, tag)` (adds a bold top-left panel tag), `save_composite(plot, base, name, width_mm, height_mm)` (writes pdf+png to `<base>/b_reports/main/{pdf,png}/`), `composite_caption(text)` (a `plot_annotation` caption with the algebra line).

- [ ] **Step 1: Write the helper**

```r
# 04_Figures_v2/functions/08_composite_layout.R
# Pipeline step 08: stitch refined standalone panels into composites. The last
# assembly step — sources nothing the panels do not already provide.
source(here::here("04_Figures_v2", "functions", "01_style_palettes_theme.R"))

library(patchwork)

add_tag <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag = element_text(face = "bold", size = BASE_TAG),
        plot.tag.position = c(0.01, 0.99))

composite_caption <- function(text)
  patchwork::plot_annotation(
    caption = text,
    theme = theme(plot.caption = element_text(size = 5, color = "grey35",
                                              hjust = 0, lineheight = 1.1)))

save_composite <- function(plot, base, name, width_mm, height_mm) {
  pdf_dev <- get_pdf_device()
  out_pdf <- file.path(base, "b_reports", "main", "pdf")
  out_png <- file.path(base, "b_reports", "main", "png")
  for (d in c(out_pdf, out_png)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_pdf, paste0(name, ".pdf")), plot, width = width_mm,
         height = height_mm, units = "mm", device = pdf_dev, limitsize = FALSE)
  ggsave(file.path(out_png, paste0(name, ".png")), plot, width = width_mm,
         height = height_mm, units = "mm", dpi = 300, limitsize = FALSE)
  invisible(file.path(out_png, paste0(name, ".png")))
}
```

- [ ] **Step 2: Verify it sources and tags a toy plot**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript -e 'library(ggplot2); source(here::here("04_Figures_v2","functions","08_composite_layout.R"));
  p <- add_tag(ggplot(mtcars,aes(mpg,wt))+geom_point(), "A");
  stopifnot(inherits(p,"ggplot")); cat("OK\n")'
```
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add 04_Figures_v2/functions/08_composite_layout.R
git commit -m "add the composite layout helper as functions step 08"
```

---

### Task 3: Refactor F01 PCA into a builder and refine it

**Files:**
- Create: `04_Figures_v2/01_PCA/a_script/_build.R`
- Modify: `04_Figures_v2/01_PCA/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: Task 1 constants, `P05`, `H9C2_GROUP_LEVELS`, `fmt_p`.
- Produces: `build_pca_panel()` returning a list `list(plot=<ggplot>, scores=<df>, permanova=<df>)`.

- [ ] **Step 1: Create `_build.R`** holding everything from the current `01_main_panels.R` lines 20-91 wrapped in a function, with the refinements applied. Key changes vs current code:
  - `GRP_LAB` becomes `c(Ctl="Ctl", Mito="Mito", PHE="PHE", PHE_Mito="PHE+Mito")` → use the role-free group names but fix the legend; keep `PHE_Mito` label as `"PHE+Mito"` is fine (group, not contrast).
  - Stats annotation: drop the bordered `annotate("label", ...)`; render top-left with `label.size = 0`, `fill = NA`, smaller `size = 1.6`, anchored to the panel via `annotation_custom`/`patchwork::inset_element` of a `geom_text` strip so it never overlaps points. Simplest robust approach: place it in the top-left using `coord_cartesian` headroom but with `fill = NA, label.size = 0` and left-justified.
  - Tighten whitespace: replace the `ytop <- yr[2] + 0.42*diff(yr)` headroom hack with `scale_x/y_continuous(expand = expansion(mult = 0.04))` and put the stats in a slim top strip via `patchwork` instead of stealing 42% of the panel.
  - Group legend: move to bottom-right, vertical stack — `legend.position = c(0.99, 0.02)`, `legend.justification = c(1, 0)`, `guide_legend(ncol = 1)`.
  - Subtitle: `sprintf("PERMANOVA Group R² = %.2f, %s | n = %d, %s proteins (imputed)", perm_R2, fmt_p(perm_p), n, fmt_proteins)`.

```r
# 04_Figures_v2/01_PCA/a_script/_build.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(vegan); library(patchwork)
})
source(here::here("04_Figures_v2", "functions", "02_data_paths_and_loaders.R"))

build_pca_panel <- function() {
  dal_imp  <- readRDS(P05$imp_rds)
  imp_mat  <- as.matrix(dal_imp$data)
  imp_meta <- as_tibble(dal_imp$metadata)
  imp_meta$Group <- factor(imp_meta$Group, levels = H9C2_GROUP_LEVELS)
  imp_mat  <- imp_mat[, imp_meta$Col_ID]

  pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
  var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
  pca_df <- as.data.frame(pca$x[, 1:2]) |>
    mutate(Col_ID = rownames(pca$x)) |>
    left_join(dplyr::select(imp_meta, Col_ID, Group), by = "Col_ID") |>
    mutate(Group = factor(Group, levels = H9C2_GROUP_LEVELS))

  set.seed(42); dist_mat <- dist(scale(t(imp_mat)))
  set.seed(42)
  perm <- adonis2(dist_mat ~ Group, data = imp_meta, permutations = 9999, by = "terms")
  perm_R2 <- perm["Group", "R2"]; perm_p <- perm["Group", "Pr(>F)"]

  pairs_to_test <- list(`Disease`    = c("Ctl", "PHE"),
                        `Transplant` = c("Ctl", "Mito"),
                        `Rescue`     = c("PHE", "PHE_Mito"))
  pair_res <- bind_rows(lapply(names(pairs_to_test), function(nm) {
    pr <- pairs_to_test[[nm]]; keep <- imp_meta$Group %in% pr
    sm <- imp_mat[, keep]; smeta <- imp_meta[keep, ]
    smeta$Group <- droplevels(smeta$Group)
    set.seed(42)
    r <- adonis2(dist(scale(t(sm))) ~ Group, data = smeta,
                 permutations = 9999, by = "terms")
    tibble(role = nm, R2 = r$R2[1], p = r$`Pr(>F)`[1])
  }))

  set.seed(42)
  bd_p <- permutest(betadisper(dist_mat, imp_meta$Group),
                    permutations = 999)$tab$`Pr(>F)`[1]
  if (!is.na(bd_p) && bd_p < 0.05)
    warning("Heterogeneous group dispersions (betadisper p < 0.05)")

  permanova_out <- bind_rows(
    tibble(role = "Group (overall)", R2 = perm_R2, p = perm_p),
    pair_res,
    tibble(role = "dispersion (betadisper)", R2 = NA_real_, p = bd_p))

  GRP_LAB <- c(Ctl = "Ctl", Mito = "Mito", PHE = "PHE", PHE_Mito = "PHE+Mito")
  GRP_SHP <- c(Ctl = 16, Mito = 17, PHE = 15, PHE_Mito = 18)
  fmt_perm <- function(role, r2, p) sprintf("%s R²=%.2f, %s", role, r2, fmt_p(p))
  stat_lines <- paste(c(
    fmt_perm("Disease", pair_res$R2[1], pair_res$p[1]),
    fmt_perm("Transplant", pair_res$R2[2], pair_res$p[2]),
    fmt_perm("Rescue", pair_res$R2[3], pair_res$p[3])), collapse = "\n")

  n  <- nrow(imp_meta)
  np <- format(nrow(imp_mat), big.mark = ",")
  p <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = Group)) +
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.08,
                 level = 0.80, show.legend = FALSE) +
    stat_ellipse(level = 0.80, linewidth = 0.3, linetype = "dashed",
                 show.legend = FALSE) +
    geom_point(size = 1.8, alpha = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = stat_lines, hjust = -0.05,
             vjust = 1.2, size = 1.6, color = "grey25", lineheight = 0.95) +
    scale_color_manual(values = GROUP_COLORS, labels = GRP_LAB, name = NULL,
                       guide = guide_legend(ncol = 1, override.aes = list(size = 2))) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_shape_manual(values = GRP_SHP, labels = GRP_LAB, name = NULL,
                       guide = guide_legend(ncol = 1)) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +
    scale_y_continuous(expand = expansion(mult = 0.04)) +
    labs(title = "Sample PCA",
         subtitle = sprintf("PERMANOVA Group R²=%.2f, %s  |  n=%d, %s proteins (imputed)",
                            perm_R2, fmt_p(perm_p), n, np),
         x = sprintf("PC1 (%.1f%%)", var_pct[1]),
         y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
    FIG_THEME +
    theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key = element_blank(), legend.key.size = unit(2.8, "mm"),
          legend.spacing.y = unit(0.3, "mm"),
          plot.margin = margin(4, 3, 2, 2))

  list(plot = p, scores = pca_df, permanova = permanova_out)
}
```

- [ ] **Step 2: Slim the driver** — replace `01_main_panels.R` body (keep the header comment) so it calls the builder and keeps the existing workbook + ggsave:

```r
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "functions", "06_supplementary_workbook.R"))

BASE <- here::here("04_Figures_v2", "01_PCA")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_dev <- get_pdf_device()

res <- build_pca_panel()
ggsave(file.path(RPT_PDF, "MAIN_F01_pca.pdf"), res$plot, width = 110, height = 95,
       units = "mm", device = pdf_dev)
ggsave(file.path(RPT_PNG, "MAIN_F01_pca.png"), res$plot, width = 110, height = 95,
       units = "mm", dpi = 300)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  figure_title = "F01 — Sample PCA + PERMANOVA on the imputed protein matrix",
  sheet_specs = list(
    list(name = "pca_scores", df = res$scores,
         role = "Panel coordinates — the PCA scatter points",
         contents = "PC1/PC2 scores per sample (Col_ID) with Group; % variance in axis titles"),
    list(name = "permanova", df = res$permanova,
         role = "Stats annotation block on the figure",
         contents = "adonis2 R2 and p for overall Group + Disease/Transplant/Rescue pairwise (3 uncorrected pairwise tests), plus betadisper dispersion p")))
message("F01 PCA rebuilt")
```

- [ ] **Step 3: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/01_PCA/a_script/01_main_panels.R && \
  test -f 04_Figures_v2/01_PCA/b_reports/main/png/MAIN_F01_pca.png && echo RENDERED
```
Expected: `RENDERED`. Then open the PNG and confirm: stats top-left without overlapping points, legend stacked bottom-right, tighter framing, recolored groups (PHE orange, Ctl blue, Mito green, PHE_Mito purple).

- [ ] **Step 4: Commit**

```bash
git add 04_Figures_v2/01_PCA/
git commit -m "refine the pca panel and split it into a reusable builder"
```

---

### Task 4: Refactor F02 DEP bars, fix the dodge and the logFC truncation

**Files:**
- Create: `04_Figures_v2/02_DEP_bars/a_script/_build.R`
- Modify: `04_Figures_v2/02_DEP_bars/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: Task 1 constants, `load_combined_wide`, `contrast_brief`/`role_label`, `DIR_COLORS`.
- Produces: `build_dep_count_panel()` and `build_dep_effect_panel()` each returning a ggplot; `dep_count_data()` returning the counts tibble for the workbook.

- [ ] **Step 1: Create `_build.R`** porting the current panel code with three fixes:
  1. **Narrative order:** `CORE <- H9C2_CONTRAST_ORDER` (Disease first), not the current Transplant-first vector.
  2. **Dodge:** delete the manual `dir_off/thr_off/xc` math (current lines 81-88, 94-100) and draw with `position_dodge2`. Use `aes(x = contrast, group = interaction(direction, threshold))` and `geom_col(position = position_dodge2(width = 0.9, padding = 0.1))`, mapping `alpha = threshold` (0.3/0.6/1.0) and `fill = direction` from `DIR_COLORS`, `linetype = direction` for the Up/Down outline.
  3. **Effect-size truncation:** keep the plotted window `abs(logFC) <= 1` for the histogram, but compute `med_abs` on the **full** `logFC` (remove the `<= 1` filter from the `summarise`). Relabel the stat `sprintf("median |log2FC| %.2f (all)", med_abs)` and the subtitle `"signed log2FC (±1 shown); median over all proteins"`.
  Drop the `geom_rect` contrast-tint backgrounds (current lines 95-97 and 150-151). Replace contrast fills in the effect panel with a neutral grey histogram (`fill = "grey75"`) so red/blue stays direction-only.

```r
# 04_Figures_v2/02_DEP_bars/a_script/_build.R   (key fragments)
CORE <- H9C2_CONTRAST_ORDER
ctr_levels <- role_label(CORE)
# counts panel:
p_counts <- ggplot(counts_df, aes(contrast, pct, fill = direction, alpha = threshold)) +
  geom_col(position = position_dodge2(width = 0.9, padding = 0.1),
           color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = DIR_COLORS[c("Up","Down")], name = NULL) +
  scale_alpha_manual(values = c(0.30, 0.60, 1.00), name = NULL) +
  ...
# effect panel stat (median over ALL proteins, not windowed):
lfc_stats <- bind_rows(lapply(CORE, \(c)
  tibble(contrast = role_label(c),
         med_abs = median(abs(dep_results[[c]]$logFC), na.rm = TRUE)))) |>
  mutate(lab = sprintf("median |log2FC| %.2f (all)", med_abs))
```

(Port the remaining structure — facet, density, workbook prep — verbatim from the current script, swapping `CTR_LAB`/`contrast_brief` results, which now return role names, and removing tint `geom_rect`s.)

- [ ] **Step 2: Slim the driver** to source `_build.R`, call both builders, `patchwork` them `p_counts + p_eff + plot_layout(widths = c(1, 0.55))`, ggsave to the existing F02 paths (178×100 mm), and write the workbook from `dep_count_data()` + per-contrast tables (unchanged logic).

- [ ] **Step 3: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/02_DEP_bars/a_script/01_main_panels.R && \
  test -f 04_Figures_v2/02_DEP_bars/b_reports/main/png/MAIN_F02_dep_bars.png && echo RENDERED
```
Expected: `RENDERED`. Confirm: Disease-first order, no panel tints, Up=red/Down=blue with faint→solid by threshold, median label says "(all)".

- [ ] **Step 4: Commit**

```bash
git add 04_Figures_v2/02_DEP_bars/
git commit -m "fix the dep-bar dodge and stop truncating the logfc median"
```

---

### Task 5: Refactor F03 Venn, use DIR_COLORS, report eulerr fit

**Files:**
- Create: `04_Figures_v2/03_Venn/a_script/_build.R`
- Modify: `04_Figures_v2/03_Venn/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: Task 1 constants, `contrast_brief`, `DIR_COLORS`.
- Produces: `build_venn_panels()` returning `list(venn=<ggplot>, strip=<ggplot>, membership=<df>, region_counts=<df>, fit_stats=<df>)`.

- [ ] **Step 1: Create `_build.R`** porting the current Venn logic with:
  - Strip fill `scale_fill_manual(values = DIR_COLORS[c("Up","Down")])` (replace the inline `c(Up="#D6604D", Down="#4393C3")` at current line 120).
  - Set colors already come from `CONTRAST_COLORS[SET_CONTRASTS]` — keep (these are the contrast hues for the three ellipses, distinct from direction).
  - After `eu_fit <- euler(...)`, capture fit quality: `fit_stats <- tibble(metric = c("stress","diagError"), value = c(eu_fit$stress, eu_fit$diagError))`.
  - Titles/subtitles use `contrast_brief()` (now role names).

- [ ] **Step 2: Slim the driver** to call the builder, `patchwork::wrap_plots(venn, strip, widths = c(1.6, 1))`, ggsave to F03 paths (140×100 mm), and add a third workbook sheet `euler_fit` from `fit_stats` (role: "Euler fit quality — stress and diagError; high values mean the area-proportional layout is approximate").

- [ ] **Step 3: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/03_Venn/a_script/01_main_panels.R && \
  test -f 04_Figures_v2/03_Venn/b_reports/main/png/MAIN_F03_venn.png && echo RENDERED
```
Expected: `RENDERED`. Confirm role names on the three sets, strip uses the canonical red/blue, workbook has the `euler_fit` sheet.

- [ ] **Step 4: Commit**

```bash
git add 04_Figures_v2/03_Venn/
git commit -m "wire the venn strip to the shared direction palette and log the euler fit"
```

---

### Task 6: Refactor F04 pathway bars, align fills, de-cram subtitle

**Files:**
- Create: `04_Figures_v2/04_Pathway_bars/a_script/_build.R`
- Modify: `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: Task 1 constants, `CANONICAL_DBS`, `MITO_PATHWAY_REGEX`, `deduplicate_enrichment`, `contrast_brief`.
- Produces: `build_pathway_bar_panel()` returning `list(plot=<ggplot>, bar_df=<df>, sig_pw=<df>)`.

- [ ] **Step 1: Create `_build.R`** porting current F04 with:
  - `CORE <- H9C2_CONTRAST_ORDER` (Disease first).
  - Fills from constants: `FILL_TOTAL <- c(Up = scales::alpha(DIR_COLORS[["Up"]], 0.55), Down = scales::alpha(DIR_COLORS[["Down"]], 0.55))` and `FILL_MITO <- DIR_COLORS_MITO` (deep crimson/navy) so the mito subset is the darker shade of the same hue.
  - Drop the `bg_df` contrast-tint `geom_rect` (current lines 76-77, 82-83).
  - Subtitle de-crammed: `"total bar = all DBs; dark = mitochondrial subset"`. Move the DB list + dedup detail into the workbook sheet `role`/`contents` text only.

- [ ] **Step 2: Slim the driver** to call the builder, keep the inline 4-key legend `inset_element`, ggsave to F04 paths (120×70 mm), write the workbook (unchanged sheet logic). Keep the `shown_pathways.csv` cleanup block.

- [ ] **Step 3: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R && \
  test -f 04_Figures_v2/04_Pathway_bars/b_reports/main/png/MAIN_F04_pathway_bars.png && echo RENDERED
```
Expected: `RENDERED`. Confirm: no tints, Up=red/Down=blue with darker mito overlay, short subtitle, Disease-first.

- [ ] **Step 4: Commit**

```bash
git add 04_Figures_v2/04_Pathway_bars/
git commit -m "align the pathway bars to the direction palette and drop the tint"
```

---

### Task 7: Refactor F05 rings, neutral background, role names

**Files:**
- Create: `04_Figures_v2/05_Enrich_Volcano/a_script/_build.R`
- Modify: `04_Figures_v2/05_Enrich_Volcano/a_script/01_main_panels.R`

**Interfaces:**
- Consumes: Task 1 constants, `make_volcano_ring`, `build_ring_180_split`, `build_nes_legend_bar`, `deduplicate_enrichment`, `CANONICAL_DBS`, `CONTRAST_MATH_BRIEF`.
- Produces: `build_ring(ctr, tag, role)` returning `list(plot=<ggplot>, tag=, ctr=, role=, terms=, full=)`; `build_nes_legend()` returning the legend ggplot.

- [ ] **Step 1: Create `_build.R`** porting current F05 `build_ring`, changing the ring background from the per-contrast tint to a neutral light grey: replace `bg_color = unname(CONTRAST_COLORS[ctr]), bg_alpha = 0.20` with `bg_color = "grey95", bg_alpha = 0.5`. Titles already use `contrast_brief()` (now role names). Keep filename tags (`ctlvphe` etc.) for backward-compatible file paths. Also define the legend wrapper the composite consumes:

```r
build_nes_legend <- function()
  build_nes_legend_bar(text_size = 9, title_size = 10,
                       bar_margin = margin(2, 6, 2, 6, "mm"))
```
(the driver in Step 2 calls `build_nes_legend()` instead of inlining `build_nes_legend_bar`.)

- [ ] **Step 2: Slim the driver** to source `_build.R`, `Map(build_ring, ...)` over the 4 contrasts, save each ring + the NES legend to F05 paths (unchanged), write the workbook + `shown_pathways.csv` (unchanged).

- [ ] **Step 3: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/05_Enrich_Volcano/a_script/01_main_panels.R && \
  ls 04_Figures_v2/05_Enrich_Volcano/b_reports/main/png/ | grep -c ring
```
Expected: `4`. Confirm rings have a neutral background and role-named titles.

- [ ] **Step 4: Commit**

```bash
git add 04_Figures_v2/05_Enrich_Volcano/
git commit -m "give the rings a neutral background and role-named titles"
```

---

### Task 8: C1 overview composite

**Files:**
- Create: `04_Figures_v2/07_Composites/a_script/01_C1_overview.R`

**Interfaces:**
- Consumes: `build_pca_panel`, `build_dep_count_panel`, `build_dep_effect_panel`, `build_venn_panels`, `add_tag`, `save_composite`, `composite_caption`.

- [ ] **Step 1: Write the assembler**

```r
#!/usr/bin/env Rscript
# C1 overview — PCA + DEP bars + Venn. Who differs.
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "01_PCA", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "02_DEP_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "03_Venn", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")

pca  <- build_pca_panel()$plot
dep  <- build_dep_count_panel()
eff  <- build_dep_effect_panel()
venn <- build_venn_panels()

design <- "
AABB
AACC
DDEE
"
fig <- add_tag(pca, "A") + add_tag(venn$venn, "B") +
  add_tag(dep, "C") + add_tag(eff, "D") + add_tag(venn$strip, "E") +
  plot_layout(design = design) +
  composite_caption(paste(
    "Disease = PHE − Ctl; Transplant = Mito − Ctl; Rescue = PHE_Mito − PHE.",
    "Significance Π < 0.05 (Xiao 2014). Up = red, Down = blue."))

save_composite(fig, BASE, "MAIN_C1_overview", width_mm = PANEL_MD, height_mm = 200)
message("C1 composite built")
```

(Adjust the `design` string after the first render so PCA reads as the dominant top-left panel; the grid above gives PCA a 2×2 block, Venn top-right, DEP counts + effect mid-right, strip bottom. Tune to taste on the eyeball pass.)

- [ ] **Step 2: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/07_Composites/a_script/01_C1_overview.R && \
  test -f 04_Figures_v2/07_Composites/b_reports/main/png/MAIN_C1_overview.png && echo RENDERED
```
Expected: `RENDERED`. Open it: five tagged panels A–E, one shared visual language, no clashing reds.

- [ ] **Step 3: Commit**

```bash
git add 04_Figures_v2/07_Composites/a_script/01_C1_overview.R 04_Figures_v2/07_Composites/b_reports/
git commit -m "stitch the c1 overview composite"
```

---

### Task 9: C2 enrichment composite

**Files:**
- Create: `04_Figures_v2/07_Composites/a_script/02_C2_enrichment.R`

**Interfaces:**
- Consumes: `build_pathway_bar_panel`, `build_ring`, `build_nes_legend`, `add_tag`, `save_composite`, `composite_caption`.

- [ ] **Step 1: Write the assembler** — pathway bars on top (tag A), the three rings Disease/Transplant/Rescue below (tags B/C/D), NES legend placed once.

```r
#!/usr/bin/env Rscript
# C2 enrichment — pathway-count bars + Disease/Transplant/Rescue rings. What biology.
source(here::here("04_Figures_v2", "functions", "08_composite_layout.R"))
source(here::here("04_Figures_v2", "04_Pathway_bars", "a_script", "_build.R"))
source(here::here("04_Figures_v2", "05_Enrich_Volcano", "a_script", "_build.R"))

BASE <- here::here("04_Figures_v2", "07_Composites")

bars   <- build_pathway_bar_panel()$plot
disease    <- build_ring("CTLvPHE",      "ctlvphe",      "Disease")$plot
transplant <- build_ring("CTLvMITO",     "ctlvmito",     "Transplant")$plot
rescue     <- build_ring("PHEvPHE_MITO", "phevphe_mito", "Rescue")$plot
nes    <- build_nes_legend()

design <- "
AAA
BCD
"
fig <- add_tag(bars, "A") + add_tag(disease, "B") +
  add_tag(transplant, "C") + add_tag(rescue, "D") +
  plot_layout(design = design, heights = c(0.8, 1)) +
  inset_element(nes, left = 0.35, right = 0.65, top = 0.52, bottom = 0.48) +
  composite_caption(paste(
    "5-DB lens (Hallmark/Reactome/KEGG/MitoCarta/GO Slim), EnrichmentMap dedup.",
    "Ring = top pathways by FDR; centre = protein volcano. NES colour bar shared."))

save_composite(fig, BASE, "MAIN_C2_enrichment", width_mm = PANEL_MD, height_mm = 175)
message("C2 composite built")
```

- [ ] **Step 2: Render and verify**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript 04_Figures_v2/07_Composites/a_script/02_C2_enrichment.R && \
  test -f 04_Figures_v2/07_Composites/b_reports/main/png/MAIN_C2_enrichment.png && echo RENDERED
```
Expected: `RENDERED`. Open it: bars on top, three role-named rings below, one NES legend, consistent fonts/tags with C1.

- [ ] **Step 3: Commit**

```bash
git add 04_Figures_v2/07_Composites/a_script/02_C2_enrichment.R 04_Figures_v2/07_Composites/b_reports/
git commit -m "stitch the c2 enrichment composite"
```

---

### Task 10: README, run-list, and the cleanliness sweep

**Files:**
- Modify: `04_Figures_v2/README.md`

- [ ] **Step 1: Update the README** — change the "No composite stitching" line (current line 6) to note `07_Composites/` builds C1 and C2 in R; add `07_Composites/` to the Layout tree and the `functions/` table (`08_composite_layout.R`); add the two composite scripts to the Run block; document the role-name + direction-palette conventions in the Conventions section.

- [ ] **Step 2: Lint and style the changed R**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
Rscript -e 'styler::style_dir("04_Figures_v2", filetype="R")'
Rscript -e 'lintr::lint_dir("04_Figures_v2")' 2>&1 | head -40
```
Expected: styler reports no/again-stable diff; lintr shows no new issues on the touched files.

- [ ] **Step 3: AI-tell grep on the changed files**

Run:
```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
git diff --name-only main...HEAD -- '*.R' | xargs grep -nE "Step [0-9]:|# ===|leverage|comprehensive|seamlessly|delve|robustly|temp_df" || echo "clean"
```
Expected: `clean` (or only legitimate matches you then remove).

- [ ] **Step 4: Full regeneration smoke test** — run all five panels + both composites end to end:

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_Mito_2026
for s in 01_PCA/a_script/01_main_panels 02_DEP_bars/a_script/01_main_panels \
         03_Venn/a_script/01_main_panels 04_Pathway_bars/a_script/01_main_panels \
         05_Enrich_Volcano/a_script/01_main_panels \
         07_Composites/a_script/01_C1_overview 07_Composites/a_script/02_C2_enrichment; do
  echo "== $s =="; Rscript "04_Figures_v2/$s.R" || exit 1
done; echo "ALL GREEN"
```
Expected: `ALL GREEN`.

- [ ] **Step 5: Commit**

```bash
git add 04_Figures_v2/README.md 04_Figures_v2/
git commit -m "document the composites and run the cleanliness sweep"
```

---

## Notes for the implementer

- Inputs are read-only; never touch `02_Normalization/`, `03_DEP/`, or `04_Figures/`.
- The supplementary `.xlsx` are timestamp-free and must keep writing for every figure.
- When porting a panel into `_build.R`, copy the existing tested logic verbatim except for the listed changes — the goal is build-once, not a rewrite.
- Composite `design`/`heights`/`inset` values are first-pass; tune them on the eyeball render. Width stays `PANEL_MD` (178 mm).
