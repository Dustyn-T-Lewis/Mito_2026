# 04_Figures_v2 coherence + composites — design

Date: 2026-06-22
Scope: figures C1 and C2 (PCA, DEP bars, Venn, pathway bars, enrichment rings).
C3 (clustering) is deferred to its own later spec.

## Goal

Make the six-figure Mito (H9c2) suite read as one coherent set, then stitch the
panels into two composites built in R. Refine each standalone panel until it is
clean on its own, verify it renders, then assemble. Validate methods and trim
code to human standards along the way (root `.claude/` rules).

## Decisions (locked)

- **Composites are built in R** (patchwork), not hand-composited. New
  `07_Composites/` holds the assembler scripts.
- **Role names everywhere**: Disease, Transplant, Rescue, Interaction, Secondary.
  Contrast algebra (`PHE - Ctl`, ...) is stated once per composite caption.
- **C1+C2 now; C3 later.** C3 opens with a deep literature sweep on clustering
  methods (over-large-cluster concern), the GO BP/CC/MF + `simplify()` enrichment
  choice, and code validation/optimization. Its own brainstorm + spec.
- **Direction owns red/blue**: `Up = #D6604D`, `Down = #4393C3`, `NS = grey70`,
  used by every bar/volcano/ring/strip; mito subset = same hue, darker shade.
- **Groups recolored** (PCA legend only): `Ctl = #3B7DB5`, `Mito = #009E73`,
  `PHE = #E08214`, `PHE_Mito = #8073AC` — no group equals direction red/blue.
- **C2 rings**: Disease + Transplant + Rescue in the main composite; Interaction
  and Secondary go supplementary.
- **Panel-background tints dropped** — clean white panels across all figures
  (the tints keyed on contrast, which cannot apply to the group-keyed PCA, so
  they could never be universal).

## A. Consistency contract — `functions/01_style_palettes_theme.R`

| Element | Rule |
|---|---|
| Naming | `role_label()` maps contrast code -> role; used by every panel and axis. |
| Direction palette | `DIR_COLORS` (Up/Down/NS) owns red/blue; mito subset darker shade. All inline hex (F03 strip, F04 fills) replaced by the constant. |
| Group palette | recolored `GROUP_COLORS`; PCA only. |
| Contrast order | narrative: Disease -> Transplant -> Rescue -> Interaction. |
| Panel tags | single bold uppercase `A,B,C...` per composite, top-left, `BASE_TAG`. |
| Title / subtitle | title bold; subtitle italic grey, one short informative line; method detail moves to caption/workbook. |
| Backgrounds | white; contrast tints removed. |
| Type scale | existing `FIG_*` + `scale_text()` kept, verified identical across panels. |
| Dead code | remove unused YvO-ported helpers (`classify_proteins_f2/f3`, `SIG_*_F2/F3`, `ORA_QUAD_COLORS_*`, `make_sigmoid_ribbon`, `composite_text_sizes`, `strip_for_composite`). Verified 0 consumers. |

## B. Build architecture (build-once)

Each panel is built by a `build_*()` function returning the ggplot object(s).
Both the standalone driver and the composite assembler consume that function, so
standalone PDFs/PNGs regenerate unchanged and composites are pure assembly.

```
functions/08_composite_layout.R       NEW shared stitcher: tag theme, caption
                                       block, NES-legend placement, save_composite()
<fig>/a_script/_build.R                NEW build_*() returns the ggplot(s)
<fig>/a_script/01_main_panels.R        thin driver: source builder, save standalone
07_Composites/a_script/01_C1_overview.R   source 3 builders, patchwork-stitch
07_Composites/a_script/02_C2_enrichment.R source F04+F05 builders, stitch
07_Composites/b_reports/main/{pdf,png}/
```

`functions/` keeps pipeline-order prefixes; `08` (composite) is the last step.

## C. Per-panel refinement (do first, one at a time, verify each renders)

| Panel | Visual | Method/code |
|---|---|---|
| F01 PCA | stats top-left, borderless + compact (no point overlap); tighten axes to cut whitespace (rein in the lone PHE outlier's pull); group legend stacked in the bottom-right corner; informative subtitle; recolored groups; clean tag. | note pairwise PERMANOVA = 3 uncorrected tests in caption; keep betadisper guard. |
| F02 DEP bars | drop tint; direction palette; role names; narrative order; clean subtitle. | replace hand-rolled `dir_off/thr_off` dodge with `position_dodge`; fix `abs(logFC) <= 1` truncation that biases `med|LFC|` (compute on full data; label the plotted view as windowed). |
| F03 Venn | inline hex -> `DIR_COLORS`; role names; clean subtitle/tag. | report eulerr `stress`/`diagError` to the workbook. |
| F04 Pathway bars | drop tint; fills aligned to `DIR_COLORS` (mito = darker shade); de-cram subtitle. | confirm dedup params; no structural change. |
| F05 rings | per-contrast tint -> consistent neutral ring bg; role names in titles; tags. | keep fgsea + dedup pipeline; no structural change. |

## D. Composite layout (`07_Composites/`, after panels are clean)

- C1 (who differs): `A` PCA (large, top-left) - `B` DEP counts - `C` effect size
  - `D` Venn - `E` direction strip. One shared direction legend.
- C2 (what biology): `A` pathway-count bars (top) - `B` Disease ring -
  `C` Transplant ring - `D` Rescue ring; NES legend placed once.

## E. Verification

Per panel and per composite: script runs clean; PNG/PDF render; supplementary
`.xlsx` still writes; `lintr::lint_dir` / `styler` clean on changed files;
AI-tell grep on changed files. Standalone outputs must not regress.

## Out of scope

- C3 clustering (own spec).
- Any change to upstream stages 00-03 or to `04_Figures/`.
- New gene-set derivation.
