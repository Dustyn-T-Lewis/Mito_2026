# 04_Figures — Mito (H9c2) figure suite

Standalone figure suite for the H9c2 mitochondrial-transplant proteomics project,
parallel to the YvO proteome-plasticity manuscript. Three self-contained figures
(F01 proteome overview, F02 enrich volcanoes, F03 clustering), each emitting its
composite, the individual panels, and one supplementary workbook.

## Study design

Clean 2×2 factorial (PHE hypertrophic stress × mitochondrial transplant), n = 6/group,
24 samples, 4,806 proteins (DIA-MS). Contrasts and their narrative role:

| Contrast (code) | Definition | Role |
|---|---|---|
| Disease (`CTLvPHE`) | PHE − Ctl | insult signature |
| Rescue (`PHEvPHE_MITO`) | PHE_Mito − PHE | **primary** — reversal |
| Transplant (`CTLvMITO`) | Mito − Ctl | intervention in healthy cells |
| Interaction | (PHE_Mito − Mito) − (PHE − Ctl) | orthogonal 2×2 (underpowered) |
| Secondary (`MITOvPHE_MITO`) | PHE_Mito − Mito | exploratory |

Significance: **Π-score < 0.05** (Xiao 2014, `Π = P.Value^|log2FC|`) primary; BH-FDR < 0.10 secondary.
Spine = **Disease → Rescue reversal**.

## Layout

```
04_Figures/
  functions/                 shared engine code (pipeline-ordered; see below)
  shared/                     fgsea + rat gene-set caches read by the figures
  F05_modules/c_data/         cached WGCNA network read by F03
  F01_Proteome_Overview/      PCA, DEP, effect-size, overlap, pathway-count panels + composite
  F02_Enrich_Volcanoes/       4 enrichment rings (enrichVolcano) + 2x2 composite
  F03_Clustering/             WGCNA module figure: counts · heatmap · trajectory · ORA
  README.md                   this file
```

Three self-contained figures. F01 holds the proteome-overview panel builders
(`_build_pca.R`, `_build_dep.R`, `_build_venn.R`, `_build_pathway.R`) and assembles
them; F02 builds the enrichment rings directly from the DE + fgsea tables via the
`enrichVolcano` package; F03 is the WGCNA module figure. Each writes its composite
to `b_reports/main/png/`, the individual panels to `b_reports/panels/`, and the
panel data to a single `c_data/F0x_supplementary.xlsx`.

Each standalone figure subdir:

```
F01_Proteome_Overview/
  a_script/01_proteome_overview.R    composite driver
  a_script/_build_{pca,dep,venn,pathway}.R   panel builders
  b_reports/main/png/                MAIN_F01_proteome_overview.png
  b_reports/panels/                  panel_a_pca ... panel_f_pathway
  c_data/F01_supplementary.xlsx      PCA scores, DE tables, overlap membership, pathway counts
F02_Enrich_Volcanoes/
  a_script/01_enrich_volcanoes.R     4 enrichment rings (Disease/Transplant/Rescue/Interaction)
  b_reports/main/png/                MAIN_F02_enrich_volcanoes.png   (2x2 composite)
  b_reports/panels/                  MAIN_F02_{disease,transplant,rescue,interaction}_ring.png
  c_data/F02_supplementary.xlsx      contrast key + every tested pathway per contrast
F03_Clustering/
  a_script/01_clustering.R           WGCNA module figure (single row-aligned block)
  b_reports/main/png/                MAIN_F03_clustering.png
  b_reports/rat_go_bpccmf_sets.rds   cached GO BP/CC/MF sets for the per-module ORA
  c_data/F03_supplementary.xlsx      module counts, eigengene-limma, trajectory, ORA, hubs
```

## functions/ — shared engines (pipeline order)

The validated enrichment/style engine, named by where each piece occurs in a
figure's pipeline. The derived-data caches the figures read travel with the
suite: `shared/rat_gene_sets.rds`, `shared/fgsea_tstat_all_h9c2.csv`, and
`F05_modules/c_data/wgcna_network.rds`.

F02's ring panels and composite are drawn by the external `enrichVolcano` package
via `volcano_ring()` / `volcano_ring_grid()`, fed the tidy DE table and the fgsea
cache; there is no in-suite ring builder.

| File | Provides |
|---|---|
| `01_style_palettes_theme.R` | `FIG_THEME`, group/contrast palettes + `GROUP_LABELS`, thresholds, `scale_text()`, `clean_pathway_name()`, `fmt_p()` |
| `02_data_paths_and_loaders.R` | `P05` input paths, `load_combined_wide()`, contrast naming maps, `contrast_brief()`; sources 01 |
| `03_pathway_enrichment_dedup_ora.R` | `deduplicate_enrichment()`, `build_harmonized_collection()`, `clean_display_label()`, `classify_database()` |
| `04_mitocarta_lens_lookup.R` | `MITO_PATHWAY_REGEX` (mito-keyword lens) |
| `06_supplementary_workbook.R` | `build_workbook(out_file, figure_title, sheet_specs)` — writes an Overview sheet first, then data sheets; deletes nothing |
| `08_composite_layout.R` | `add_tag()`, `composite_caption()`, `save_composite()` — assembles patchwork composites and writes the PNG |

## Inputs (read-only)

From upstream pipeline stages:

- `02_Normalization/imputation/c_data/DAList_imputed_missforest.rds` — `$data` 4806×24 (rownames = uniprot_id), `$annotation` (uniprot_id→gene), `$metadata` (Col_ID, Group, Replicate).
- `03_DEP/a_non_imputed/c_data/combined_results_pi.csv` — long DE table (per protein × contrast).

Caches that travel with the suite:

- `shared/fgsea_tstat_all_h9c2.csv` — precomputed fgsea (limma moderated-t) for all 5 contrasts × 5 databases.
- `shared/rat_gene_sets.rds` — rat gene sets (Hallmark / Reactome / KEGG / MitoCarta / GO Slim).
- `F05_modules/c_data/wgcna_network.rds` — signed WGCNA network (eigengenes, kME, module colours) for F03.

## Conventions

- Root resolved via `here::here()`. Width invariant 178 mm (`PANEL_MD`); figures export a single PNG at 300 dpi, `units = "mm"`.
- Each figure's tabular output lives **only** in its `F0x_supplementary.xlsx` (Overview sheet documents every other sheet: name, role in the figure, contents).
- No edits to upstream stages (00–03) or to the cached inputs. No new gene-set derivation.
- **Group colours** (Ctl=#3B7DB5, Mito=#009E73, PHE=#E08214, PHE_Mito=#8073AC) match `GROUP_COLORS` in `01_style_palettes_theme.R`. The four conditions display as CTL / MitoTx / Phe-only / Phe+MitoTx via `GROUP_LABELS` (legends, trajectory axis, contrast math); the internal tokens (Ctl/Mito/PHE/PHE_Mito) stay unchanged for data joins. Contrast **role names** (Disease / Transplant / Rescue / Secondary / Interaction) come from `contrast_brief()` (backed by `CONTRAST_DISPLAY_MAP`) and are used on axis labels, bar fills, and panel titles throughout.
- **Direction palette**: `DIR_COLORS` (Up=#D6604D / Down=#4393C3 / NS=grey70) for all DEP direction plots.

## Run

```sh
Rscript 04_Figures/F01_Proteome_Overview/a_script/01_proteome_overview.R
Rscript 04_Figures/F02_Enrich_Volcanoes/a_script/01_enrich_volcanoes.R
Rscript 04_Figures/F03_Clustering/a_script/01_clustering.R
```

The three figures are independent; run order does not matter. Benign X11/cairo
warnings on stderr are expected on machines without XQuartz and do not affect
the PNG output.
