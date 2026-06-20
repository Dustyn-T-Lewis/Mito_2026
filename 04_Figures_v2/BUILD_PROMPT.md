# BUILD_PROMPT — regenerate the `04_Figures_v2` suite from scratch

This prompt re-specifies all six figures (inputs, engine calls, style invariants,
outputs) precisely enough for an agent to rebuild the suite. Reuse the existing
engines in `04_Figures_v2/functions/` — do **not** re-derive style, palettes, gene
sets, or enrichment logic. No edits to stages 00–03 or to `04_Figures/`.

## Global contract

- **Root**: `here::here()` → project root. There is no `setup.R`; use `here::here()`.
- **Source engines** (pipeline order) from `04_Figures_v2/functions/`:
  `01_style_palettes_theme.R`, `02_data_paths_and_loaders.R`,
  `03_pathway_enrichment_dedup_ora.R`, `04_mitocarta_lens_lookup.R`,
  `05_volcano_ring_plot_builder.R`, `06_supplementary_workbook.R`.
  `02` transitively sources `01`. Each script sources only what it needs.
- **Data caches** stay read from `04_Figures/shared/`: `fgsea_tstat_all_h9c2.csv`,
  `rat_gene_sets.rds`, `protein_localization_rat.csv`. Do not copy them into v2.
- **Inputs**: `P05$imp_rds` (DAList_imputed_missforest.rds: `$data` 4806×24 rownames=uniprot_id,
  `$annotation` uniprot_id→gene, `$metadata$Group` ∈ {Ctl,Mito,PHE,PHE_Mito}, cols keyed by Col_ID);
  `P05$comb` (combined_results_pi.csv, long: uniprot_id, gene, logFC, t, P.Value, adj.P.Val,
  pi_score, sig_pi, contrast; 5 contrasts). `load_combined_wide()` → wide
  (`logFC_<ctr>`, `t_<ctr>`, `P.Value_<ctr>`, `pi_score_<ctr>`, …).
- **fgsea cache** columns: pathway, pval, padj, log2err, ES, NES, size, leadingEdge (";"-joined),
  database (Hallmark/Reactome/KEGG/MitoCarta/GO Slim), contrast. Precomputed on limma
  moderated-t — never recompute.
- **Contrasts**: Disease=`CTLvPHE`, Rescue=`PHEvPHE_MITO` (primary), Transplant=`CTLvMITO`,
  `Interaction`, Secondary=`MITOvPHE_MITO`. Core set for most figures = the first four.
- **Thresholds**: `H9C2_PI_THRESH`=0.05 (primary), `H9C2_FDR_EXPLOR`=0.10.
- **Style invariants**: `FIG_THEME`; width 178 mm (`PANEL_MD`); pdf `device = get_pdf_device()`;
  png `dpi = 300`; `units = "mm"`. Palettes `H9C2_PAL_GROUP`, `H9C2_PAL_DIR`,
  `CONTRAST_COLORS`. Title 7pt bold, subtitle bold-italic grey30. Unicode minus in math.
- **Outputs**: standalone `MAIN_F0x_*.{pdf,png}` to `b_reports/main/{pdf,png}/`
  (diagnostics to `b_reports/supp/`). Tabular results go to a single
  `c_data/F0x_supplementary.xlsx` via `build_workbook(out_file, figure_title, sheet_specs)`,
  where each spec is `list(name=, df=, role=, contents=)` and the builder writes an
  **Overview** sheet first. Emit **no loose CSVs** except `shown_pathways.csv` in
  figure 05 (written for reference; F06 does **not** read it; F04 does not write it).

## F01 — 01_PCA

- Read `P05$imp_rds`; `imp_mat <- dal$data[, dal$metadata$Col_ID]`.
- `prcomp(t(imp_mat), center=TRUE, scale.=TRUE)`; scatter PC1/PC2, points colored+shaped
  by Group (`H9C2_PAL_GROUP`; shapes Ctl=16,Mito=17,PHE=15,PHE_Mito=18), 80% group ellipses.
- `vegan::adonis2` on `dist(scale(t(imp_mat)))`, 9999 perms, `by="terms"`: overall `~Group`
  + pairwise Disease (Ctl|PHE), Transplant (Ctl|Mito), Rescue (PHE|PHE_Mito). `vegan::betadisper`
  + `permutest` (999) dispersion check; warn if p<0.05. `set.seed(42)` before each.
- Annotate the PERMANOVA block (overall + 3 pairwise R²/p + betadisper p) in plot headroom.
- **Outputs**: `MAIN_F01_pca.{pdf,png}` (~120×110 mm); xlsx sheets `pca_scores`, `permanova`.

## F02 — 02_DEP_bars

- Long `P05$comb`; core 4 contrasts.
- Panel A: DEP counts as **% of proteome** at three independent thresholds (p<0.05 /
  FDR<0.10 / Π<0.05), Up/Down dodged, contrast-colored accent backgrounds, in-bar threshold cues.
- Panel B: signed log₂FC histogram+density faceted by contrast, per-facet median |log₂FC|.
- **Outputs**: `MAIN_F02_dep_bars.{pdf,png}` (178 mm wide); xlsx sheets `dep_counts` +
  one DEP table per contrast.

## F03 — 03_Venn

- 3 Π<0.05 membership sets by uniprot_id: Disease / Transplant / Rescue.
- Area-proportional 3-set Venn via `eulerr::euler()` (fills = contrast colors α≈0.5,
  quantities shown). Render through `ggplotify::as.ggplot()` if available, else to the device
  directly. Companion dodged Up/Down direction strip (Up `#D6604D`, Down `#4393C3`).
- **Outputs**: `MAIN_F03_venn.{pdf,png}` (~140×100 mm); xlsx sheets `membership`, `region_counts`.

## F04 — 04_Pathway_bars (Panel-D pathway count summary)

- Pool fgsea cache rows for `CANONICAL_DBS` = `{Hallmark,Reactome,KEGG,MitoCarta,GO Slim}`,
  `padj<0.05`, `size≥10`, drop `MITO_DROP_SETS` (MitoCarta compartment aggregates).
- Per contrast: `deduplicate_enrichment(df, pathways=set_pool, cutoff=0.375, cross_db=TRUE)`
  (EnrichmentMap-style combined Jaccard+overlap at 0.375 threshold).
- **Single panel, no NES bars.** Per contrast, two dodged bars (Up/Down) drawn from 0.
  Light fill = total significant pathways; dark fill = mito subset overdrawn.
  Mito flag: `database == "MitoCarta"` OR `grepl(MITO_PATHWAY_REGEX, pathway)`. Sqrt y.
  Contrast-tinted background bands. 4-key inline legend (Up total / Up mito / Down total / Down mito).
  Numeric total label above each bar. Names via `clean_display_label()`.
- **Outputs**: `MAIN_F04_pathway_bars.{pdf,png}` (~120×70 mm); xlsx `dep_pathway_counts` +
  per-contrast `<role>_sig_pathways` sheets. **Do not write** `shown_pathways.csv`
  (F06 no longer excludes prior-figure pathways; F04 removes any legacy copy).

## F05 — 05_Enrich_Volcano

- `dep_df <- load_combined_wide()`; canonical pool `{Hallmark,Reactome,KEGG,MitoCarta}`.
- For each of the 4 contrasts (Disease, Transplant, Rescue, Interaction): filter cache
  (`padj<0.05`, `size≥10`, no compartment sets) → `deduplicate_enrichment(jaccard=0.5, cross_db=TRUE)`
  → `pick_symmetric(n_each=6, per_db_cap=2)` → `build_ring_180_split()` →
  `make_volcano_ring(..., bg_color=CONTRAST_COLORS[ctr], show_legend=FALSE)`. Empty-ring fallback
  for contrasts with no FDR-significant terms.
- Save each ring standalone `MAIN_F05_<ctr_lower>_ring.{pdf,png}` (~110 mm square) plus shared
  `MAIN_F05_nes_legend.{pdf,png}` from `build_nes_legend_bar()`.
- **Outputs**: 4 rings + legend; xlsx `contrast_map` + per-contrast ring-term sheets;
  **keep** `c_data/shown_pathways.csv` (union of displayed ring terms).

## F06 — 06_Cluster (6-pilot cluster framework)

Six self-contained pilot figures. All share per-cluster row layout: **trajectory left |
Hallmark-only ORA bars right**, patchwork widths `c(1, 1.4)`. **Do not exclude** F04/F05
pathways — each pilot's ORA surfaces unfiltered biology.

- **Pilot 1 `pilot_p`** — fuzzy c-means on 4-D group means; gate `p < 0.05` in ≥1 core contrast.
  Engine: `e1071::cmeans`, `standardise` = per-gene z, Schwämmle–Jensen fuzzifier m,
  **fixed c = 6**, `set.seed(42)`. Supp: Dmin-vs-c sweep (c = 2..12, 5 seeds).
- **Pilot 2 `pilot_pi`** — same engine; gate `Π < 0.05` in ≥1 core contrast.
- **Pilot 3 `pilot_fdr`** — same engine; gate `adj.P.Val < 0.10` in ≥1 core contrast.
- **Pilot 4 `pilot_wgcna`** — WGCNA modules from `04_Figures/F05_modules/c_data/wgcna_network.rds`
  (drop grey). Correlate MEs with contrast indicator vectors. Row order: Reversal first, then
  Concordant, then Other. Trajectory: group-mean z-profiles (4 conditions); supp: ME×trait heatmap.
- **Pilot 5 `pilot_logfc`** — k-means (c = 6, nstart = 50) on per-protein 4-D logFC vector.
  Cluster label from `(mean_logFC_Disease, mean_logFC_Rescue)` sign quadrant. Trajectory: bar
  of mean logFC per contrast (not group-mean line). Supp: WSS elbow k = 2..10.
- **Pilot 6 `pilot_rrho2`** — RRHO2 threshold-free Disease↔Rescue map. Heatmap saved separately
  (base graphics); per-quadrant row layout (UU/DD/UD/DU) with group-mean trajectory + ORA.

ORA for all pilots: `run_hallmark_ora()` (Hallmark only, `fgsea::fora`). Top 6 by padj, bar
fill = cluster color, x = `-log10(padj)`.

- **Outputs**: `MAIN_F06_<pilot_key>.{pdf,png}` per pilot; supp diagnostics to `b_reports/supp/`;
  one `c_data/F06_supplementary.xlsx` — Overview sheet + per-pilot `_membership` and `_ora` sheets.
  No loose CSVs.

## Verify each script

`Rscript 04_Figures_v2/<fig>/a_script/01_main_panels.R`; confirm images regenerate and the
single Overview-first xlsx exists. F04 and F06 are independent — no ordering constraint.
