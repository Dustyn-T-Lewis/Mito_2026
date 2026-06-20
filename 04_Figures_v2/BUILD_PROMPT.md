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
  figures 04 & 05 (consumed by 06).

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

## F04 — 04_Pathway_bars (pathway-changes main report)

- Pool fgsea cache rows for `database ∈ {Hallmark,Reactome,KEGG,MitoCarta}`, `padj<0.05`,
  `size≥10`, drop MitoCarta compartment sets (`MITOCARTA_{IMM,IMS,MATRIX,OMM,ALL}`).
- Per contrast: `deduplicate_enrichment(df, pathways=set_pool, jaccard_cutoff=0.5, cross_db=TRUE)`
  where `set_pool <- c(rat_gene_sets$Hallmark, $Reactome, $KEGG, $MitoCarta)`.
- Top ≤6 Up / ≤6 Down by padj per contrast; horizontal bars faceted by contrast
  (`free_y`), fill = NES gradient `c("#08306B","#4393C3","white","#D6604D","#67000D")`,
  `values=rescale(c(-3,-1.5,0,1.5,3))`, `limits=c(-3,3)`, `oob=squish`. Names via `clean_display_label()`.
- **Outputs**: `MAIN_F04_pathway_bars.{pdf,png}` (178 mm wide); xlsx `pathway_bars` +
  per-contrast enrichment sheets; **keep** `c_data/shown_pathways.csv` (contrast, database, pathway).

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

## F06 — 06_Cluster (fuzzy c-means → per-cluster ORA)

- **Do not use the Mfuzz package** (it needs tcltk/X11). Reimplement on `e1071::cmeans`:
  `standardise` = per-gene z (`t(scale(t(mat)))`); `mestimate` = Schwämmle & Jensen 2010
  (`m = 1 + (1418/N+22.05)*D^-2 + (12.33/N+0.243)*D^(-0.0406*log(N)-0.1134)`, N=genes, D=dims);
  `mfuzz` = `cmeans(z, centers=c, m=m, method="cmeans")`, `set.seed(42)`.
- Gene-level matrix via `limma::avereps(dal$data, ID=gene)` (drop NA/empty genes).
- Three variants: **Group_Pi** (4 condition means, genes with Π<0.05 in ≥1 contrast),
  **Group_All** (4 condition means, all genes), **Sample_All** (24-sample profiles, all genes).
- **Cluster selection (report it)**: for c=2..12, mean Dmin (= `min(dist(centers))`) over ~5 seeds;
  pick the curvature knee (point where the per-step decrease flattens), default c=9 if no knee
  clears ~5% of the Dmin range. Emit a Dmin-vs-c diagnostic to `b_reports/supp/`. Report m and chosen c.
- Trajectory plots per cluster (membership-weighted alpha/linewidth, faceted). ORA on core
  members (membership>0.5) via `run_ora_deduplicated(genes, universe=all genes,
  pathways=build_harmonized_collection(), jaccard_cutoff=0.5)`. **Exclude** pathways listed in
  `04_Pathway_bars/c_data/shown_pathways.csv` and `05_Enrich_Volcano/c_data/shown_pathways.csv`
  (union of `pathway` IDs) so clusters surface new biology. Dot/bar of surviving terms per cluster.
- **Outputs**: `MAIN_F06_<variant>_traj.{pdf,png}` ×3, `MAIN_F06_<variant>_ora.{pdf,png}` ×3,
  `b_reports/supp/.../MAIN_F06_<variant>_cluster_selection.{pdf,png}` ×3; one
  `c_data/F06_supplementary.xlsx` with per-variant `_selection`/`_membership`/`_ora` sheets. No loose CSVs.

## Verify each script

`Rscript 04_Figures_v2/<fig>/a_script/01_main_panels.R`; confirm images regenerate and the
single Overview-first xlsx exists. Run 04 and 05 before 06.
