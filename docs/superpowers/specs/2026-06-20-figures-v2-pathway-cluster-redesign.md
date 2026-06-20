# 04_Figures_v2 — Pathway bars + Cluster pilots redesign

Replaces:
- `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R` (current NES-bar 2×2)
- `04_Figures_v2/06_Cluster/a_script/01_main_panels.R` (current Group_Pi / Group_All / Sample_All trio)

Upstream stats unchanged. Reads existing artifacts: `combined_results_pi.csv`, `DAList_imputed_missforest.rds`, `fgsea_tstat_all_h9c2.csv`, `rat_gene_sets.rds`, `wgcna_network.rds`.

## Doc-rot fix (independent)

`04_Figures_v2/README.md` and `04_Figures_v2/BUILD_PROMPT.md` say `DAList_imputed_imp4p.rds`. The actual loader (`functions/02_data_paths_and_loaders.R:48`) uses `DAList_imputed_missforest.rds`. Update both docs to match the code.

## F04 — Pathway bars (Panel-D-style count summary)

Single panel, no NES bars. One figure file.

- **Pool**: `CANONICAL_DBS` = Hallmark + Reactome + KEGG + MitoCarta + GO Slim.
- **Filter**: `padj < 0.05`, `size ≥ 10`, drop `MITO_DROP_SETS` (compartment aggregates).
- **Dedup**: `deduplicate_enrichment(jaccard_cutoff = 0.5*overlap + 0.5*jaccard ≥ 0.375, cross_db = TRUE)` per contrast (existing engine).
- **Mito flag**: pathway is mito if `database == "MitoCarta"` OR `grepl(MITO_PATHWAY_REGEX, pathway)` (existing `04_mitocarta_lens_lookup.R`).
- **Form**: per contrast, dodged Up bar + Down bar. Each bar drawn from 0 in light fill (= total), then mito subset overdrawn in dark fill (= mito). Numeric total label above bar. sqrt y. Contrast-tinted background bands. 4-key inline legend (Up total / Up mito / Down total / Down mito).
- **Output**: `b_reports/main/{pdf,png}/MAIN_F04_pathway_bars.{pdf,png}` (~120 × 70 mm).
- **Workbook** (`c_data/F04_supplementary.xlsx`): Overview sheet + `dep_pathway_counts` + per-contrast `<role>_sig_pathways` (post-dedup pathway, database, NES, padj, is_mito).
- **Delete** the existing `c_data/shown_pathways.csv` writeback — F06 no longer excludes prior figures' pathways, so the file has no consumer.

## F06 — Cluster pilots (6 pilots)

Six self-contained pilot figures. All share the per-cluster row layout: **trajectory left | Hallmark ORA bars right**, color-keyed to cluster identity, one row per cluster/module/quadrant. Patchwork widths `c(1, 1.4)`.

### Shared rules

- Gene-level matrix via `limma::avereps(dal$data, ID = gene)`; drop NA/empty gene rows.
- Group means matrix = column-collapsed by `H9C2_GROUP_LEVELS` (Ctl, Mito, PHE, PHE_Mito).
- ORA: per-DB `fgsea::fora` with per-DB BH (existing `run_ora_deduplicated`), then cross-DB dedup at 0.375. **Hallmark only on the main panels.** GO Slim swap → supp. No F04/F05 exclusion.
- ORA bar: top 6 by `padj` per cluster, bar fill = cluster identity color, x = `-log10(padj)`, label = `clean_display_label()`.
- Trajectory: standardized expression (z per gene), one line per gene with `alpha`/`linewidth` ~ membership (for c-means); for WGCNA/k-means use `alpha = 0.3` fixed; mean-by-cluster overlay as a thick line in cluster color. x-axis = group levels for group-based pilots, RRHO2 quadrants use a different layout (see below).
- Cluster identity colors: `viridis::turbo` interpolated to whatever cluster count comes out.
- Output naming: `b_reports/main/{pdf,png}/MAIN_F06_<pilot_key>.{pdf,png}` per pilot.
- Workbook: one `c_data/F06_supplementary.xlsx`, Overview sheet enumerates pilots, then per-pilot `<key>_membership` and `<key>_ora` sheets. Cluster-selection diagnostics for the c-means pilots go to `b_reports/supp/`.

### Pilot 1 — `pilot_p` (fuzzy c-means, group means)

- Filter: `p < 0.05` in ≥1 core contrast (`CTLvMITO`, `CTLvPHE`, `PHEvPHE_MITO`, `Interaction`).
- Engine: `e1071::cmeans` on standardized (per-gene z) 4-D group-mean matrix. Fuzzifier `m` = Schwämmle–Jensen 2010. **Fixed c = 6**, `set.seed(42)`.
- Cluster ordering: by mean group-profile shape similarity (hierarchical on cluster centroids).
- Supp: Dmin-vs-c sweep (c = 2..12, 5 seeds), `b_reports/supp/MAIN_F06_pilot_p_cluster_selection.png`.

### Pilot 2 — `pilot_pi` (same as Pilot 1, gate = `Π < 0.05`)

Identical engine and layout. Filter = `pi_score < 0.05` in ≥1 core contrast.

### Pilot 3 — `pilot_fdr` (same as Pilot 1, gate = `adj.P.Val < 0.10`)

Identical engine and layout. Filter = `adj.P.Val < 0.10` in ≥1 core contrast.

### Pilot 4 — `pilot_wgcna` (modules from existing build)

- Source: `04_Figures/F05_modules/c_data/wgcna_network.rds` (signed Pearson, n = 24, already built).
- "Clusters" = modules (drop the `grey` module = unassigned).
- Module eigengenes (MEs) already in the network object; correlate each ME with four binary contrast indicator vectors (Disease, Transplant, Rescue, Interaction) using Pearson r + Student p (samples = 24).
- Per module, derive a sign pattern from (Disease r, Rescue r): `Reversal` (opposite signs, both \|r\| ≥ 0.4), `Concordant up` (both positive), `Concordant down` (both negative), `Other`. Rows ordered: Reversal modules first, then Concordant, then Other.
- Row layout: trajectory (group means; sample-level profiles add noise that obscures module shape — implemented as 4-condition group-mean z-profiles for visual consistency with Pilots 1–3), ORA bars right. A small row header annotates module label, n_genes, ME×Disease r, ME×Rescue r, sign pattern.
- Supp: full ME×trait correlation heatmap (modules × 4 contrasts × indicators).

### Pilot 5 — `pilot_logfc` (k-means on per-protein 4-D logFC vector)

- Vector per protein: `(logFC_CTLvPHE, logFC_PHEvPHE_MITO, logFC_CTLvMITO, logFC_Interaction)` from `load_combined_wide()`. Drop rows with any NA in those 4 columns.
- Standardization: none (logFCs are already on a common scale).
- Engine: `kmeans(centers = 6, nstart = 50, iter.max = 100)`, `set.seed(42)`.
- Cluster labeling: derive geometric label from each cluster's `(mean_logFC_Disease, mean_logFC_Rescue)` quadrant — Reversed Down (D+, R−), Reversed Up (D−, R+), Concordant Up (D+, R+), Concordant Down (D−, R−), Neutral (close to origin), Other.
- Row layout: trajectory is **bar of mean logFC per contrast** (not group-mean line — because the basis is logFC vector, not group means), ORA bars right. Cluster header annotates geometric label + n_genes.
- Supp: silhouette plot for k = 2..10 + within-cluster sum of squares elbow.

### Pilot 6 — `pilot_rrho2` (threshold-free Disease↔Rescue concordance map)

- Inputs: two per-protein signed ranking vectors — Disease (`CTLvPHE`) and Rescue (`PHEvPHE_MITO`) — using signed `-log10(P.Value) * sign(logFC)` (the standard RRHO ranking; gives equal weight to direction and significance).
- Engine: `RRHO2::RRHO2_initialize(list1, list2, stepsize = ..., labels = c("Disease","Rescue"), boundary = 0.025, log10.ind = TRUE)`. Stepsize auto-set to ~`sqrt(n)` (RRHO2 default). `set.seed(42)`.
- **Main figure layout (differs from pilots 1–5 because RRHO2's headline is a heatmap):**
  - Panel A: RRHO2 rank-rank heatmap (`RRHO2_heatmap()`), with the four corner quadrants labeled UU (concordant up), DD (concordant down), UD (reversed: Disease up / Rescue down), DU (reversed: Disease down / Rescue up).
  - Panels B–E (right column, per-quadrant rows): for each of the four significant-overlap gene sets returned by RRHO2 (`hypermat.by.signs`-derived `genelist_uu/dd/ud/du`), one row showing trajectory (group means) of those genes + Hallmark ORA bars right. Color per quadrant: reversal quadrants = warm/cool reversal pair; concordant = green tones.
- Output: `MAIN_F06_pilot_rrho2.{pdf,png}` (~178 × 200 mm).
- Workbook sheets: `pilot_rrho2_genelist_<quadrant>` (the 4 gene lists), `pilot_rrho2_ora` (per-quadrant ORA).

## What gets deleted

- All NES-bar code in `04_Pathway_bars/a_script/01_main_panels.R` (no supp move).
- `Group_All` and `Sample_All` variant code in `06_Cluster/a_script/01_main_panels.R`.
- `06_Cluster/c_data/MAIN_F06_*` orphans from prior runs (cleaned in the new script's startup).

## Run order

```sh
Rscript 04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R
Rscript 04_Figures_v2/06_Cluster/a_script/01_main_panels.R   # produces 6 pilot figures + supp + workbook
```

F06 no longer reads `shown_pathways.csv` from either F04 or F05. F05's writeback stays untouched (F05 is out of scope).

## Out of scope

- 01_PCA, 02_DEP_bars, 03_Venn, 05_Enrich_Volcano (audited clean).
- Upstream stages 00–03 (unchanged).
- Pipeline-level imputation choice (missForest stays).
- The `04_Figures/F05_modules/` knowledge-driven WGCNA panel (separate figure family; this redesign only consumes its network artifact).

## Verification

After implementation:
- All 7 figure files (F04 + 6 pilots) regenerate without warnings (`Rscript` exit 0).
- Each figure visibly differs by pilot (cluster counts and sign patterns differ across gates).
- Each pilot's workbook sheet count matches the pilot description above.
