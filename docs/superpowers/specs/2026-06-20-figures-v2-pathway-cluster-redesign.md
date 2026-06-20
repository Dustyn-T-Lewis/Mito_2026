# 04_Figures_v2 — Pathway bars + Cluster pilots redesign

> **Spec status (last revised 2026-06-20 post-audit):** sections marked
> *Updated:* reflect a divergence between the original 2026-06-20 spec and the
> implemented code; the code is now authoritative and the spec has been
> rewritten to match. Each *Updated* note carries the commit hash that
> introduced the change so the rationale is recoverable.

Replaces:
- `04_Figures_v2/04_Pathway_bars/a_script/01_main_panels.R` (current NES-bar 2×2)
- `04_Figures_v2/06_Cluster/a_script/01_main_panels.R` (current Group_Pi / Group_All / Sample_All trio)

Upstream stats unchanged. Reads existing artifacts: `combined_results_pi.csv`, `DAList_imputed_missforest.rds`, `fgsea_tstat_all_h9c2.csv`, `rat_gene_sets.rds`, `wgcna_network.rds`.

## File layout (post-cleanup, 2026-06-20)

The figure scripts each consume a small set of `functions/` helpers via a
compatibility-stub pattern so single-responsibility files stay ≤200 lines
without breaking callers:

- `functions/03_pathway_enrichment_dedup_ora.R` — stub sourcing the three
  pieces below
  - `functions/03a_dedup_engine.R` — dedup primitives + `fora_odds_ratio()`
  - `functions/03b_enrichment_runners.R` — gene-set collection, fgsea/fora
  - `functions/03c_pathway_display.R` — pathway display names and category
    classifiers
- `functions/05_volcano_ring_plot_builder.R` — stub sourcing 05a-05e
  (label cleaner / ring data / volcano layers / ring layers / composite)
- `functions/07_cluster_row_layout.R` — shared cluster-row helpers; uses
  `fora_odds_ratio()` from 03a (transitive source order:
  panel script → 02 → 03 (→ 03a/b/c) → 06 → 07)
- `06_Cluster/a_script/01_main_panels.R` — driver for F06; runs setup +
  pilots 1-3 (c-means), then sources `02_pilot_wgcna.R`, `03_pilot_logfc.R`,
  `04_pilot_rrho2.R` siblings before assembling the workbook

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

Six self-contained pilot figures. All share the per-cluster row layout:
**trajectory left | Hallmark ORA middle | MitoCarta ORA right**, color-keyed to
cluster identity, one row per cluster/module/quadrant. Patchwork widths
`c(1, 1, 1)` (was `c(1, 1.4)` when only Hallmark was rendered).

### Shared rules

- Gene-level matrix via `limma::avereps(dal$data, ID = gene)`; drop NA/empty gene rows.
- Group means matrix = column-collapsed by `H9C2_GROUP_LEVELS` (Ctl, Mito, PHE, PHE_Mito).
- ORA: per-DB `fgsea::fora` with per-DB BH (existing `run_ora_deduplicated`), then cross-DB dedup at 0.375. *Updated (commit `d0eb731`):* **Hallmark + MitoCarta side-by-side on the main panels**, paired per the two-source guidance in Reimand 2019 (PMID 30664679); MitoCarta gives compartment-level resolution (Malate-Aspartate Shuttle, Mitoribosome, etc.) that Hallmark alone misses for a mito-transplant study. GO Slim swap → supp. No F04/F05 exclusion.
- ORA bar: top 6 by `padj` per cluster per DB, bar fill = cluster identity color, x = `-log10(padj)`, label = `clean_display_label()`.
- Trajectory: standardized expression (z per gene), one line per gene with `alpha`/`linewidth` ~ membership (for c-means); for WGCNA/k-means use `alpha = 0.3` fixed; mean-by-cluster overlay as a thick line in cluster color. x-axis = group levels for group-based pilots, RRHO2 quadrants use a different layout (see below).
- Cluster identity colors: `viridis::turbo` interpolated to whatever cluster count comes out.
- Output naming: `b_reports/main/{pdf,png}/MAIN_F06_<pilot_key>.{pdf,png}` per pilot. The RRHO2 heatmap saves as a separate file `MAIN_F06_pilot_rrho2_heatmap.{pdf,png}` (base-graphics output that cannot share a patchwork canvas with the ggplot per-quadrant rows).
- Workbook: one `c_data/F06_supplementary.xlsx` with two index sheets — `Overview` (auto-generated data dictionary listing every sheet + its role + its columns) and `Pilot_summary` (per-pilot method, gate, gene count, fuzzifier m, chosen c, selection basis) — then per-pilot data sheets. *Updated (commit `b47e2b5`):* the per-pilot summary was originally also named `Overview`, which collided with the auto-generated dictionary; renamed to `Pilot_summary`. Cluster-selection diagnostics for the c-means and k-means pilots go to `b_reports/supp/`.

### Pilot 1 — `pilot_p` (fuzzy c-means, group means)

- Filter: `p < 0.05` in ≥1 core contrast (`CTLvMITO`, `CTLvPHE`, `PHEvPHE_MITO`, `Interaction`).
- Engine: `e1071::cmeans` on standardized (per-gene z) 4-D group-mean matrix. Fuzzifier `m` = Schwämmle–Jensen 2010 (PMID 20880957). `set.seed(42)`.
- *Updated (commit `8392a4e`):* **per-pilot c via Dmin elbow over c = 2..⌊√(N/2)⌋**. Drop-rule: first c where the marginal drop to the next c falls below 10% of the Dmin curve range. The √(N/2) cap is the Mardia, Kent & Bibby (1979) heuristic for cluster-count upper bound. Fallback to `DEFAULT_C = 6L` only if no flattening is detected. Observed: pilot_p c=3 (cap=29), pilot_pi c=4 (cap=9), pilot_fdr c=3 (cap=13). The original fixed `c=6` over-clustered the smaller gates.
- Cluster ordering: by mean group-profile shape similarity (hierarchical on cluster centroids).
- Supp: Dmin-vs-c sweep (c = 2..12, 5 seeds), `b_reports/supp/MAIN_F06_pilot_p_selection.{pdf,png}`.

### Pilot 2 — `pilot_pi` (same as Pilot 1, gate = `Π < 0.05`)

Identical engine and layout. Filter = `pi_score < 0.05` in ≥1 core contrast.

### Pilot 3 — `pilot_fdr` (same as Pilot 1, gate = `adj.P.Val < 0.10`)

Identical engine and layout. Filter = `adj.P.Val < 0.10` in ≥1 core contrast.

### Pilot 4 — `pilot_wgcna` (modules from existing build)

- Source: `04_Figures/F05_modules/c_data/wgcna_network.rds` (signed Pearson, n = 24, already built — Langfelder & Horvath 2008, PMID 19114008).
- "Clusters" = modules (drop the `grey` module = unassigned).
- Module eigengenes (MEs) already in the network object; correlate each ME with four binary contrast indicator vectors (Disease, Transplant, Rescue, Interaction) using Pearson r + Student p (samples = 24). *Updated (commit `b47e2b5`):* the Interaction indicator is `c(Ctl = +1, Mito = -1, PHE = -1, PHE_Mito = +1)` to match the limma contrast `(PHE_Mito − Mito) − (PHE − Ctl) = +Ctl − Mito − PHE + PHE_Mito` (the original spec had the sign reversed in the supp ME-trait sheet — Reversal classification was unaffected because it keys only off Disease and Rescue, but the supp column was sign-flipped).
- Per module, derive a sign pattern from (Disease r, Rescue r): `Reversal` (opposite signs, both \|r\| ≥ 0.25), `Concordant up` (both positive), `Concordant down` (both negative), `Other`. Rows ordered: Reversal modules first, then Concordant, then Other. *Updated (commits `6e113d4` → `d7b2f3b`):* the floor was originally `|r| ≥ 0.40`; at n = 24 that floor collapses every module into "Other" (Student p for r = 0.40 at n = 24 is ~p = 0.05; for r = 0.25 it is ~p = 0.24). The 0.25 floor surfaces the documented 4 Reversal + 1 Concordant Down cohort (blue, brown, black, turquoise + red). **Treat the sign_pattern column as a ranking aid rather than a per-module inferential claim**, and propagate this hedge into any figure caption that quotes module assignments.
- Row layout: trajectory (group means; sample-level profiles add noise that obscures module shape — implemented as 4-condition group-mean z-profiles for visual consistency with Pilots 1–3), Hallmark + MitoCarta ORA bars to the right of the trajectory. A small row header annotates module label, n_genes, ME×Disease r, ME×Rescue r, ME×Transplant r, sign pattern.
- Supp: full ME×trait correlation heatmap (modules × 4 contrasts × indicators).

### Pilot 5 — `pilot_logfc` (k-means on per-protein 4-D logFC vector)

- Vector per protein: 4 columns of `load_combined_wide()` keyed by `H9C2_CORE_CONTRASTS` order = `(logFC_CTLvMITO, logFC_CTLvPHE, logFC_PHEvPHE_MITO, logFC_Interaction)`. *Updated (canonical-order fix):* original spec listed the 4 contrasts in a different order; the code uses the canonical `H9C2_CORE_CONTRASTS` ordering used everywhere else in the pipeline. k-means is order-invariant on the vector so cluster geometry is unaffected, but centroid columns are labeled in CORE order.
- Drop rows with any NA in those 4 columns.
- Standardization: none (logFCs are already on a common scale).
- Engine: `kmeans(centers = k, nstart = 50, iter.max = 100)`, `set.seed(42)`.
- *Updated (commit `8392a4e`):* **k chosen by gap statistic firstSEmax** (Tibshirani, Walther & Hastie 2001, doi:10.1111/1467-9868.00293) via `cluster::clusGap(B = 50)` — smallest k whose gap is within 1 SE of the maximum. Fallback to `DEFAULT_C = 6L` if `clusGap` fails. The original fixed `k = 6` happens to match the current gap-stat output for this dataset, but the mechanism is now data-driven rather than hard-coded.
- Cluster labeling: derive geometric label from each cluster's `(mean_logFC_Disease, mean_logFC_Rescue)` quadrant — Reversed Down (D+, R−), Reversed Up (D−, R+), Concordant Up (D+, R+), Concordant Down (D−, R−), Neutral (close to origin), Other. The (Disease, Rescue) axes use `logFC_CTLvPHE` and `logFC_PHEvPHE_MITO` regardless of where they fall in the contrast vector.
- Row layout: trajectory is **bar of mean logFC per contrast** (not group-mean line — because the basis is logFC vector, not group means), Hallmark + MitoCarta ORA bars right. Cluster header annotates geometric label + n_genes.
- Supp: gap-statistic diagnostic — gap(k) ± 1 SE bars across k = 2..10 with the chosen k highlighted.

### Pilot 6 — `pilot_rrho2` (threshold-free Disease↔Rescue concordance map)

- Inputs: two per-protein signed ranking vectors — Disease (`CTLvPHE`) and Rescue (`PHEvPHE_MITO`). *Updated (commit `ed28a1e`):* **ranking metric = signed limma moderated-t** (`cw$t_CTLvPHE`, `cw$t_PHEvPHE_MITO`) rather than the original signed `-log10(P.Value) * sign(logFC)`. Justifications: (a) variance-stabilized at n = 24 (Smyth 2004, doi:10.2202/1544-6115.1027), (b) matches the statistic fgsea consumes upstream so RRHO2 and fgsea agree on which gene "leads" each contrast, (c) RRHO2's hypergeometric tail depends on rank order — not metric value — so the choice is essentially equivalent on signal quadrants but better-behaved for the small-overlap tail.
- Engine: `RRHO2::RRHO2_initialize(list1, list2, stepsize = ceiling(sqrt(n)), labels = c("Disease (signed mod-t)", "Rescue (signed mod-t)"), boundary = 0.025, log10.ind = TRUE)`. `set.seed(42)` (RRHO2_initialize is deterministic; seed is defensive).
- **Main figure layout (differs from pilots 1–5 because RRHO2's headline is a heatmap):**
  - Panel A: RRHO2 rank-rank heatmap (`RRHO2_heatmap()` from base graphics, saved as a **separate file** `MAIN_F06_pilot_rrho2_heatmap.{pdf,png}` — base graphics cannot share a patchwork canvas with the ggplot per-quadrant rows). Outer-margin labels at the four corners: UU (concordant up, top-left), DD (concordant down, bottom-right), UD (reversed: Disease up / Rescue down, top-right), DU (reversed: Disease down / Rescue up, bottom-left). The labels are drawn with `mtext(outer = TRUE)` because `RRHO2_heatmap` leaves the device active on its color-bar panel.
  - Panels B–E (`MAIN_F06_pilot_rrho2.{pdf,png}`, per-quadrant rows): for each of the four quadrants, one row with trajectory (group means) of the quadrant gene set + Hallmark ORA + MitoCarta ORA bars. Color per quadrant: UU green, DD blue, UD warm-red (reversed), DU coral (reversed).
- *Updated (added at implementation time, commit `b47e2b5`):* **Top-20% percentile fallback for sparse quadrants.** When the RRHO2 peak-overlap set returned by `genelist_<q>$gene_list_overlap_<q>` is < 5 genes (common with sparse concordance in 4-group n = 24 designs), the quadrant gene set is filled from the intersection of the top-20% rank fraction of both lists in the quadrant's implied direction. The fallback is tagged per gene in the workbook (`source ∈ {rrho2_peak, pct_fallback}`) and flagged in the figure subtitle so reviewers can tell which rows bypass the hypergeometric significance test. **The fallback yields a biologically meaningful gene set but is not an RRHO2-significant overlap**; manuscript captions should inherit this hedge.
- Output: `MAIN_F06_pilot_rrho2.{pdf,png}` (~178 × 200 mm) + `MAIN_F06_pilot_rrho2_heatmap.{pdf,png}` (separate).
- Workbook sheets: `pilot_rrho2_genelists` (single long-format sheet with `quadrant, role, gene, source` for all 4 quadrants), `pilot_rrho2_ora` (per-quadrant Hallmark ORA), `pilot_rrho2_ora_mito` (per-quadrant MitoCarta ORA).

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

Verifier scripts live under `/tmp/` (developer-local; the verifier set is
small and self-contained, so it doesn't need to be in the repo):

- `/tmp/verify_F04.R` — checks F04 PDF/PNG present, workbook has the 6 expected sheets (`Overview`, `dep_pathway_counts`, `<contrast>_sig_pathways` × 4).
- `/tmp/verify_F06_rrho2.R` — sources `verify_F06_logfc.R` which sources `verify_F06_wgcna.R` which sources `verify_F06_cmeans.R`, so a single `Rscript /tmp/verify_F06_rrho2.R` exercises all four F06 sub-pilots in one call. Each verifier checks PDF/PNG presence + minimum size + the pilot's expected workbook sheets.

After implementation, the following must hold (last confirmed 2026-06-20 on commit `0157237`):

- All 8 figure files regenerate cleanly (`Rscript` exit 0):
  - F04 main PNG/PDF
  - F06 × 6 main PNG/PDF (one per pilot: `pilot_p`, `pilot_pi`, `pilot_fdr`, `pilot_wgcna`, `pilot_logfc`, `pilot_rrho2`)
  - F06 RRHO2 heatmap (separate PNG/PDF: `MAIN_F06_pilot_rrho2_heatmap`)
- All 5 verifier `PASS` lines printed.
- F04 produces 83 total post-dedup significant pathways across 4 contrasts (Rescue 52 > Transplant 14 > Disease 10 > Interaction 7, mito subset 14 / 5 / 0 / 2).
- F06 produces a 22-sheet workbook: `Overview` + `Pilot_summary` + per-pilot `*_membership` + per-pilot `*_ora` + per-pilot `*_ora_mito` (plus pilot-specific extras: `pilot_wgcna_me_traits`, `pilot_logfc_centroids`, `pilot_rrho2_genelists`).
- RRHO2 quadrant sizes for this dataset: UU = 35 (pct_fallback), DD = 36 (pct_fallback), UD = 1524 (rrho2_peak), DU = 1160 (rrho2_peak).
- WGCNA at `r_floor = 0.25` surfaces 4 Reversal modules (blue, brown, black, turquoise) + 1 Concordant down (red); remaining modules classified Other.
- `pilot_logfc` cluster 3 (n = 528) carries the "Reversed Up" label with Hallmark OXPHOS + Adipogenesis at the top and MitoCarta Malate-Aspartate Shuttle + Mitoribosome + SLC25A small-molecule transport leading the compartment ORA.
