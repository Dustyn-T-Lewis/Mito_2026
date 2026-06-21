# Design — `04_Figures_v2` + parallel Mito manuscript

Date: 2026-06-19
Status: approved (design); spec under review before implementation

## 1. Goal

Produce, for the Mito (H9c2) proteomics project, a figure suite and an IMRAD prose
draft that parallel the YvO proteome-plasticity manuscript in narrative arc, statistical
treatment, and visual style. Figures are reorganized into one-plot-per-subdir under a new
`04_Figures_v2/` tree, reusing the existing `04_Figures/shared/` engines (no re-derivation
of style, palettes, gene sets, or enrichment logic).

## 2. Study framing (the parallel to YvO)

YvO: *does aging constrain the training response, and does training reverse aging?* on human muscle.
Mito: structurally identical pair on **H9c2 cardiomyoblasts**, a clean **2×2 factorial**
(PHE hypertrophic stress × mitochondrial transplant), **6 replicates/group, 24 samples,
4,806 proteins**.

| Mito contrast | Definition | YvO analogue | Role |
|---|---|---|---|
| Disease (`CTLvPHE`) | PHE − Ctl | Aging (insult signature) | insult |
| Rescue (`PHEvPHE_MITO`) | PHE_Mito − PHE | Training-in-Old (reversal) | **primary** |
| Transplant (`CTLvMITO`) | Mito − Ctl | intervention in healthy | intervention |
| Interaction | (PHE_Mito − PHE) − (Mito − Ctl) | Age×Training | orthogonal 2×2 |
| Secondary (`MITOvPHE_MITO`) | PHE_Mito − Mito | — | exploratory |

Significance: **Π-score < 0.05** (Xiao 2014, `Π = P.Value^|logFC|`) primary; **BH-FDR < 0.10** secondary.
Manuscript spine = **Disease → Rescue reversal**, mirroring YvO Aging → Reversal.

## 3. One-sentence narrative backbone

- **Question (technique):** DIA-MS quantifies the H9c2 proteome (4,806 proteins × 24 samples)
  to ask which proteins/pathways PHE hypertrophic stress remodels and whether mitochondrial
  transplant reverses that remodeling.
- **Test (hypothesis):** a 2×2 design tests a disease signature (PHE−Ctl), a rescue signature
  (PHE_Mito−PHE), the transplant effect (Mito−Ctl), and their interaction.
- **Prediction:** PHE imposes a coordinated mitochondrial/metabolic stress signature, and Mito
  transplant directionally opposes (rescues) a substantial fraction of disease-altered
  proteins/pathways rather than acting on an unrelated program.
- **Method:** filter → cyclic-loess normalize → imp4p impute; limma with replicate block via
  `duplicateCorrelation`, Π-score (FDR<0.10 secondary); fgsea/ORA on rat gene sets
  (Hallmark/Reactome/KEGG/GO Slim/MitoCarta); mfuzz trajectory soft-clustering.
- **Results:** a defined PHE disease proteome signature; measurable Mito reversal (rescue
  opposing disease); mfuzz programs whose per-cluster enrichment localizes stressed vs restored biology.
- **Discussion:** Mito transplant selectively counteracts the hypertrophic-stress proteome
  (not mere added mitochondrial mass), supporting targeted organelle transplantation and naming
  the rescued and resistant programs.

## 4. Directory layout

```
04_Figures_v2/
  shared/   -> reuse 04_Figures/shared/ (style.R, config.R, pathway_utils.R,
               volcano_ring.R, mitocarta_utils.R, go_slim_categories.R,
               rat_gene_sets.rds, goslim_rat_gene_sets.rds, fgsea_tstat_all_h9c2.csv,
               fgsea_goslim_h9c2.csv, protein_localization_rat.csv)
  01_PCA/            a_script/  b_reports/{main,supp}/{pdf,png}  c_data/
  02_DEP_bars/       "
  03_Venn/           "
  04_Pathway_bars/   "
  05_Enrich_Volcano/ "
  06_Cluster/        "
  BUILD_PROMPT.md
  README.md
```

`shared/` reuse: scripts resolve project root via `rprojroot::find_root(has_file("setup.R"))`
and `source()` the existing `04_Figures/shared/*` files directly (single source of visual/analytic truth).
Each figure: `a_script/01_main_panels.R` builds + exports standalone `MAIN_F0x_*.{pdf,png}`
(178 mm wide, `FIG_THEME`, `get_pdf_device()`), writes `c_data/*.csv` + `F0x_supplementary.xlsx`.
**No `90_stitch`, no composite — standalones only.**

## 5. Per-figure spec

### 01_PCA
- Input: `02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds` (`$data` 4806×24, `$metadata`).
- `prcomp(t(data), scale.=TRUE)`; group ellipses; `H9C2_PAL_GROUP`; per-group shapes.
- Annotate `vegan::adonis2` (Euclidean, 9999 perms): overall + pairwise Disease/Transplant/Rescue;
  `vegan::betadisper` dispersion check.
- Outputs: `MAIN_F01_pca.{pdf,png}`, `panel_A_pca.csv`, `panel_A_pairwise_permanova.csv`.
- Lifts `04_Figures/F01_QC_overview` panel A.

### 02_DEP_bars
- Input: `03_DEP/a_non_imputed/c_data/combined_results_pi.csv` (long).
- Per-contrast DEP counts as **% of proteome** at p<0.05 / FDR<0.10 / Π<0.05; Up/Down dodged;
  contrast-colored; nested thresholds. Companion |log₂FC| effect-size distribution per contrast.
- Outputs: `MAIN_F02_dep_bars.{pdf,png}`, `panel_dep_counts.csv`.
- Lifts F01 panel C.

### 03_Venn (replaces UpSet)
- Input: `combined_results_pi.csv`; Π<0.05 DEP membership.
- Area-proportional **3-set Venn** (`eulerr`): Disease ∩ Transplant ∩ Rescue. Companion
  Up/Down directional breakdown strip. Palette from `CONTRAST_COLORS`.
- Outputs: `MAIN_F03_venn.{pdf,png}`, `panel_venn_membership.csv`.

### 04_Pathway_bars
- **Same DB pool & logic as volcano rings**: Hallmark + Reactome + KEGG + MitoCarta;
  fgsea on limma moderated-t (cache `fgsea_tstat_all_h9c2.csv`); FDR<0.05;
  `deduplicate_enrichment(jaccard=0.5, cross_db=TRUE)`.
- Horizontal bars, top-N enriched pathways per contrast, fill = NES gradient
  (`#08306B→#4393C3→white→#D6604D→#67000D`, limits ±3), faceted by contrast. Cleaned names via
  `clean_pathway_name()`.
- Outputs: `MAIN_F04_pathway_bars.{pdf,png}`, `panel_pathway_bars.csv`. This is the "pathway changes" main report.

### 05_Enrich_Volcano
- Reuse `make_volcano_ring()` / `build_ring_180_split()` for **4 contrasts**: Disease, Transplant,
  Rescue, Interaction (secondary MITOvPHE_MITO dropped). Standalone ring plots + shared NES legend (`build_nes_legend_bar()`).
- DB pool = canonical pooled (Hallmark+Reactome+KEGG+MitoCarta), per-DB cap 2, ≤6 arcs/side.
- Outputs: `MAIN_F05_<contrast>_ring.{pdf,png}`, `MAIN_F05_nes_legend.{pdf,png}`, per-panel `*_ring_terms.csv`.

### 06_Cluster (mfuzz → per-cluster ORA)
Three soft-labeled variants, each its own output set:
- `Group_Pi`  = 4 condition means (Ctl, Mito, PHE, PHE_Mito), proteins Π<0.05 in ≥1 contrast.
- `Group_All` = 4 condition means, all 4,806 proteins.
- `Sample_All`= full 24-sample profiles per protein.

Per variant pipeline:
1. Gene-level matrix via `limma::avereps` on `DAList_imputed_imp4p$data` by gene symbol;
   build the matrix (group means or samples).
2. `ExpressionSet` → `Mfuzz::standardise`.
3. **Cluster-number selection (must be reported):**
   - fuzzifier `m` via `Mfuzz::mestimate` (report value).
   - cluster count `c`: evaluate a range (c = 2..12) by **minimum centroid distance** (Dmin)
     elbow and within-cluster trajectory coherence; pick the elbow, default `c = 9` if no clear elbow.
   - Report `m`, the Dmin-vs-c curve (written to `c_data/<variant>_cluster_selection.csv` and a
     diagnostic PNG), and the chosen `c` with justification in figure subtitle + supplement.
4. `Mfuzz::mfuzz` soft clustering; trajectory line plots per cluster (membership-weighted alpha).
5. **Core members** (membership > 0.5) → `run_ora_deduplicated()` (`fgsea::fora`, per-DB BH,
   odds ratio) against `rat_gene_sets`, universe = all measured genes; Jaccard-0.5 dedup;
   **then drop pathway terms already displayed in 04_Pathway_bars / 05_Enrich_Volcano** so clusters
   surface new biology (dedup-of-prior-figures via a shared "already shown" term list written by 04/05).
6. ORA dot/bar per cluster.
- Outputs per variant: `MAIN_F06_<variant>_traj.{pdf,png}`, `MAIN_F06_<variant>_ora.{pdf,png}`,
  `<variant>_cluster_selection.{csv,png}`, `<variant>_membership.csv`, `<variant>_ora.csv`,
  `F06_supplementary.xlsx`.

## 6. Manuscript prose

Markdown under `docs/` (parallel to YvO IMRAD; same numbers-forward voice, citation pattern):
- **Methods**: cell model & PHE/Mito rationale → DIA-NN acquisition → filtering (5002→4806) +
  cycloess + imp4p → limma + `duplicateCorrelation` replicate block (ρ≈0.02) + Π-score + 5 contrasts
  → fgsea/ORA (rat gene sets) → mfuzz (with selection reporting).
- **Results**: PCA/PERMANOVA → DEP landscape (bars + Venn) → pathway enrichment (bars) →
  volcano rings → mfuzz clusters + ORA; each paragraph effect-size-forward.
- **Discussion**: disease signature, Mito rescue/reversal, cluster programs, limitations
  (n=6, single cell line, single timepoint, in vitro).

## 7. Master build-prompt

`04_Figures_v2/BUILD_PROMPT.md` re-specifies all six figures (inputs, engine calls, style
invariants, outputs) so the suite regenerates from scratch by an agent.

## 8. Defaults chosen (overridable)

- mfuzz default `c = 9` when no clear Dmin elbow; selection always reported.
- 05_Enrich_Volcano covers 4 contrasts (drops secondary MITOvPHE_MITO).
- Venn = 3-set Disease/Transplant/Rescue.

## 9. Out of scope

- No edits to upstream stages (00–03) or to existing `04_Figures/`.
- No new gene-set derivation; reuse `rat_gene_sets.rds`.
- No stitched composite figure (hand-composited downstream, as in YvO Box).
