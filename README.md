# Mito Proteomics Analysis

Label-free DIA-MS proteomics of H9c2 rat cardiomyoblasts under a 2x2 factorial
design crossing phenylephrine (`PHE`) hypertrophic stress with mitochondrial
transplantation (`Mito`), n = 6 per group:

- `Ctl`: untreated control
- `Mito`: mitochondrial transplant, no stress
- `PHE`: phenylephrine stress
- `PHE_Mito`: phenylephrine stress with mitochondrial transplant

The `Replicate` field indexes the six within-group replicates and is carried into
the differential-abundance fit as a `duplicateCorrelation` block.

## Design and Contrasts

DEP fits the means model `~ 0 + group` (one mean per group) with
`duplicateCorrelation` blocking on `Replicate`, and reports six contrasts. Each
is a linear combination of the four group means, so all are estimable.

- `CTLvPHE = PHE - Ctl` — disease: PHE stress remodelling
- `CTLvMITO = Mito - Ctl` — transplant: mitochondria into healthy cells
- `PHEvPHE_MITO = PHE_Mito - PHE` — rescue: transplant under stress
- `CTLvPHE_MITO = PHE_Mito - Ctl` — recovery: residual disease vs control
- `Interaction = (PHE_Mito - PHE) - (Mito - Ctl)` — rescue minus transplant
- `MITOvPHE_MITO = PHE_Mito - Mito` — secondary: PHE effect in transplanted cells

Discovery uses the Pi-score (`Pi = P.Value^|log2FC|`, Xiao et al. 2014) at
`Pi < 0.05`. Pi is an effect-weighted cutoff on the **raw** p-value — a discovery
cutoff, not an FDR. Benjamini-Hochberg FDR is reported alongside as a separate,
FDR-controlled column; the two are never conflated. The Disease (`CTLvPHE`) and
Rescue (`PHEvPHE_MITO`) contrasts share the PHE group with opposite sign, so their
correlation is coupled — rescue is tested with the `Interaction` term or a direct
per-sample group comparison, never the coupled contrast correlation. The
`duplicateCorrelation` consensus is low (rho ~ 0.02) and is retained in the fit
as a blocked generalized-least-squares term rather than dropped.

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | Raw DIA-NN intensity matrix, sample sheet, contaminant FASTA, HPA secreted-to-blood list |
| `01` | `01_Filtering/` | Contaminant and HPA secreted-to-blood removal, missingness filter -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization of the filtered matrix; `imputation/` holds three arms (`imp4p`, MsCoreUtils hybrid, `missForest`), each writing a method-tagged `DAList_imputed_<method>.rds`. `missForest` is the arm the figures and WGCNA read |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, six factorial contrasts, Pi-score. `b_imputed/`: robustness DEP on the imputed matrices |
| `04` | `04_Figures/` | F01 proteome overview; F02 enrichment ring-volcanoes (reads the shared fgsea cache); F03 WGCNA modules |

Stages `00`–`03` are the canonical pipeline (filter → normalize → DE) and hold
only that. Anything beyond canonical characterisation — WGCNA, single-sample
scoring, reversal/return — lives with its figure in `04_Figures/`, never in a
pipeline stage.

### `04_Figures/F03_WGCNA` — the clustering analysis

WGCNA is the last figure and self-contained under `04_Figures/F03_WGCNA/a_script`.
`00_build_wgcna.R` builds a signed WGCNA network with biweight midcorrelation
(bicor, `maxPOutliers = 0.05`) on the missForest-imputed matrix, extracts module
eigengenes and hub proteins, and tests module preservation (Zsummary) of the
control-reference modules in each treatment group. `01_module_stats.R` runs
module-level statistics on the same `~ 0 + group` design and `Replicate` block as
the protein DEP: eigengene limma per contrast, an omnibus moderated-F, and
`fry`/`camera` set tests. At n = 6 per group these are exploratory; the protein
DEP in `a_non_imputed` is the confirmatory result.

## Setup

Install every package the pipeline uses, once after cloning:

```sh
Rscript setup.R
```

## Canonical Run Order

Run the stage scripts from the project root, in order. Each is self-contained:
it loads its own packages, creates its output dirs, and reads the previous
stage's `c_data`.

### Core stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R
Rscript 02_Normalization/imputation/a_script/c_missforest.R   # arm read by figures + WGCNA

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R            # primary DEP
Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R        # robustness DEP
```

The other two imputation arms (`a_imp4p.R`, `b_mscoreutils.R`) are exploratory
alternatives; run them only for the imputed-DEP concordance comparison.

### Figures

```sh
Rscript 04_Figures/shared/a_script/01_fetch_rat_gene_sets.R   # once; builds gene-set cache
Rscript 04_Figures/shared/a_script/02_build_fgsea_cache.R     # once; fgsea cache
Rscript 04_Figures/shared/a_script/03_fetch_go_sets.R         # once; GO sets

# each figure is self-contained; run its composite from the project root
Rscript 04_Figures/F01_Proteome_Overview/a_script/01_proteome_overview.R
Rscript 04_Figures/F02_Enrich_Volcanoes/a_script/01_enrich_volcanoes.R
Rscript 04_Figures/F03_WGCNA/a_script/00_build_wgcna.R
Rscript 04_Figures/F03_WGCNA/a_script/01_module_stats.R
Rscript 04_Figures/F03_WGCNA/a_script/02_clustering.R
```

Each figure's composite sources its panels, runs the shared fit or load, and
writes the figure + workbook. F01, F02, and F03 end with an optional local Box
mirror block (see Conventions).

## Repository Conventions

Every stage and figure uses the same `a_/b_/c_` triad:

- `a_script/`: the scripts. One panel builder per file; a `00_`/`01_` composite
  sources the panels, runs the shared fit or load, and writes the figure plus its
  supplementary workbook. Tests and supplements sit in named subfolders.
- `b_reports/`: generated renders (PDF/PNG) and QC reports.
- `c_data/`: the tables downstream steps read.

Shared code and shared data are kept apart:

- `04_Figures/functions/`: every reusable helper the figure scripts source. The
  filename states its scope. `shared_*` is used across figures —
  `shared_theme_palettes.R` (palettes, theme, thresholds), `shared_data_loaders.R`
  (input paths and loaders), `shared_enrichment_ora.R` (fgsea/ORA and pathway
  dedup), `shared_gene_set_helpers.R` (fgsea cache and GO-Slim set builders),
  `shared_nes_bars.R` (the NES bar panel), `shared_workbook.R`, and
  `shared_composite_layout.R`. A `f0N_` prefix marks a helper built for one
  figure — `f01_pca_stats.R`, `f01_mitocarta_lens.R`.
- `04_Figures/shared/`: the shared gene-set data and fgsea cache. Built once by
  `shared/a_script/01`–`03` and read by the enrichment figures. Only the build
  scripts live here; the helpers they call sit in `functions/`.

Every path resolves from the project root through `here::here()`, so any script
runs standalone regardless of the working directory it is launched from. Box sync
is optional and local: the F01, F02, and F03 scripts end with a `mirror_to_box`
block that copies their outputs to the author's Box folder and no-ops silently
when that folder is not mounted. Skip it by not running the trailing block.

## Gene Sets and Provenance

The enrichment figures (F01, F02) and the WGCNA module ORA (F03) share one rat
gene-set collection and one fgsea cache, built once under `04_Figures/shared/`:

- **Hallmark, KEGG, Reactome** — `msigdbr(species = "Rattus norvegicus")`. KEGG
  uses `CP:KEGG_LEGACY`, the older licensed pathway set, for continuity; both KEGG
  curations ortholog-map to rat.
- **GO Slim** — the 62 Biological Process terms of the GO Consortium Generic
  GO-Slim, expanded over `GO.db::GOBPOFFSPRING`, mapped to symbols via
  `org.Rn.eg.db`, and size-filtered to 10–500.
- **MitoCarta** — the Broad MitoCarta3.0 rat inventory, tracked as
  `shared/mitocarta3_rat.rds` because the upstream `.xls` is not committed.

Each database is collapsed within itself by Jaccard ≥ 0.5 (larger set kept);
cross-database redundancy is resolved at runtime by `shared_enrichment_ora.R`. The
`msigdbr`, `GO.db`, and `org.Rn.eg.db` versions are stamped into the collection's
`provenance` attribute. With `msigdbr 26.1.0` the rebuilt cache is byte-for-byte
identical to the committed one; a newer release may shift set membership, so diff
against the prior cache before adopting a rebuild.

## Reproducibility Rules

- Path resolution uses `here::here()` from the project root
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed normalized matrix; imputation feeds only
  figures and the WGCNA modules
- the fit blocks on the `Replicate` metadata via duplicateCorrelation (low consensus rho)
- species is *Rattus norvegicus*; gene sets are the rat collections cached under
  `04_Figures/shared/`
- packages install into the system library via `Rscript setup.R`; each script
  loads what it needs with `pacman::p_load()`

## enrichVolcano

The figures depend on `enrichVolcano` (`Dustyn-T-Lewis/enrichVolcano`), installed
by `setup.R`. F02 passes `x_scale`/`y_scale`, so it needs `d6a229e` or later; to
adopt a newer dev build, `remotes::install_github("Dustyn-T-Lewis/enrichVolcano@<sha>")`
and re-run F02.
