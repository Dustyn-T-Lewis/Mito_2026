# Mito Proteomics Analysis

Label-free DIA-MS proteomics of H9c2 rat cardiomyoblasts under a 2x2 factorial
design crossing phenylephrine (`PHE`) hypertrophic stress with mitochondrial
transplantation (`Mito`), n = 6 per group:

- `Ctl`: untreated control
- `Mito`: mitochondrial transplant, no stress
- `PHE`: phenylephrine stress
- `PHE_Mito`: phenylephrine stress with mitochondrial transplant

The design carries a paired `Replicate` block (plate/day/passage) across the four
groups, handled as a repeated-measures term in the differential-abundance fit.

## Design and Contrasts

DEP fits the means model `~ 0 + group` (one mean per group) with
`duplicateCorrelation` blocking on `Replicate`, and reports five contrasts. Each
is a linear combination of the four group means, so all are estimable.

- `CTLvPHE = PHE - Ctl` — disease: PHE stress remodelling
- `CTLvMITO = Mito - Ctl` — transplant: mitochondria into healthy cells
- `PHEvPHE_MITO = PHE_Mito - PHE` — rescue: transplant under stress
- `Interaction = (PHE_Mito - PHE) - (Mito - Ctl)` — rescue minus transplant
- `MITOvPHE_MITO = PHE_Mito - Mito` — secondary: PHE effect in transplanted cells

Significance is the Pi-score (`Pi = P.Value^|log2FC|`, Xiao et al. 2014) at
`Pi < 0.05`, with Benjamini-Hochberg FDR < 0.10 as a secondary criterion. The
`duplicateCorrelation` consensus is low (rho ~ 0.02) and is retained in the fit
as a blocked generalized-least-squares term rather than dropped.

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | Raw DIA-NN intensity matrix, sample sheet, contaminant FASTA, HPA secreted-to-blood list |
| `01` | `01_Filtering/` | Contaminant and HPA secreted-to-blood removal, missingness filter -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization of the filtered matrix; `imputation/` holds three arms (`imp4p`, MsCoreUtils hybrid, `missForest`), each writing a method-tagged `DAList_imputed_<method>.rds`. `missForest` is the arm the figures and WGCNA read |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, five contrasts, Pi-score. `b_imputed/`: robustness DEP on the imputed matrices. `c_modules/`: WGCNA co-abundance network and module-level statistics (see below) |
| `04` | `04_Figures/` | F01 proteome overview; F02 enrichment ring-volcanoes (builds the shared fgsea cache); F03 clustering, which renders the `c_modules` WGCNA |

### `03_DEP/c_modules` — the clustering analysis

`c_modules` is the canonical clustering result, rendered by `F03_Clustering`.
`00_build_wgcna.R` builds a signed WGCNA network with biweight midcorrelation
(bicor, `maxPOutliers = 0.05`) on the missForest-imputed matrix, extracts module
eigengenes and hub proteins, and tests module preservation (Zsummary) of the
control-reference modules in each treatment group. `01_module_stats.R` runs
module-level statistics on the same `~ 0 + group` design and `Replicate` block as
the protein DEP: eigengene limma per contrast, an omnibus moderated-F, and
`fry`/`camera` set tests. At n = 6 per group these are exploratory; the protein
DEP in `a_non_imputed` is the confirmatory result.

## Canonical Run Order

### Core stages

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R
Rscript 02_Normalization/imputation/a_script/c_missforest.R   # arm read by figures + WGCNA

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R            # primary DEP
Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R        # robustness DEP
Rscript 03_DEP/c_modules/a_script/00_build_wgcna.R            # WGCNA network
Rscript 03_DEP/c_modules/a_script/01_module_stats.R          # module statistics
```

The other two imputation arms (`a_imp4p.R`, `b_mscoreutils.R`) are exploratory
alternatives; run them only for the imputed-DEP concordance comparison.

### Figures

```sh
Rscript 04_Figures/shared/a_script/01_fetch_rat_gene_sets.R   # once; builds gene-set cache
Rscript 04_Figures/shared/a_script/03_fetch_go_sets.R         # once; GO sets
Rscript 04_Figures/run_all.R                                  # fgsea cache + F01, F02, F03
```

Each figure's main script writes its workbook; the supplement scripts append to
it, so `run_all.R` runs each main before its supplement. Each script ends with an
optional local Box mirror block (see Conventions).

## Repository Conventions

- `a_script/`: scripts and narrative notebooks
- `b_reports/`: generated figure renders (PDF/PNG) and QC reports
- `c_data/`: stage outputs read by downstream steps
- `functions/` and `shared/`: reusable helpers that scripts source
- Figure-script style is documented in `04_Figures/CONVENTIONS.md`
- Box sync is optional and local: each figure script ends with a `mirror_to_box`
  block that copies its outputs to the author's Box folder and no-ops silently
  when that folder is not mounted. Skip it by not running the trailing block.

## Reproducibility Rules

- Path resolution uses `here::here()` from the project root
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed normalized matrix; imputation feeds only
  figures and the WGCNA modules
- repeated-measures blocking uses the paired `Replicate` metadata
- species is *Rattus norvegicus*; gene sets are the rat collections cached under
  `04_Figures/shared/`
