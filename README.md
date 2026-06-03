# Mito_2026

H9c2 mito-transplantation proteomics pipeline (Nick Kontos dissertation data).

Four cell-culture groups — **Ctl**, **Mito**, **PHE**, **PHE_Mito** — six wells
per group (n = 24 total), paired by Replicate (each replicate index shares
passage / plating day / plate across the four groups). PHE = phenylephrine
(α1-adrenergic hypertrophy proxy for cardiomyocyte remodelling in right-HF).
Mito = isolated mitochondria delivered to the wells. The pipeline asks:

1. Does PHE stress remodel the proteome? (CTLvPHE)
2. Does mito transplant shift it in healthy cells? (CTLvMITO)
3. Does mito transplant rescue the PHE response? (PHEvPHE_MITO)
4. Is rescue additive vs. PHE alone? (Interaction)

Stages 01–03 produce the normalized matrix, imputed matrix, and limma DEP
table. Figures F01–F06 read those frozen outputs and build journal-format
composites.

---

## Requirements

- **R ≥ 4.4** (developed against 4.6)
- ~3 GB free disk for cached gene-set lookups + intermediate artifacts
- About 20 minutes of install time on a fresh machine

### Install R packages (one-time)

Open R in any folder and paste:

```r
install.packages(c(
  "BiocManager", "remotes", "here",
  # tidyverse + I/O
  "dplyr", "tidyr", "tibble", "purrr", "readr", "readxl", "openxlsx",
  "stringr", "scales",
  # plotting
  "ggplot2", "ggrepel", "ggforce", "ggtext", "ggnewscale",
  "patchwork", "cowplot", "gridExtra", "png",
  # stats / modelling
  "missForest", "boot", "pwr", "lme4", "lmerTest", "emmeans",
  "vegan", "WGCNA",
  # gene sets / annotation
  "msigdbr", "babelgene"
))

BiocManager::install(c(
  "limma", "GSVA", "singscore", "fgsea",
  "ComplexHeatmap", "circlize",
  "AnnotationDbi", "GO.db", "org.Rn.eg.db"
))

remotes::install_github("ByrumLab/proteoDA")
remotes::install_github("RischanLab/RRHO2")
```

A note on the stats packages, since it isn't obvious why a few of them are
here when proteoDA is doing the heavy lifting:

- **proteoDA** runs the main DEP workflow — filtering, cycloess normalization,
  the `lmFit` → `eBayes` → contrasts pipeline, and the QC reports.
- **limma** comes in as a proteoDA dependency, but the robustness refits in
  `03_run_robustness.R` call it directly (`duplicateCorrelation`, `lmFit`,
  `eBayes`, `camera`) because they need explicit control over the design
  matrix that proteoDA's wrapper doesn't expose. F06 also uses
  `limma::camera` for the correlation-aware pathway test.
- **lme4 + lmerTest + emmeans** are for the figure-level mixed models — F05
  set scores and F06 mitonuclear / content / stoichiometry — which fit
  `lmer(y ~ PHE * Mito + (1|Replicate))` on per-sample summaries, not on
  protein-by-protein rows. Those LMMs live outside proteoDA's territory.

---

## Get the code

```bash
git clone https://github.com/Dustyn-T-Lewis/Mito_2026.git
cd Mito_2026
```

Then open R in this folder (RStudio: `File > Open File > 01_normalization/a_script/01_normalize.R` and let RStudio set the working directory automatically; VS Code / Positron: open the folder). Every script uses `here::here()` and finds the project root via the `.git/` marker, so cwd is forgiving.

---

## Run the pipeline

Each script is independently re-runnable — re-running Stage 02 doesn't force
a Stage 01 re-run, just reads its frozen `c_data/` artifact. In RStudio, just
open each file and run line-by-line with Ctrl+Enter, or `source()` the whole
file at once.

### Stage 01 — Normalization (~30 s)

```r
source(here::here("01_normalization", "a_script", "01_normalize.R"))
source(here::here("01_normalization", "a_script", "02_generate_reports.R"))
```

Reads `00_input/H9c2_raw.xlsx` (DIA-NN protein groups) and `H9c2_meta.csv`.
Removes contaminants (keratin + FBS serum proteins), filters by missingness
(≥ 4 measurements in ≥ 1 group), flags outliers via 4-method consensus
(missingness, PCA-Mahalanobis, MAD intensity, inter-sample correlation;
≥ 3/4 flags = consensus outlier), and applies cycloess normalization via
proteoDA. Writes:

- `c_data/02_normalized.csv` — protein × sample matrix
- `c_data/03_DAList_normalized.rds` — proteoDA DAList
- `c_data/01_normalization.xlsx` — summary workbook (7 sheets)
- `b_reports/04_diagnostics.pdf` — QC plots

### Stage 02 — Imputation (~2 min)

```r
source(here::here("02_Imputation", "a_script", "01_impute.R"))
source(here::here("02_Imputation", "a_script", "02_generate_reports.R"))
```

3-method MAR/MNAR consensus (k-means on intensity × %missing, global
logistic, left-tail proximity; ≥ 2/3 votes = MNAR) for reporting, then
missForest imputation on *all* proteins. Writes:

- `c_data/01_DAList_imputed.rds` — imputed DAList (figures only)
- `c_data/02_mar_mnar_classification.csv`
- `b_reports/01_missingness_report.pdf`, `02_imputation_report.pdf`

**Note:** Stage 03 DEP uses the NON-imputed cycloess matrix (limma handles
NAs per protein). The imputed matrix feeds figures and WGCNA only.

### Stage 03 — Differential expression (~3 min)

```r
source(here::here("03_DEP", "a_script", "01_run_dep.R"))
source(here::here("03_DEP", "a_script", "02_generate_reports.R"))
source(here::here("03_DEP", "a_script", "03_run_robustness.R"))
```

limma fit with `~ 0 + group + (1|Replicate)` — proteoDA picks up the block
and routes through `duplicateCorrelation` (Smyth, Michaud & Scott 2005).
Five contrasts: CTLvPHE, CTLvMITO, PHEvPHE_MITO, Interaction, MITOvPHE_MITO.
Significance via BH-FDR and Pi-score (`P.Value^|logFC|`, Xiao 2014). Writes:

- `c_data/03_combined_results.csv` — wide-format DEP table
- `c_data/03_DEP_results.xlsx` — per-contrast sheets, DA_summary,
  bootstrap CI, power analysis, outlier / reinjection / imputation /
  robust-eBayes sensitivity sheets

### Build the figures

Each `05_Figures/F0X_*/` directory builds one composite. Order doesn't
matter — they all read frozen Stage 03 outputs.

```r
# F01 — QC + proteome overview (PCA, DEPs, enrichment, UpSet, rank)
source(here::here("05_Figures", "F01_QC_overview", "a_script", "01_main_panels.R"))

# F02 — All-contrast volcano-in-ring (3 enrichment lenses)
source(here::here("05_Figures", "F02_volcano",     "a_script", "01_main_panels.R"))

# F03 — Concordance: Disease vs Intervention
source(here::here("05_Figures", "F03_concordance", "a_script", "01_main_panels.R"))

# F04 — Rescue: Disease vs Rescue
source(here::here("05_Figures", "F04_rescue",      "a_script", "01_main_panels.R"))

# F05 — Set-score modules (GSVA + singscore + 2x2 LMM)  (~3 min)
source(here::here("05_Figures", "F05_modules",     "a_script", "01_set_scores.R"))
source(here::here("05_Figures", "F05_modules",     "a_script", "02_main_panels.R"))

# F06 — Mitochondrial complexes, content, balance, camera
source(here::here("05_Figures", "F06_complex_mito", "a_script", "01_analysis.R"))
source(here::here("05_Figures", "F06_complex_mito", "a_script", "02_lmm_dynamics.R"))
source(here::here("05_Figures", "F06_complex_mito", "a_script", "03_figure.R"))
```

Each composite lands in `05_Figures/F0X_*/b_reports/main/{pdf,png}/`.
PDFs print at 178 mm width (Nature double-column).

### Optional supplements

```r
# F01 QC supplement (normalization + imputation diagnostics)
source(here::here("05_Figures", "F01_QC_overview", "a_script", "02_supp_panels.R"))

# F05 WGCNA supplement — build the network first (~3 min)
source(here::here("05_Figures", "shared",      "build_wgcna_network.R"))
source(here::here("05_Figures", "F05_modules", "a_script", "91_wgcna_supp.R"))

# F06 covariation supplement (ProteomeHD-flavoured)
source(here::here("05_Figures", "F06_complex_mito", "a_script", "04_covariation.R"))
```

### Optional: rebuild gene-set caches

The figures consume frozen fGSEA caches in `04_Figures/shared/`. To rebuild
from current MSigDB / CORUM / MitoCarta releases:

```r
source(here::here("04_Figures", "shared", "build_fgsea_cache.R"))
source(here::here("04_Figures", "shared", "build_fgsea_corum_h9c2.R"))
source(here::here("04_Figures", "shared", "build_fgsea_goslim_h9c2.R"))
source(here::here("04_Figures", "shared", "build_localization_lookup.R"))
```

---

## Layout

```
Mito_2026/
├── 00_input/                          raw data + design constants
│   ├── H9c2_raw.xlsx                  DIA-NN protein-group output
│   ├── H9c2_meta.csv                  sample metadata (Col_ID, Group, Replicate, ...)
│   ├── H9c2_collaborator_DAP_lists.xlsx
│   └── h9c2_design.R                  groups, contrasts, palettes, thresholds, helpers
│
├── 01_normalization/
│   ├── a_script/                      pipeline scripts (committed)
│   ├── b_reports/                     PDFs (gitignored, regenerated)
│   └── c_data/                        .rds/.csv outputs (gitignored, regenerated)
│
├── 02_Imputation/                     same a_script/b_reports/c_data layout
├── 03_DEP/                            same
│
├── 04_Figures/shared/                 cross-figure utilities + frozen caches
│   ├── style.R                        palettes, themes, sizing
│   ├── pathway_utils.R                fGSEA + Jaccard dedup
│   ├── mitocarta_utils.R              MitoCarta lookups
│   ├── volcano_ring.R                 F02 ring engine
│   ├── *.rds, *.csv                   frozen gene-set caches (rat-specific)
│   └── comparison_panels/             reusable panel scripts (B/C/D/E)
│
├── 05_Figures/                        one folder per composite (F01..F06)
│   ├── F01_QC_overview/
│   ├── F02_volcano/
│   ├── F03_concordance/
│   ├── F04_rescue/
│   ├── F05_modules/
│   ├── F06_complex_mito/
│   └── shared/                        config.R, build_wgcna_network.R
│
├── .gitignore
└── README.md                          (you are here)
```

`b_reports/` and `c_data/` directories ship empty (with a `.gitkeep`
placeholder); the scripts populate them on first run.

---

## Notes

- **Run from the project root.** `here::here()` finds the `.git/` marker, so
  the project root resolves correctly even if your R session was started
  somewhere else. Just don't `setwd()` away mid-session.

- **Stage order matters.** Stage 02 reads Stage 01 outputs; Stage 03 reads
  Stage 01 outputs directly (DEP uses the non-imputed matrix); figures read
  Stage 03 outputs (plus the Stage 02 imputed matrix where set scores or
  per-sample heatmaps need a dense matrix).

- **Replicate is a real block.** Wells sharing a Replicate index also share
  passage / plate / day. The design formula is
  `~ 0 + group + (1|Replicate)`, and proteoDA detects the random-effect
  term and routes through `duplicateCorrelation` automatically. F05/F06
  LMMs use the same `(1|Replicate)` block.

- **N = 24 — the orthogonal Interaction is underpowered.** Composite
  footers carry this caveat. Interaction is reported descriptively only.

- **The eBayes "> 20% DA" warning is intentional.** CTLvMITO and
  PHEvPHE_MITO legitimately move > 20% of the proteome because mito
  transplant delivers organelles. The robust-eBayes sensitivity sheet in
  `03_DEP_results.xlsx` confirms the standard fit is stable (Spearman ρ
  ≈ 0.99).

- **Reinjection imbalance is flagged.** A few wells were re-injected and
  the imbalance across groups is documented in `03_run_robustness.R`. The
  sensitivity refit shows core contrasts robust (Spearman ρ ≈ 0.95–0.99);
  MITOvPHE_MITO is weaker (≈ 0.82). Review before publication.

- **Cycloess normalization is the canonical norm.** Stage 01 caches a
  pre-outlier cycloess matrix into `00_report_intermediates.rds` so Stage
  03's outlier-sensitivity refit doesn't repeat the normalization.

- **macOS / cairo warnings are cosmetic.** If you see
  `failed to load cairo DLL` on macOS, the figure scripts fall back to
  Quartz and outputs are identical.

---

## Background

- PHE-induced H9c2 hypertrophy: Hahn 2014 (Cell Signal); Jeong 2009
  (Exp Mol Med).
- Mitochondrial transplant cardiac rescue: Masuzawa 2013 (AJP-Heart);
  Doulamis 2024 (Sci Rep).
- Pi-score significance: Xiao 2014 (Bioinformatics).
- Block design `(1|Replicate)`: Smyth, Michaud & Scott 2005 (Bioinformatics).
- Method references live in `03_DEP/a_script/01_run_dep.R` head comment and
  `05_Figures/F02_volcano/c_data/F02_citations.csv`.
