# functions/ — shared figure code

The reusable engine every figure `source()`s. Code only; no data, no figures.

| File | Provides |
|---|---|
| `F01-F03_style_palettes_theme.R` | `FIG_THEME`, palettes, contrast labels, `H9C2_PI_THRESH` |
| `F01-F03_data_paths_and_loaders.R` | canonical input paths (`P05`), `load_combined_wide`, `load_fgsea_cache`, `contrast_brief` |
| `F01-F03_pathway_enrichment_dedup_ora.R` | `significant_pathways`, EnrichmentMap `deduplicate_enrichment`, `run_fora_by_db` |
| `F01-F03_supplementary_workbook.R` | `build_workbook` / `append_workbook` |
| `F01-F03_composite_layout.R` | `add_tag`, `save_composite` |
| `F03-F04_comparison_panels.R` | contrast-pair panels: RRHO2, NES scatter, fry barcode, quadrant ORA |

Reference **data** (the gene-set collections and the fGSEA cache) lives in `shared/`, not here.

Note: the `F01-F03_` tags predate the current figure numbering (F01, F02, F03_pilot, F04_WGCNA) and are now understated — these helpers are shared across all figures. Retag when convenient.
