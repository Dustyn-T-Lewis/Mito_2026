# shared/ — shared reference data

Builds the canonical databases every figure enriches against, once, and caches them. Data and its build scripts; no figures.

- `a_script/01_fetch_rat_gene_sets.R` → `c_data/rat_gene_sets.rds` (Hallmark, Reactome, KEGG, GO-Slim, MitoCarta; within-DB Jaccard dedup)
- `a_script/03_fetch_go_sets.R` → `c_data/rat_go_bpccmf_sets.rds` (GO:BP/CC/MF for module ORA)
- `a_script/02_build_fgsea_cache.R` → `c_data/fgsea_tstat_all_h9c2.csv` (moderated-t fGSEA per contrast × database)
- `a_script/00_gene_set_helpers.R` — builders the above scripts share
- `mitocarta3_rat.rds` — tracked MitoCarta3.0 source; `PROVENANCE.md` documents it

Reusable **code** (loaders, enrichment engine, theme) lives in `functions/`, not here.
