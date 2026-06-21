# Compatibility stub: pathway enrichment utilities split across three siblings.
# 03a — Jaccard/overlap dedup engine + Fisher 2x2 odds ratio
# 03b — gene-set collection + fgsea + fora runners (depends on 03a)
# 03c — display-name dictionary + functional category classifier
source(here::here("04_Figures_v2", "functions", "03a_dedup_engine.R"))
source(here::here("04_Figures_v2", "functions", "03b_enrichment_runners.R"))
source(here::here("04_Figures_v2", "functions", "03c_pathway_display.R"))
