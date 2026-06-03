# Mito_2026

H9c2 mito-transplantation proteomics pipeline (Nick Kontos dissertation data).
Four groups (Ctl, Mito, PHE, PHE_Mito), n = 6/group, paired by Replicate
(passage / plating day / plate). Stages 01–03 run normalization, imputation,
and DEP (limma + duplicateCorrelation on the Replicate block). Figures F01–F06
under `05_Figures/` build on the frozen Stage 01–03 outputs in each stage's
`c_data/`.

Pipeline outputs (`b_reports/`, `c_data/`) are gitignored — regenerate locally
by running each stage's `a_script/*.R` in order.

