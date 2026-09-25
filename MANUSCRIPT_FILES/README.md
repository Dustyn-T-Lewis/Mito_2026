# Manuscript supplementary files

S1-S4 Figure and S1-S3 Table, named as the manuscript's Supplemental Appendix cites them.
Each S Figure is one multipage PDF with one panel per page, each page stamped with its
label. Each S Table is the supplementary workbook of the matching figure; its `Overview`
sheet lists the sheets.

| File | Pages / content | Source panels or workbook |
| --- | --- | --- |
| `S1_Figure.pdf` | A PCA by group with PERMANOVA; B PCA by transplant; C fgsea on PC1 loadings; D ORA of significant proteins; E MitoCarta compartment shares | `04_Figures/F01_Proteome_Overview/b_reports/supp/` |
| `S2_Figure.pdf` | A top pathways per contrast; B disease-specificity of the transplant response; C disease-selected sets read across contrasts | `04_Figures/F02_Enrich_Volcanoes/b_reports/supp/` |
| `S3_Figure.pdf` | A network construction; B module preservation; C hub proteins per module | `04_Figures/F03_WGCNA/b_reports/supp/` |
| `S4_Figure.pdf` | Contrast-coupling permutation nulls | `04_Figures/F03_WGCNA/b_reports/supp/SUPP_F03_orthogonal_axes.png` |
| `S1_Table.xlsx` | Proteome overview, DE tables for all six contrasts | `04_Figures/F01_Proteome_Overview/c_data/F01_supplementary.xlsx` |
| `S2_Table.xlsx` | Every tested pathway per contrast with deduplication status | `04_Figures/F02_Enrich_Volcanoes/c_data/F02_supplementary.xlsx` |
| `S3_Table.xlsx` | Module landscape, eigengene statistics, preservation, hubs | `04_Figures/F03_WGCNA/c_data/F03_supplementary.xlsx` |

Rebuild after the figure scripts, from the project root:

```sh
Rscript MANUSCRIPT_FILES/build_manuscript_files.R
```

The panels themselves come from the figure scripts: `F01`/`F02` `a_script/supp/`,
`F03_WGCNA/a_script/02_clustering.R` (construction, preservation, hubs) and
`F03_WGCNA/a_script/supp/orthogonal_axes.R`.
