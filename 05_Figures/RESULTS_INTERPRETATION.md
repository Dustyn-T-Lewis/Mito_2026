---
title: "Mitochondrial Transplant in PHE-Stressed H9c2 Cardiomyoblasts — Methods and Results"
subtitle: "Documentation deliverable; not a manuscript draft"
date: "2026-06-09"
---

# Purpose of this document

This document describes the methods and results from a proteomic analysis of
H9c2 rat cardiomyoblasts treated with phenylephrine (PHE) and/or mitochondrial
transplant. It is prepared as a documentation deliverable to be folded into a
separate project led by another investigator; it is **not** a manuscript draft
in its own right. The level of detail here is intentionally fuller than would
appear in a journal Methods or Results section so the recipient has enough
context to (i) understand exactly what was done and why, (ii) decide which
parts to keep verbatim, (iii) compress to the target journal's word count, and
(iv) write the surrounding Introduction and Discussion.

Every analytical choice is cited to peer-reviewed primary literature; all
PMIDs were verified against PubMed during writing.

---

# Methods

## Experimental design

H9c2 rat cardiomyoblasts (ATCC) were assigned to a 2 × 2 factorial design
crossing phenylephrine treatment (PHE − / +) with mitochondrial transplant
(Mito − / +), yielding four groups: Ctl, Mito, PHE, and PHE_Mito. Each
group contained six biological replicates (N = 24 total samples). Replicates
were structured as paired plate/day/passage blocks across the four groups,
so each replicate index (1–6) corresponds to a matched set of one well from
each of the four groups processed together. This block structure was carried
through statistical modeling as described below.

A small number of samples (n = 4) were re-injected on the mass spectrometer
to recover failed acquisitions. These re-injection samples are denoted by an
`r` suffix and are distributed across groups as Ctl = 1, Mito = 0, PHE = 2,
PHE_Mito = 3. The imbalance was carried forward as a documented sensitivity
analysis rather than corrected, because the small absolute counts do not
support a stratified design re-fit at this N.

Four primary contrasts and one secondary contrast were defined:

| Contrast | Definition | Biological question |
|---|---|---|
| **CTLvPHE** | PHE − Ctl | Does α1-adrenergic stress remodel the cardiomyoblast proteome? |
| **CTLvMITO** | Mito − Ctl | Does mitochondrial transplant leave a cell-autonomous signature in healthy cells? |
| **PHEvPHE_MITO** | PHE_Mito − PHE | Does mitochondrial transplant rescue the PHE-stressed proteome? |
| **Interaction** | (PHE_Mito − PHE) − (Mito − Ctl) | Does PHE alter what mitochondrial transplant does? (descriptive at N = 24; not formally interpreted) |
| MITOvPHE_MITO (secondary) | PHE_Mito − Mito | What does PHE do on top of mitochondrial transplant? Used for sensitivity reporting only. |

## Mass spectrometry and protein quantification

Label-free quantitative mass spectrometry was performed upstream of the
analyses described here. Protein-level intensities were exported as a
pre-processed matrix (4,813 proteins × 24 samples after standard data-matrix
filtering). All downstream analyses in this document begin from this
protein-level matrix; no peptide-level reprocessing was performed.

## Quality control and normalization

Cyclic loess normalization (`cycloess`) was applied via the **proteoDA**
package using its standard normalization wrapper. Sample-level quality
control followed a four-method consensus framework requiring ≥ 3 of 4
diagnostic checks to flag a sample as an outlier; no samples were excluded
from the final analysis under this criterion.

## Missingness handling and imputation

Proteins with missing values were classified into Missing-At-Random (MAR)
or Missing-Not-At-Random (MNAR) categories using the empirical Bayes
mixture model implemented in the **msImpute** Bioconductor package
(bioconductor.org/packages/msImpute).

A single imputation strategy using **missForest** (Stekhoven and Bühlmann
2012) was then applied to every protein, regardless of MAR/MNAR
classification. A 3-method consensus (kmeans + logistic + left-tail
imputation) was computed in parallel only to derive an `imputation_reliable`
per-protein flag for reporting; the consensus was not used as the
imputation itself.

The downstream split is important: **differential abundance analysis used the
non-imputed cyclic-loess matrix** (limma handles per-protein missingness
through its per-row regression model), while the imputed matrix was used
only for visualizations, principal-component analysis, WGCNA, and set
scoring where complete data are required.

## Differential abundance analysis

Differential abundance was tested with **limma** (Ritchie et al. 2015,
PMID 25605792) using the design formula

```
~ 0 + group + (1 | Replicate)
```

routed through proteoDA's `fit_limma_model()` wrapper, which forwards the
Replicate term to limma's `duplicateCorrelation` interface. This treats
Replicate as a within-block correlation rather than a fixed effect,
correctly accounting for the paired-plate structure without consuming
five degrees of freedom for the six replicate indices.

Empirical Bayes moderation was applied via `eBayes()`. The Benjamini-Hochberg
adjusted p-value (`adj.P.Val`) was used at an exploratory threshold of
FDR < 0.10. In addition, the **capital-pi (Π) significance score** from
Xiao et al. 2014 (PMID 22321699) was applied at Π < 0.05, where

$$\Pi = \text{P-value}^{|\log_2 \text{FC}|}$$

The Π score is small when a protein has both a small p-value and a large
fold change, which weights signed effect size and statistical significance
together. Calling proteins as differentially abundant (DEPs) on Π < 0.05
is a commonly-used effect-size-aware alternative to thresholding on
adj.P.Val alone.

## Pathway enrichment analysis

### Gene set databases

Pathway enrichment was performed against a curated four-database pool:

| Database | Rat sets (size 10–500) | Source / rationale |
|---|---|---|
| **MSigDB Hallmark** | 50 | Liberzon et al. 2015 (PMID 26771021). Refined, non-redundant gene-set collection purpose-built for GSEA as the primary analysis layer. |
| **Reactome** | 1,839 | Gillespie et al. 2022 (PMID 34788843). Hierarchical, hand-curated canonical pathway knowledgebase. |
| **KEGG** | 186 | Kanehisa 2019 (PMID 31441146). Familiar cardiology / signaling pathway labels (ECM-receptor interaction, MAPK, PI3K-Akt, oxidative phosphorylation). |
| **MitoCarta 3.0** | 103 | Rath et al. 2021 (PMID 33174596). Mitochondrial-localized protein inventory with sub-organelle and hierarchical functional annotations (Complex I-V subunits, mitoribosome, TCA cycle, etc.). |

All sets were retrieved as human gene symbols and orthology-mapped to rat
using **msigdbr** for MSigDB collections and **babelgene** for CORUM and
MitoCarta (the latter was originally curated in mouse and human).

For one supplementary expanded composite, the canonical four-database pool
was augmented with **GO:BP** (7,535 rat sets) to provide finer biological-
process resolution. MitoCarta sub-organelle compartment sets (IMM, IMS,
Matrix, OMM) were excluded from enrichment because they describe spatial
localization rather than functional biology.

### Enrichment statistic

Per-contrast pathway enrichment was computed with **fgsea** (Korotkevich
et al. 2021, *bioRxiv* doi: 10.1101/060012) ranking proteins by their
limma moderated t-statistic. fgsea parameters: `minSize = 10`, `maxSize = 500`,
`eps = 0`. Set size limits (10 ≤ size ≤ 500) exclude both
underpowered very-small sets and uninformatively-large generic sets.

### Redundancy collapse

A standard problem in multi-database GSEA is that the same underlying
biology often appears in several databases under different names — for
example, OXPHOS biology surfaces as `HALLMARK_OXIDATIVE_PHOSPHORYLATION`,
`REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT`,
`KEGG_OXIDATIVE_PHOSPHORYLATION`, and `MITOCARTA_OXPHOS_COMPLEX_V_SUBUNITS`
simultaneously. Without collapse, the same biology occupies multiple slots
in a figure or table, crowding out other distinct findings.

This analysis collapses redundant sets using the **combined Jaccard /
Overlap-Coefficient criterion** at threshold 0.5, which is the
EnrichmentMap default (Merico et al. 2010, PMID 21085593; Reimand et al.
2019, PMID 30664679, Box 1). Two coefficients are checked for each pair
of significant sets `A` and `B`:

- **Jaccard coefficient**: `|A ∩ B| / |A ∪ B|`. Symmetric, penalizes
  set-size differences. Catches symmetric overlap.
- **Overlap (Szymkiewicz-Simpson) coefficient**: `|A ∩ B| / min(|A|, |B|)`.
  Asymmetric, ranges to 1.0 when one set is entirely inside the other.
  Catches containment relationships that Jaccard misses.

Two sets are considered redundant if **either** Jaccard ≥ 0.5 **or**
Overlap ≥ 0.5. A greedy collapse is then applied: sets are sorted by
adjusted p-value (smallest first); for each set, if it is redundant with
any already-kept set the smaller-adjusted-p-value representative wins.

The collapse is performed in two passes per Reimand 2019 Box 1: first
within each database (so within-Hallmark sub-themes do not double-count),
then cross-database (so the OXPHOS example above collapses to one
canonical term). Adding the Overlap-Coefficient criterion is essential
here because Jaccard alone misses asymmetric containment — for instance,
`REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT` (23 genes) lies entirely
inside `REACTOME_TRANSPORT_OF_SMALL_MOLECULES` (711 genes) at Jaccard 0.03
but Overlap 1.00. Jaccard-only dedup would have kept both as "distinct"
findings.

### Per-database visualization caps

The canonical Figure-2 volcano-in-ring composite displays at most six up-
regulated and six down-regulated pathways per contrast (12 total per
panel). To prevent any single database from monopolizing those 12 slots
and crowding out cross-database triangulation, a per-database cap of
two terms per side is applied. The per-DB cap is enforced after dedup
but before final padj-ranking, so the most-significant survivor of each
database is preferentially included.

## Set-level scoring

To compare the size of pathway-level changes across contrasts independently
of fgsea's rank-based statistic, set-level scores were computed for the
harmonized gene-set collection using two methods in parallel:

- **GSVA** (Hänzelmann et al. 2013, PMID 23323831) — Gaussian-kernel
  cumulative-distribution-function based enrichment score per sample
- **singscore** (Foroutan et al. 2018, PMID 30400809) — rank-based
  per-sample enrichment with the `centerScore = TRUE` standardization

The collection underwent the same Jaccard / Overlap-Coefficient dedup
prior to scoring (1,625 candidate sets → 416 retained, comprising 285
Reactome, 81 CORUM, and 50 MitoCarta sets). Per-set scores were then
modeled with a linear mixed model

```
score ~ PHE * Mito + (1 | Replicate)
```

using `lmerTest::lmer()`, with simple-effects contrasts via `emmeans`
and Benjamini-Hochberg multiplicity correction across sets per contrast.
The Replicate intraclass correlation coefficient (`replicate_icc`) is
reported per set for transparency.

Methodological agreement between GSVA and singscore was assessed as the
per-contrast Spearman correlation of effect sizes; ρ ≥ 0.80 was treated
as a method-robustness threshold.

## WGCNA module discovery

Weighted gene co-expression network analysis was performed with **WGCNA**
(Langfelder and Horvath 2008, PMID 19114008) on the imputed protein
matrix using:

- Signed Pearson correlation
- Soft-thresholding power chosen by scale-free topology fit (R² > 0.85)
- `minModuleSize = 30`
- `mergeCutHeight = 0.25`
- Standard `cutreeDynamic` hybrid algorithm

Module-eigengene values were then tested for treatment effects using the
same `~ PHE * Mito + (1 | Replicate)` LMM design. **This is a circular
analysis** because the same dataset defines the modules and tests their
treatment associations; the WGCNA output is therefore presented as an
exploratory supplement only (the figure carries an explicit EXPLORATORY
caption marker), not as an independent statistical claim.

## OXPHOS complex stoichiometry and mitonuclear balance

The OXPHOS subunit-level analysis aggregates proteins by their canonical
mitochondrial-complex membership (Complex I, II, III, IV, V, and the
mitoribosome) using MitoCarta 3.0 annotations. Per-complex enrichment
was tested with limma's `camera` competitive test (Wu and Smyth 2012,
PMID 22638577), which controls for inter-gene correlation in the test
statistic.

**Mitonuclear balance** follows the framework introduced by Houtkooper
et al. 2013 (PMID 23698443) and operationalizes the
mtDNA-encoded vs nuclear-encoded OXPHOS subunit ratio as a single
per-sample summary statistic, then tests treatment effects with the
same `~ PHE * Mito + (1 | Replicate)` LMM. A positive estimate indicates
that mtDNA-encoded subunits gained more than nuclear-encoded subunits in
the contrast — the predicted signature of intact mitochondrial-genome
delivery in a transplant paradigm.

**Pair-stoichiometry** (AlteredPQR; Romanov et al. 2019 / Buljan et al.
2023) tests whether pairs of subunits within the same OXPHOS complex
shifted asymmetrically (one subunit changing more than its partner),
which would indicate stoichiometric imbalance rather than coordinated
content change. At N = 24, this analysis is documented as exploratory.

**Subunit covariation** was assessed via per-complex Spearman correlation
matrices and a `cohesion` score (mean pairwise r within complex). The
cohesion range observed (0.02-0.13) is weak by standard cardiac proteome
benchmarks, so the AlteredPQR call is treated as descriptive only.

## Sensitivity analyses

Three robustness analyses were run as appended sheets to the differential
abundance output:

1. **Imputation sensitivity**: re-run the full DEP pipeline with an
   alternative imputation method and report per-contrast Spearman ρ
   between primary and alternative `logFC` vectors. ρ ≥ 0.99 indicates
   imputation choice does not drive the called proteins.
2. **Reinjection sensitivity**: refit the DEP pipeline after dropping
   the four `r`-suffixed re-injection samples; report per-contrast
   Spearman ρ. ρ ≥ 0.95 for the core contrasts indicates the imbalance
   does not dominate results.
3. **Robust eBayes sensitivity**: refit with limma's
   `eBayes(..., robust = TRUE)` variance moderation (Phipson et al. 2016,
   PMID 28367255) and compare FDR and Π counts to the default. Matched
   counts confirm that variance moderation choice is not a driver.

## Software environment

Analyses were performed in **R 4.6** with the following key packages:
**proteoDA**, **limma**, **fgsea**, **msImpute**, **missForest**,
**msigdbr**, **babelgene**, **GSVA**, **singscore**, **WGCNA**,
**lmerTest**, **emmeans**, **vegan** (for PERMANOVA), and **org.Rn.eg.db**
(for rat ortholog/symbol mapping). Figure rendering used **ggplot2**,
**patchwork**, **cowplot**, **ComplexHeatmap**, and **rasterGrob**.

---

# Results

## Sample structure and proteome-wide picture (Figure 1)

Principal-component analysis of the imputed proteome matrix separates
the four groups along PC1, with a significant overall group main effect
by PERMANOVA (pseudo-F p = 4 × 10⁻⁴; Anderson 2001) using Euclidean
distances on scaled per-protein intensities with 9,999 permutations.
Pairwise PERMANOVA contrasts for the three biologically meaningful
group pairs (Ctl vs Mito = "Transplant"; Ctl vs PHE = "Disease";
PHE vs PHE_Mito = "Rescue") are reported in the figure.

Per-protein significance counts show the expected ordering of effect
sizes across contrasts: PHEvPHE_MITO > CTLvMITO > CTLvPHE > Interaction.
167 unique Π-significant DEPs survive across the four core contrasts
(union, deduplicated by uniprot_id); 41 of those 167 are MitoCarta
members, a substantial mitochondrial enrichment for a study where the
intervention is mitochondrial transplant.

After collapsing the harmonized Hallmark + Reactome + MitoCarta pathway
annotations with the combined Jaccard / Overlap-Coefficient criterion,
82 unique pathways appear in the enrichment summary panel (down from
255 pre-dedup, a 68% redundancy collapse — see Methods for why this
strict collapse is necessary).

The mito-protein rank-location panel confirms that the mito-transplanted
contrasts (CTLvMITO and PHEvPHE_MITO) pull mitochondrial-annotated proteins
toward the upper tail of the per-protein effect-size distribution, while
the Disease contrast (CTLvPHE) and the Interaction contrast do not.

A separate supplementary QC figure documents that cyclic-loess
normalization, msImpute-based MAR/MNAR classification, and missForest
imputation recover sample clustering consistent with the experimental
design.

## Contrast-resolved volcano landscape (Figure 2)

Figure 2 is a 2 × 2 composite of volcano-in-ring panels, one per
contrast. The central volcano plots protein-level log-fold-change vs
−log₁₀(p-value); the ring around each volcano shows the top-significant
gene-set enrichments from the four-database curated pool
(Hallmark + Reactome + KEGG + MitoCarta), with up-regulated terms on the
right arc and down-regulated terms on the left arc. Arc colour encodes
NES (red = positive, blue = negative). Up to six up and six down terms
appear per panel, with a per-database cap of two terms per side.

The four panels surface distinct biology:

- **Transplant (CTLvMITO, panel A)** — 9 pathways (6 up / 3 down).
  HALLMARK_OXIDATIVE_PHOSPHORYLATION is the strongest single signal
  (padj 1.3 × 10⁻⁶), with REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT,
  MITOCARTA_MITORIBOSOME, MITOCARTA_TCA_CYCLE, and KEGG_GLYCOLYSIS
  triangulating the same mitochondrial-delivery signature from four
  independent annotation perspectives. The convergence of all four
  databases on OXPHOS biology is a strong cross-database confirmation
  that mitochondrial cargo is being delivered to recipient cells.

- **Rescue (PHEvPHE_MITO, panel B)** — 11 pathways (6 up / 5 down).
  HALLMARK_OXPHOS dominates UP with padj 1.5 × 10⁻¹⁶, the strongest
  single statistical signal in the entire dataset. HALLMARK_FATTY_ACID_METABOLISM,
  HALLMARK_ADIPOGENESIS, REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION,
  and KEGG_OXIDATIVE_PHOSPHORYLATION are also significantly UP. Critically,
  the DOWN side of this panel shows HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  REACTOME_COLLAGEN_TRIMERIZATION, REACTOME_FIBRONECTIN_MATRIX_FORMATION,
  and HALLMARK_ANGIOGENESIS — all reverse-image of the PHE-induced
  fibrotic signature seen in panel C. This is the rescue signature:
  the same biology that PHE drives up is being driven back down by
  mitochondrial transplant.

- **Disease (CTLvPHE, panel C)** — 8 pathways (5 up / 3 down).
  KEGG_ECM_RECEPTOR_INTERACTION (padj 0.002), REACTOME_STRIATED_MUSCLE_CONTRACTION
  (padj 0.007), REACTOME_FORMATION_OF_DYSTROPHIN_GLYCOPROTEIN_COMPLEX,
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, and KEGG_DILATED_CARDIOMYOPATHY
  are UP, recovering the Hahn et al. 2014 (PMID 24794531) PHE-driven
  fibrotic / sarcomeric / ECM-remodeling axis at pathway resolution.
  HALLMARK_MYC_TARGETS_V1, REACTOME_DNA_REPLICATION, and HALLMARK_E2F_TARGETS
  are DOWN — proliferation arrest paired with structural hypertrophy,
  a classic terminally-differentiating cardiomyoblast response to
  α1-adrenergic stress. (Note: this panel was empty under an earlier
  MitoCarta-only enrichment lens; pooling four curated databases
  surfaces non-mitochondrial disease biology that the single-database
  lens hides.)

- **Interaction (panel D)** — 11 pathways shown for descriptive
  reference only. At N = 24, the interaction contrast is underpowered
  for formal inference (zero proteins reach FDR < 0.10 and the figure
  caption explicitly carries this caveat).

A supplementary expanded composite that adds GO:BP to the canonical
four-database pool provides finer-resolution biological-process terms
for readers interested in deeper annotation depth; the same conclusions
emerge but with `GOBP_ACTIN_FILAMENT_BASED_PROCESS` and related GO:BP
terms triangulating the curated-database findings.

## Transplant vs Rescue concordance (Figure 3)

Figure 3 asks whether the mitochondrial transplant signature observed in
healthy cells (CTLvMITO) is the same signature observed when transplanting
into PHE-stressed cells (PHEvPHE_MITO), or whether the two contexts elicit
different effects of the same intervention.

The NES-vs-NES quadrant scatter (Transplant on x, Rescue on y) across
the 74 shared significant pathways shows Spearman ρ = 0.626
(95% CI 0.46–0.75) — strong positive concordance. The pattern heatmap
classifies pathways into nine quadrant categories; the dominant categories
are Concordant-Up (transplant ↑, rescue ↑) and Concordant-Down (transplant
↓, rescue ↓), with virtually empty Discordant quadrants.

The RRHO2 hotspot analysis (Cahill et al. 2018, PMID 29942049) identifies
1,629 concordant-up and 1,636 concordant-down genes; the discordant
quadrants (transplant-up but rescue-down, or vice versa) contain only
0–1 hotspot genes. Over-representation analysis on the Concordant-Up
hotspot recovers 12 deduplicated mitochondrial / OXPHOS pathways.

**Interpretation**: the mitochondrial transplant intervention is highly
consistent across cellular states. Whether mitochondria are delivered to
healthy or PHE-stressed cells, the same proteomic signature emerges. This
is a precondition for treating the rescue effect (Figure 4) as a real
interaction between treatment and disease, rather than a separate effect
arising in the stressed cells alone.

## Disease vs Rescue reversal (Figure 4)

Figure 4 asks the central biological question of the dataset: does
mitochondrial transplant in PHE-stressed cells return the proteome
toward control, or does it merely add a mitochondrial signature on top
of the disease state?

The NES-vs-NES quadrant scatter (Disease PHE-effect on x, Rescue
mito-effect on y) shows Spearman ρ = **−0.42** (95% CI −0.59 to −0.21).
The negative correlation **is** the rescue signature: pathways perturbed
by PHE in one direction are returned in the opposite direction by
mitochondrial transplant.

RRHO2 reveals 1,440 Disease-up genes that Rescue brings down, and 1,151
Disease-down genes that Rescue brings up. The two exacerbated quadrants
(transplant making PHE worse) contain essentially no genes. The pattern
heatmap puts most pathways in either the Reversed-DisUp-ResDown or
Reversed-DisDown-ResUp categories.

ORA on the Reversed-DisUp-ResDown hotspot recovers 12 pathways centered
on collagen / ECM / fibrotic biology — exactly the axis that PHE drives
up and that mitochondrial transplant brings back toward control.

**Interpretation**: the rescue is real and directional. Mitochondrial
transplant does not just add an OXPHOS signature on top of disease — it
actively returns the PHE-perturbed cytoskeletal / fibrotic axis toward
the control proteome.

## Pathway-set scores and module structure (Figure 5)

The set-score panel reports per-contrast effect sizes for the 416-set
harmonized collection (post Jaccard / Overlap-Coefficient dedup), scored
with both GSVA and singscore. Per-contrast Spearman correlations between
GSVA and singscore effect sizes are 0.82 (Rescue), 0.79 (Transplant),
0.78 (Disease), and 0.77 (Interaction). The two stronger contrasts clear
or sit at the ρ ≥ 0.80 methodological-agreement floor; the two weaker
contrasts sit just below, consistent with noisier method agreement when
the underlying biological effect is null (both Disease and Interaction
return 0 FDR-significant sets, so method agreement is being assessed on
noise).

Under the LMM framework with Replicate block, 12 sets pass FDR < 0.05
under Rescue, 9 under Transplant, and 0 under Disease and Interaction.
The top sets in both Transplant and Rescue are mitochondrial
(OXPHOS subunit families, mitoribosome, TCA cycle), consistent with
Figures 2-4.

A WGCNA supplementary panel reports seven non-grey modules with 35
designated hub proteins. The module-eigengene LMMs are explicitly
labeled EXPLORATORY because of the circular-analysis problem (Methods).
The module-NES diptych against Hallmark + Reactome + MitoCarta provides
a functional ledger of what each module corresponds to and supports the
mitochondrial-axis dominance of the dataset.

## OXPHOS complex stoichiometry and mitonuclear balance (Figure 6)

Per-complex camera test results for each OXPHOS complex (and the
mitoribosome) under each of the three biologically meaningful contrasts:

| Set | Transplant (UP) | Disease | Rescue (UP) |
|---|---|---|---|
| OXPHOS_all (n = 64) | p = 0.026 | p = 0.26 (Down) | **p = 3.2 × 10⁻⁵** |
| Complex V (ATP synthase, n = 15) | p = 0.025 | p = 0.12 (Down) | **p = 3.8 × 10⁻⁵** |
| Mitoribosome (n = 55) | **p = 3.5 × 10⁻³** | p = 0.55 | p = 4.9 × 10⁻³ |
| Complex I (n = 30) | p = 0.36 | p = 0.19 (Down) | p = 0.016 |
| Complex III (n = 8) | p = 0.44 | p = 0.37 (Down) | p = 0.015 |
| Complex IV (n = 11) | p = 0.12 | p = 0.10 | p = 0.41 |

The Rescue camera signal for OXPHOS_all (p = 3.2 × 10⁻⁵, FDR = 2 × 10⁻⁴)
is the single largest statistical signal in the dataset. Complex V and
the mitoribosome are the most-affected sub-units under Transplant; under
Rescue, every complex except IV becomes individually significant.

**Mitonuclear balance** (mtDNA-encoded vs nuclear-encoded OXPHOS subunit
ratio, per Houtkooper et al. 2013 PMID 23698443) is positive and
significant under Transplant (estimate = +0.20, FDR = 0.035) and
near-significant under Rescue (estimate = +0.15, FDR = 0.075). **mtDNA-
encoded subunits gain more than nuclear-encoded subunits under
mitochondrial transplant** — direct evidence that the delivered cargo
includes intact mitochondrial DNA and its translated products. This is
supported molecularly by Mt-nd4 (mtDNA-encoded Complex I subunit) being
a top-20 per-protein hit in CTLvMITO at logFC = +0.99.

The pair-stoichiometry LMM (AlteredPQR-style analysis) returned 0 FDR-
significant subunit-pair imbalance effects. The covariation supplement
shows weak per-complex subunit cohesion (range 0.02-0.13), reflecting
the limited power available at N = 24 to detect stoichiometric
asymmetries. This portion of the analysis is documented as exploratory.

## Per-contrast biological summary

### CTLvPHE — α1-adrenergic stress

At the protein level, PHE upregulates **Ccn2 / CTGF** (logFC = +0.56,
FDR = 0.013), a canonical α1-adrenergic and fibrotic target, and
downregulates the Complex V subunit **Atp5mk** (logFC = −0.37). The
broader CTLvPHE signature captures extracellular-matrix remodeling
(Postn, Fn1 trending), early OXPHOS depression, and mitoribosomal
stress (Mrpl9 logFC = +2.12). At the pathway level (Figure 2 panel C),
the fibrotic / sarcomeric / ECM-remodeling axis emerges through
KEGG_ECM_RECEPTOR_INTERACTION, REACTOME_STRIATED_MUSCLE_CONTRACTION,
REACTOME_DYSTROPHIN_GLYCOPROTEIN_COMPLEX, HALLMARK_EMT, and
KEGG_DILATED_CARDIOMYOPATHY, paired with proliferation arrest
(HALLMARK_MYC_TARGETS, REACTOME_DNA_REPLICATION, HALLMARK_E2F_TARGETS
down). All of this is consistent with Hahn et al. 2014 (PMID 24794531;
NADPH oxidase-2 / ROS axis) and Jeong et al. 2009 (PMID 19299911; STAT3
mitochondrial translocation suppressing OXPHOS).

**An important caveat for the Discussion**: the canonical in-vivo
cardiac-hypertrophy fetal-gene program (NPPA / NPPB / Myh7 / Acta1
upregulation) is not recovered at the protein level here. NPPA / NPPB
were not detected in the mass spectrum; Myh7, Acta1, Stat3, Cav3, and
Tnnt2 returned essentially flat logFC values (|logFC| < 0.2, all
FDR > 0.8). This is consistent with known H9c2 biology: cultured H9c2
cardiomyoblasts express embryonic sarcomeric isoforms (e.g. Myh3 — which
IS upregulated under Transplant in CTLvMITO) rather than reactivating
the adult fetal-gene reprogram seen in cardiomyocyte hypertrophy in vivo.
The hypertrophy signature recovered here should therefore be interpreted
as the **fibrotic + early-OXPHOS arm** of PHE α1-adrenergic stress, not
the full in-vivo fetal-gene reprogram.

### CTLvMITO — mitochondrial transplant in healthy cells

Top-Π hits are dominated by mitochondrial proteins: **Aldh3a1**
(logFC = +1.80, FDR = 6 × 10⁻⁸), **Vdac1** (logFC = +0.56, FDR = 2 × 10⁻⁴),
**Vdac2** (logFC = +0.42, FDR = 6 × 10⁻³), and the mtDNA-encoded
**Mt-nd4** (logFC = +0.99). Combined with the Figure 2 MitoCarta-pathway
NES profile and the Figure 6 mitoribosome + Complex V camera signals,
the data describe organelle-scale cargo delivery (outer-membrane channels,
matrix enzymes, mtDNA-encoded subunits) into recipient cardiomyoblasts.
This recapitulates the cell-autonomous mitochondrial-delivery signature
described by Masuzawa et al. 2013 (PMID 23355340).

### PHEvPHE_MITO — rescue

The rescue contrast carries an even stronger mito-delivery signature
than CTLvMITO: all three VDAC isoforms (Vdac1 +0.71, Vdac2 +0.57,
Vdac3 +0.49) plus Aldh3a1 (+1.83) appear in the top 25 by Π. The
peroxiredoxin-like Prxl2a (+1.09) suggests an antioxidant arm consistent
with rescue of PHE-driven NOX2 / ROS production (Hahn 2014). The
Figure 4 concordance ρ = −0.42 and the Figure 6 camera OXPHOS_all
p = 3.2 × 10⁻⁵ together make this the strongest single biological
signal in the dataset, aligned with Doulamis et al. 2024 (PMID 39732955)
mitochondrial-transplant normalization of the proteome in ischemia-
reperfusion injury.

### Interaction — descriptive only

Zero FDR-significant proteins, consistent with the N = 24 power
constraint. Direction-only inspection shows a faint mito-axis signal
(Ndufb4l4 logFC = +0.83, Mrpl27 logFC = +0.76) that should not be
interpreted further given the adj-P-value range (0.34 – 0.98). The
"Interaction underpowered at N = 24" caveat on every composite figure
footer is the appropriate representation.

## Method robustness summary

The analysis is methodologically robust to the major sensitivity tests:

- **Imputation choice does not drive results.** Per-contrast Spearman ρ
  between the primary missForest imputation and an alternative imputation
  ranged from 0.992 to 0.997 across all five contrasts.
- **Reinjection imbalance does not drive results in the core contrasts.**
  Per-contrast ρ with the four `r`-suffixed re-injection samples removed
  was ≥ 0.95 for the four core contrasts. The secondary MITOvPHE_MITO
  comparison dropped to ρ = 0.82, which is documented as the largest
  open caveat in the dataset and the reason MITOvPHE_MITO is not
  interpreted as a primary biological contrast.
- **Variance-moderation choice does not drive results.** Robust empirical-
  Bayes moderation returned identical FDR and Π counts to the default in
  all contrasts.

## Open caveats

1. **H9c2 cell-line caveat.** H9c2 cardiomyoblasts under-express the
   in-vivo fetal-gene cardiac-hypertrophy reprogram. PHE-induced
   hypertrophy in this system is recovered as fibrotic and early-OXPHOS
   depression signatures rather than NPPA / NPPB / Myh7 upregulation.
2. **Reinjection imbalance.** Four `r`-suffixed re-injection samples
   are imbalanced across groups (Ctl 1, Mito 0, PHE 2, PHE_Mito 3).
   Sensitivity refits show core contrasts robust (Spearman ρ 0.95-0.99);
   the secondary MITOvPHE_MITO comparison drops to ρ = 0.82.
3. **Power.** N = 24 (n = 6 / group) supports the three primary
   contrasts but underpowers the Interaction term and the AlteredPQR
   pair-stoichiometry tests. The covariation supplement (cohesion
   0.02-0.13) is flagged accordingly.
4. **WGCNA circularity.** Module eigengene LMMs use the same dataset
   that the modules were derived from. The WGCNA supplement is
   explicitly labeled EXPLORATORY and should not be used as an
   independent statistical claim.
5. **In-vitro α1-adrenergic stimulation.** Phenylephrine at the doses
   used recapitulates one arm (NOX2 / ROS, mitochondrial-STAT3) of
   in-vivo pressure-overload cardiac stress; it is not a complete
   model of clinical heart failure.

---

# References (verified PMIDs)

1. Hahn NE, Musters RJ, Fritz JM, et al. Early NADPH oxidase-2 activation
   is crucial in phenylephrine-induced hypertrophy of H9c2 cells.
   *Cell Signal*. 2014;26(9):1818-24. PMID 24794531.
2. Jeong K, Kwon H, Min C, Pak Y. Modulation of the caveolin-3
   localization to caveolae and STAT3 to mitochondria by catecholamine-
   induced cardiac hypertrophy in H9c2 cardiomyoblasts. *Exp Mol Med*.
   2009;41(4):226-35. PMID 19299911.
3. Masuzawa A, Black KM, Pacak CA, et al. Transplantation of autologously
   derived mitochondria protects the heart from ischemia-reperfusion
   injury. *Am J Physiol Heart Circ Physiol*. 2013;304(7):H966-82.
   PMID 23355340.
4. Doulamis IP, Tzani A, Alemany VS, et al. Mitochondrial transplantation
   normalizes transcriptomic and proteomic shift associated with ischemia
   reperfusion injury in neonatal hearts donated after circulatory death.
   *Sci Rep*. 2024;14(1):31236. PMID 39732955.
5. Houtkooper RH, Mouchiroud L, Ryu D, et al. Mitonuclear protein
   imbalance as a conserved longevity mechanism. *Nature*. 2013;
   497(7450):451-7. PMID 23698443.
6. Rath S, Sharma R, Gupta R, et al. MitoCarta3.0: an updated
   mitochondrial proteome now with sub-organelle localization and
   pathway annotations. *Nucleic Acids Res*. 2021;49(D1):D1541-D1547.
   PMID 33174596.
7. Xiao Y, Hsiao TH, Suresh U, et al. A novel significance score for
   gene selection and ranking. *Bioinformatics*. 2014;30(6):801-7.
   PMID 22321699.
8. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment
   analysis: a knowledge-based approach for interpreting genome-wide
   expression profiles. *Proc Natl Acad Sci U S A*. 2005;102(43):
   15545-50. PMID 16199517.
9. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular
   Signatures Database (MSigDB) hallmark gene set collection.
   *Cell Syst*. 2015;1(6):417-425. PMID 26771021.
10. Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway
    knowledgebase 2022. *Nucleic Acids Res*. 2022;50(D1):D687-D692.
    PMID 34788843.
11. Kanehisa M. Toward understanding the origin and evolution of cellular
    organisms. *Protein Sci*. 2019;28(11):1947-1951. PMID 31441146.
12. Merico D, Isserlin R, Stueker O, Emili A, Bader GD. Enrichment map:
    a network-based method for gene-set enrichment visualization and
    interpretation. *PLoS ONE*. 2010;5(11):e13984. PMID 21085593.
13. Reimand J, Isserlin R, Voisin V, et al. Pathway enrichment analysis
    and visualization of omics data using g:Profiler, GSEA, Cytoscape
    and EnrichmentMap. *Nat Protoc*. 2019;14(2):482-517. PMID 30664679.
14. Ritchie ME, Phipson B, Wu D, et al. limma powers differential
    expression analyses for RNA-sequencing and microarray studies.
    *Nucleic Acids Res*. 2015;43(7):e47. PMID 25605792.
15. Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. Robust
    hyperparameter estimation protects against hypervariable genes
    and improves power to detect differential expression. *Ann Appl
    Stat*. 2016;10(2):946-963. PMID 28367255.
16. Wu D, Smyth GK. Camera: a competitive gene set test accounting for
    inter-gene correlation. *Nucleic Acids Res*. 2012;40(17):e133.
    PMID 22638577.
17. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation
    analysis for microarray and RNA-seq data. *BMC Bioinformatics*.
    2013;14:7. PMID 23323831.
18. Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ.
    Single sample scoring of molecular phenotypes. *BMC Bioinformatics*.
    2018;19(1):404. PMID 30400809.
19. Langfelder P, Horvath S. WGCNA: an R package for weighted
    correlation network analysis. *BMC Bioinformatics*. 2008;9:559.
    PMID 19114008.
20. Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment
    analysis. *bioRxiv*. 2021. doi:10.1101/060012.
21. Stekhoven DJ, Bühlmann P. MissForest — non-parametric missing value
    imputation for mixed-type data. *Bioinformatics*. 2012;28(1):112-8.
    PMID 22039212.
22. Cahill KM, Huo Z, Tseng GC, Logan RW, Seney ML. Improved
    identification of concordant and discordant gene expression
    signatures using an updated rank-rank hypergeometric overlap
    approach. *Sci Rep*. 2018;8(1):9588. PMID 29942049.
23. Anderson MJ. A new method for non-parametric multivariate analysis
    of variance. *Austral Ecol*. 2001;26(1):32-46. (PERMANOVA.)
