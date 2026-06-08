---
title: "A_Mito_2026 — Results and interpretation (candidate manuscript text)"
date: "2026-06-04"
---

# Overview

H9c2 rat cardiomyoblasts treated with phenylephrine (PHE) and/or co-incubated
with isolated mitochondria, in a 2 × 2 design (Ctl, Mito, PHE, PHE_Mito;
n = 6/group, N = 24). The four core comparisons map onto the rescue paradigm
introduced by Masuzawa and colleagues (PMID 23355340) and most recently
extended by Doulamis and colleagues (PMID 39732955):

| Contrast | Biological question |
|---|---|
| CTLvPHE | Does α1-adrenergic stress remodel the cardiomyoblast proteome? |
| CTLvMITO | Does mitochondrial transplant leave a cell-autonomous signature in healthy cells? |
| PHEvPHE_MITO | Does mitochondrial transplant rescue the PHE-stressed proteome? |
| Interaction | Does PHE alter what mitochondrial transplant does? (descriptive; underpowered at N = 24) |

A secondary MITOvPHE_MITO comparison was retained for sensitivity reporting
but is not interpreted as a primary biological contrast.

Stage 03 differential abundance used limma with `~ 0 + group + (1 | Replicate)`,
routed through proteoDA's duplicateCorrelation interface. The capital-pi
significance score (Π = P-value^|logFC|; Xiao et al. 2014, Bioinformatics
30:801) was applied at Π < 0.05 alongside Benjamini-Hochberg FDR < 0.10.
Pathway enrichment (fgsea) was harmonized across the curated four-database
pool Hallmark + Reactome + KEGG + MitoCarta (rat-mapped). Redundancy
across databases was collapsed using the combined Jaccard / Overlap-
Coefficient criterion at threshold 0.5 (Merico et al. 2010, PLoS ONE
PMID 21085593; Reimand et al. 2019, Nat Protoc PMID 30664679, Box 1) —
the Jaccard pass handles symmetric overlap, and the Overlap-Coefficient
pass catches asymmetric containment where a small specialized set is
entirely inside a larger general set.

# Method robustness

Imputation sensitivity (per-contrast Spearman ρ between the primary missForest
matrix and an alternative imputation) was ≥ 0.992 for all five contrasts.
Reinjection sensitivity (per-contrast Spearman ρ with the four r-suffixed
samples removed) was ≥ 0.95 for the four core contrasts and 0.82 for the
secondary MITOvPHE_MITO comparison; the imbalance (Ctl 1, Mito 0, PHE 2,
PHE_Mito 3) is the largest open caveat in the dataset. Robust empirical-Bayes
moderation returned identical FDR and Π counts to the default in all
contrasts, indicating that variance moderation choice is not a driver of the
called proteins.

# Figure-by-figure interpretation

## F01 — Sample structure and proteome-wide picture

Principal-component analysis separates the four groups along PC1
(group main effect: PERMANOVA pseudo-F p = 4 × 10⁻⁴). Per-protein
significance shows the expected ordering of effect sizes:
PHEvPHE_MITO > CTLvMITO > CTLvPHE > Interaction. 167 unique
Π-significant DEPs survive across the four core contrasts; 41 are
MitoCarta members. After harmonizing the Hallmark + Reactome + MitoCarta
pathway annotations and collapsing redundancy with the combined
Jaccard / Overlap-Coefficient criterion (Methods), 82 unique pathways
appear in the enrichment panel. The mito-protein rank panel confirms
that the mito-transplanted contrasts pull mitochondrial proteins to
the upper tail of the effect-size distribution.

The pipeline supplement (F01_QC_supplementary.xlsx) shows that normalization
(cycloess), missingness classification (msImpute EBM into MAR/MNAR), and
imputation (single missForest applied to all proteins) recover sample
clustering consistent with the design.

## F02 — Contrast-resolved volcano landscape (curated 4-DB pooled canonical, +GO:BP supplement)

The 2 × 2 volcano-in-ring composite pools four curated databases
(Hallmark + Reactome + KEGG + MitoCarta, rat-mapped) for the canonical
main figure, with an expanded supplement adding GO:BP for depth. Each
ring caps individual database contribution at two terms to prevent any
single database from dominating visual real-estate (Methods). Post-dedup
counts shown in the rings:

- **Transplant (CTLvMITO)**: 9 pathways (6 up / 3 down). HALLMARK_OXIDATIVE_PHOSPHORYLATION is the strongest single signal (padj 1.3 × 10⁻⁶), with REACTOME_AEROBIC_RESPIRATION, MITOCARTA_MITORIBOSOME, MITOCARTA_TCA_CYCLE, and KEGG_GLYCOLYSIS triangulating the mitochondrial-delivery signature from four independent annotation perspectives.
- **Rescue (PHEvPHE_MITO)**: 11 pathways (6 up / 5 down). HALLMARK_OXPHOS dominates UP (padj 1.5 × 10⁻¹⁶, the strongest single signal in the dataset), with HALLMARK_FATTY_ACID_METABOLISM and REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION also UP. HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, REACTOME_COLLAGEN_TRIMERIZATION, and FIBRONECTIN_MATRIX are DOWN — the reverse-image of the PHE-induced fibrotic signature (rescue of disease).
- **Disease (CTLvPHE)**: 8 pathways (5 up / 3 down) — populated for the first time under the pooled lens. KEGG_ECM_RECEPTOR_INTERACTION (padj 0.002), REACTOME_STRIATED_MUSCLE_CONTRACTION (padj 0.007), REACTOME_FORMATION_OF_DYSTROPHIN_GLYCOPROTEIN_COMPLEX, HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, and KEGG_DILATED_CARDIOMYOPATHY are UP, recovering the Hahn 2014 (PMID 24794531) PHE-driven fibrotic / sarcomeric / ECM-remodeling axis at pathway resolution. HALLMARK_MYC_TARGETS_V1, REACTOME_DNA_REPLICATION, and HALLMARK_E2F_TARGETS are DOWN — proliferation arrest paired with hypertrophy, classic for terminally-differentiating cardiomyoblasts under α1-adrenergic stress.
- **Interaction**: 11 pathways (descriptive only at N = 24).

The expanded supplement (Figure S2, adds GO:BP) provides depth for
readers interested in finer-resolution biological-process terms;
GO:BP-tagged actin / sarcomere / muscle-organ-development terms
triangulate the curated-DB findings for the Disease panel.

## F03 — Transplant vs Rescue concordance

The NES-vs-NES quadrant (Transplant on x, Rescue on y) shows Spearman
ρ = 0.626 (95% CI 0.46–0.75) across the 74 shared pathways. The RRHO2
matrix identifies 1,629 concordant-up and 1,636 concordant-down genes
with effectively zero discordant content. Concordant-up ORA (after 45.7%
Jaccard dedup) recovers 12 mitochondrial / OXPHOS pathways. **The
mitochondrial transplant signature observed in healthy cells (CTLvMITO) is
largely the same signature that is also observed when transplanting into
PHE-stressed cells (PHEvPHE_MITO).** The intervention is consistent across
cellular states.

## F04 — Disease vs Rescue reversal

The NES quadrant (Disease PHE-effect on x, Rescue mito-effect on y)
shows Spearman ρ = −0.42 (CI −0.59 to −0.21). The negative correlation
is the **rescue signature**: pathways perturbed by PHE in one direction
are returned in the opposite direction by mitochondrial transplant.
RRHO2 reveals 1,440 Disease-up genes that Rescue brings down, and 1,151
Disease-down genes that Rescue brings up. The exacerbated quadrants
(Mito making PHE worse) are essentially empty. This is the central
biological claim of the dataset: mitochondrial transplant does not just
add a mito signature on top of disease — it actively returns the
PHE-perturbed proteome toward control.

## F05 — Pathway-set scores and module structure

Set-level scoring with both GSVA and singscore (per-contrast Spearman ρ
0.77–0.82) shows that 12 sets pass FDR < 0.05 under PHE+Mito Rescue,
9 under Transplant, 0 under Disease and Interaction (after stricter
combined Jaccard / Overlap-Coefficient dedup of the harmonized
collection, 1,625 → 416 sets). The two stronger contrasts clear or
sit at the methodological-agreement floor (ρ ≥ 0.80); the two weaker
contrasts sit just below (0.78 and 0.77), consistent with noisier
method agreement when underlying effects are null.

The WGCNA supplement is exploratory only (the module-eigengene LMMs use
the same data the modules were derived from; the circular-analysis caveat
on the figure should not be removed). Seven non-grey modules with 35
designated hub proteins emerged, and the module-NES diptych against
Hallmark + Reactome + MitoCarta provides a coarse functional ledger
of what each module corresponds to.

## F06 — OXPHOS complex stoichiometry and mitonuclear balance

OXPHOS subunit content (`camera` competitive test on grouped subunits):

| Set | Transplant (UP) | Disease | Rescue (UP) |
|---|---|---|---|
| OXPHOS_all (n = 64) | p = 0.026 | p = 0.26 (Down) | **p = 3.2 × 10⁻⁵** |
| Complex V (n = 15) | p = 0.025 | p = 0.12 (Down) | **p = 3.8 × 10⁻⁵** |
| Mitoribosome (n = 55) | **p = 3.5 × 10⁻³** | p = 0.55 | p = 4.9 × 10⁻³ |
| Complex I (n = 30) | p = 0.36 | p = 0.19 (Down) | p = 0.016 |
| Complex III (n = 8) | p = 0.44 | p = 0.37 (Down) | p = 0.015 |
| Complex IV (n = 11) | p = 0.12 | p = 0.10 | p = 0.41 |

The Rescue camera signal for OXPHOS_all (p = 3.2 × 10⁻⁵, FDR = 2 × 10⁻⁴)
is the single largest statistical signal in the dataset. Complex V and
the mitoribosome are the most affected sub-units under Transplant; under
Rescue, every complex except IV becomes significant.

Mitonuclear balance (mtDNA-encoded vs nuclear-encoded OXPHOS subunits;
Houtkooper et al. 2013, Nature, PMID 23698443) is positive and significant
under Transplant (estimate = +0.20, FDR = 0.035) and near-significant under
Rescue (estimate = +0.15, FDR = 0.075). **mtDNA-encoded subunits gain more
than nuclear-encoded subunits under mitochondrial transplant**, providing
direct evidence that the delivered cargo includes intact mitochondrial DNA
and its translated products (Mt-nd4, top-20 per-protein hit in CTLvMITO at
logFC = +0.99, supports this molecularly).

The pair-stoichiometry LMM returned 0 FDR-significant subunit-pair
imbalance effects. The covariation supplement shows weak cohesion
(range 0.02-0.13) and is documented as underpowered at N = 24; this is
why the AlteredPQR proteostatic-rebalance call is treated as exploratory.

# Per-contrast biological summary

## CTLvPHE — α1-adrenergic stress

PHE upregulates Ccn2 / CTGF (logFC = +0.56, FDR = 0.013), a canonical
α1-adrenergic and fibrotic target, and downregulates the Complex V
subunit Atp5mk (logFC = −0.37). The broader CTLvPHE signature
captures extracellular-matrix remodeling (Postn, Fn1 trending), early
OXPHOS depression, and mitoribosomal stress (Mrpl9 logFC = +2.12), all
consistent with Hahn et al. 2014 (PMID 24794531; NADPH oxidase-2 / ROS axis)
and Jeong et al. 2009 (PMID 19299911; STAT3 mitochondrial translocation
suppressing OXPHOS).

**Limitation that should be carried into Discussion.** The canonical
in-vivo cardiac-hypertrophy fetal-gene program (NPPA / NPPB / Myh7 / Acta1
upregulation) is not recovered at the protein level here: NPPA / NPPB were
not detected, and Myh7 / Acta1 / Stat3 / Cav3 / Tnnt2 returned essentially
flat logFC values (|logFC| < 0.2, all FDR > 0.8). This is consistent with
known H9c2 biology: cultured H9c2 cardiomyoblasts express embryonic
sarcomeric isoforms (e.g. Myh3 — which IS upregulated under Transplant)
rather than reactivating the adult fetal-gene reprogram seen in cardiomyocyte
hypertrophy in vivo. The hypertrophy signature recovered here should therefore
be interpreted as the **fibrotic + early-OXPHOS arm** of PHE α1-adrenergic
stress, not the full in-vivo fetal-gene reprogram.

## CTLvMITO — mitochondrial transplant in healthy cells

Top-Π hits are dominated by mitochondrial proteins: Aldh3a1 (logFC = +1.80,
FDR = 6 × 10⁻⁸), Vdac1 (logFC = +0.56, FDR = 2 × 10⁻⁴), Vdac2 (logFC = +0.42,
FDR = 6 × 10⁻³), and the mtDNA-encoded Mt-nd4 (logFC = +0.99). Combined with
F02's MitoCarta-pathway NES profile and F06's Mitoribosome + Complex V camera
signal, the data describe organelle-scale cargo delivery — outer membrane
channels, matrix enzymes, mtDNA-encoded subunits — into recipient cardiomyoblasts.
This recapitulates the cell-autonomous mitochondrial-delivery signature
described by Masuzawa et al. 2013 (PMID 23355340).

## PHEvPHE_MITO — rescue

The rescue contrast carries an even stronger mito-delivery signature than
CTLvMITO: all three VDAC isoforms (Vdac1 + 0.71, Vdac2 + 0.57, Vdac3 + 0.49)
plus Aldh3a1 (+1.83) appear in the top 25 by Π. The peroxiredoxin-like
Prxl2a (+1.09) suggests an antioxidant arm consistent with rescue of
PHE-driven NOX2 / ROS production (Hahn 2014). The F04 concordance ρ = −0.42
and F06 camera OXPHOS_all p = 3.2 × 10⁻⁵ make this the strongest single
biological signal in the dataset and align with Doulamis et al. 2024
(PMID 39732955) mitochondrial-transplant normalization of the proteome
in ischemia-reperfusion injury.

## Interaction — descriptive only

Zero FDR-significant proteins, consistent with the N = 24 power constraint.
Direction-only inspection shows a faint mito-axis signal (Ndufb4l4 logFC = +0.83,
Mrpl27 logFC = +0.76) that should not be interpreted further given the
adj-P-value range (0.34 - 0.98). The "Interaction underpowered at N = 24"
caveat (M4) on every composite figure footer is the appropriate
representation.

# Open caveats for the manuscript Discussion / Limitations section

1. **H9c2 cell-line caveat.** H9c2 cardiomyoblasts under-express the in-vivo
   fetal-gene hypertrophy reprogram. PHE-induced hypertrophy is recovered
   as fibrotic and early-OXPHOS depression signatures rather than NPPA / NPPB
   / Myh7 upregulation.

2. **Reinjection imbalance.** Four r-suffixed re-injection samples are
   imbalanced across groups (Ctl 1, Mito 0, PHE 2, PHE_Mito 3). Sensitivity
   refits show core contrasts robust (Spearman ρ 0.95-0.99); the secondary
   MITOvPHE_MITO comparison drops to ρ = 0.82.

3. **Power.** N = 24 (n = 6 / group) supports the three primary contrasts
   but underpowers the Interaction term and the AlteredPQR pair-stoichiometry
   tests. The covariation supplement (cohesion 0.02-0.13) is flagged
   accordingly.

4. **WGCNA circularity.** Module eigengene LMMs use the same dataset that
   the modules were derived from. The WGCNA supplement is explicitly
   labeled EXPLORATORY and should not be used as an independent statistical
   claim in the main results.

5. **In-vitro α1-adrenergic stimulation.** Phenylephrine at the doses used
   recapitulates one arm (NOX2 / ROS, mitochondrial-STAT3) of in-vivo
   pressure-overload cardiac stress; it is not a complete model of clinical
   heart failure.

# References (verified PMIDs)

1. Hahn NE, Musters RJ, Fritz JM, et al. Early NADPH oxidase-2 activation is
   crucial in phenylephrine-induced hypertrophy of H9c2 cells.
   Cell Signal. 2014;26(9):1818-24. PMID 24794531.
2. Jeong K, Kwon H, Min C, Pak Y. Modulation of the caveolin-3 localization
   to caveolae and STAT3 to mitochondria by catecholamine-induced cardiac
   hypertrophy in H9c2 cardiomyoblasts. Exp Mol Med. 2009;41(4):226-35.
   PMID 19299911.
3. Masuzawa A, Black KM, Pacak CA, et al. Transplantation of autologously
   derived mitochondria protects the heart from ischemia-reperfusion injury.
   Am J Physiol Heart Circ Physiol. 2013;304(7):H966-82. PMID 23355340.
4. Doulamis IP, Tzani A, Alemany VS, et al. Mitochondrial transplantation
   normalizes transcriptomic and proteomic shift associated with ischemia
   reperfusion injury in neonatal hearts donated after circulatory death.
   Sci Rep. 2024;14(1):31236. PMID 39732955.
5. Houtkooper RH, Mouchiroud L, Ryu D, et al. Mitonuclear protein imbalance
   as a conserved longevity mechanism. Nature. 2013;497(7450):451-7.
   PMID 23698443.
6. Rath S, Sharma R, Gupta R, et al. MitoCarta3.0: an updated mitochondrial
   proteome now with sub-organelle localization and pathway annotations.
   Nucleic Acids Res. 2021;49(D1):D1541-D1547. PMID 33174596.
7. Xiao Y, Hsiao TH, Suresh U, et al. A novel significance score for gene
   selection and ranking. Bioinformatics. 2014;30(6):801-7. PMID 22321699.
8. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis:
   a knowledge-based approach for interpreting genome-wide expression
   profiles. Proc Natl Acad Sci U S A. 2005;102(43):15545-50. PMID 16199517.
9. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures
   Database (MSigDB) hallmark gene set collection. Cell Syst.
   2015;1(6):417-425. PMID 26771021.
10. Gillespie M, Jassal B, Stephan R, et al. The Reactome pathway
    knowledgebase 2022. Nucleic Acids Res. 2022;50(D1):D687-D692.
    PMID 34788843.
11. Kanehisa M. Toward understanding the origin and evolution of cellular
    organisms. Protein Sci. 2019;28(11):1947-1951. PMID 31441146.
12. Merico D, Isserlin R, Stueker O, Emili A, Bader GD. Enrichment map:
    a network-based method for gene-set enrichment visualization and
    interpretation. PLoS ONE. 2010;5(11):e13984. PMID 21085593.
13. Reimand J, Isserlin R, Voisin V, et al. Pathway enrichment analysis
    and visualization of omics data using g:Profiler, GSEA, Cytoscape
    and EnrichmentMap. Nat Protoc. 2019;14(2):482-517. PMID 30664679.
