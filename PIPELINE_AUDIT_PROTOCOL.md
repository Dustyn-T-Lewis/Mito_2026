# Pipeline Stage Audit & Tightening Protocol

A gated, **one-stage-at-a-time** protocol for auditing, validating, and tightening the
A_Mito_2026 proteomics pipeline (rat H9c2 mitochondrial-transplant DIA-MS). Each stage is
a self-contained pass that ends at a human gate; context is cleared between stages so the
fresh session re-reads this file and continues with no memory of the previous one.

---

## How to drive this

Run stages **in order** (below). For each stage:

1. Start a fresh session (`/clear`).
2. Paste the launcher:
   > Read `PIPELINE_AUDIT_PROTOCOL.md` and execute it for **STAGE = `<id>`**. Work only on
   > that stage. Do the read/validate pass first, show me findings, and stop at the GATE
   > before changing anything.
3. Approve / adjust at the gate. Let it tighten + verify + update memory + clean outputs.
4. When it reports the stage complete, `/clear` and move to the next stage.

**Stage order (do not skip ahead — later stages consume earlier outputs):**

| id | scope |
|----|-------|
| `00-01` | `00_input/` + `01_Filtering/` |
| `02`    | `02_Normalization/` (+ `imputation/`) |
| `03`    | `03_DEP/a_non_imputed/` (reported) + `b_imputed/` (exploratory) |
| `04-shared` | `04_Figures/shared/` engine (config, style, builders, pathway/mitocarta utils, comparison_panels) |
| `F01` … `F06` | one figure per stage, in number order |

A stage may read earlier outputs but must **never edit files outside its own scope**.

---

## Per-stage gated loop

Turn this into a todo list at the start of every stage.

### 1. Map (read-only)
- Read every `.R` and `.qmd` in the stage. State the **input→output contract**: each read
  points to a file a prior stage writes; each write is consumed downstream. Flag orphans,
  stale/broken paths, `setwd`, hardcoded absolute paths.
- In one paragraph, state the **scientific intent** of the stage and of each script — what
  question it answers and why it exists. If intent is unclear from the code, that is itself
  a finding (the code must make its purpose self-evident).

### 2. Validate against the literature
- For every analytical choice (filter rule, normalization, imputation mechanism, model,
  contrast, significance rule, enrichment, dedup), confirm it against the literature and
  cite PMID/DOI (anchors below). Mark each `AFFIRMED` / `CONCERN` / `FLAG`.
- Confirm package **API correctness** against installed versions (`?fn`, vignette) — right
  function, right arguments, right defaults. No hand-rolled reimplementations of functions
  that already exist in proteoDA / limma / fgsea / etc.
- Make **pathways and methods explicit**: where a non-obvious method or threshold is used,
  the code should carry a one-line *why* with a citation at that line — not a wall of text.

### 3. GATE — present findings, wait for approval
- Severity-graded table: `Severity | file:line | issue | fix | citation`.
- **Do not change anything before approval.** Statistical thresholds, model formulas,
  contrast definitions, normalization/imputation/DE methods are **never** changed without
  explicit sign-off — document them, don't alter them.

### 4. Tighten (after approval) — minimal, parsimonious, human
Apply the standards below. The goal is **fewer, clearer lines**, not more. If a change grows
total line count by >10%, you are over-engineering — revert and rethink.

### 5. Verify — no silent regression
- `Rscript -e 'invisible(parse("<file>"))'` every edited file.
- Re-run the stage's compute scripts. Re-snapshot outputs (`md5`) and **diff against the
  pre-edit snapshot**. Comment-only / rename-only edits must reproduce **byte-identical**
  outputs. Any numeric change must be explained (and approved) — a missing output file is a
  hard stop; revert.
- Stochastic steps must be deterministic on re-run (seed fixed). Prove it where practical
  (build twice, compare md5).

### 6. Update memory
- Write/refresh one memory file per durable fact this stage established (intent, method +
  citation, a non-obvious decision, the input→output contract). Follow the project memory
  format; add a one-line pointer to `MEMORY.md`. Update existing entries rather than
  duplicating; delete entries proven wrong. Do **not** record what the code already says.

### 7. Clean outputs
- Remove stale orphan outputs this stage produced that no current script generates (figures
  with no producer, CSVs with no consumer, caches from deleted builders). `git rm` tracked
  files; `rm` untracked. **Never** delete a live output a downstream stage reads. List every
  deletion with its justification.

### 8. Close the stage
- Summarize: findings fixed, outputs unchanged vs changed (with reason), memory updated,
  files removed, what still needs regeneration downstream.
- Commit on the working branch (one stage = one commit; message per the repo's git rules).
- Tell me to `/clear` and launch the next stage.

---

## Standards

### R best practices
- `here::here()` for all paths; never `setwd()`.
- `set.seed()` immediately before any stochastic step; make it order-independent where a
  loop calls the RNG repeatedly.
- Separate **compute from plotting**: builders cache; panels read caches and plot. A figure
  script must not re-run the DE or re-derive reported statistics — it reads the DEP tables.
- Tidyverse idiom, but **import only what's used** — no blanket `library(tidyverse)` in a
  shared file sourced by many scripts. Prefer the **canonical** package for the task
  (proteoDA, limma, fgsea, GSVA, singscore, WGCNA); do not introduce exotic or redundant
  dependencies when a canonical one already does the job.
- One clear function for a job, defined once and reused — not three near-duplicate copies.
  But do **not** invent new abstractions/wrappers for single-use code; inline it.

### Anti-AI-artifact (remove on sight)
- Banner comment lines (`# =====`, `# -----`, decorative section dividers).
- Comments that restate the next line (`# calculate mean` above `mean(x)`); keep *why*, drop
  *what*. One line max where a comment earns its place.
- Defensive bloat: `tryCatch` / `suppressWarnings` / `suppressMessages` around code that
  fails or warns informatively; remove so real diagnostics surface.
- `message()` narrating every step ("processing contrast 1 of 5…").
- Emoji, and "AI voice" prose in comments or quarto.
- **Naming ladders** (`df`, `df2`, `df_final`, `temp_`, `result_`, `_clean`): use one
  domain-meaningful name. **Do not rename** variables, files, or columns that are already
  clear — gratuitous renaming churns diffs and breaks intuition. Rename only when the
  current name is actively misleading, and say why.
- Config that earns its keep stays shared (used 3+ places, or a path that changes per
  machine); single-use constants get inlined.

### Quarto cleanliness
- The `.qmd` must mirror the production `.R` logic and **represent the actual workflow** —
  no drift, no double-applied transforms, no diagnostics computed on the wrong scale.
- Narrate the *why* in flowing prose between chunks; let the code show the *what*.
- It must render cleanly start to finish. Inputs it reads must exist. No dead chunks
  referencing intermediates the pipeline no longer produces.

---

## Literature anchors by stage

| stage | check | citation |
|-------|-------|----------|
| 01 | group-wise missingness filter | Välikangas 2018 *Brief Bioinform* PMID:27694351 |
| 01 | proteoDA DAList workflow | Thurman 2023 *JOSS* 8:5184 |
| 02 | cyclic loess for DIA-MS | Ritchie 2015 *NAR* PMID:25605792; Arend 2025 PMID:40231271 |
| 02 | log2 before normalization | Välikangas 2018 PMID:27694351 |
| 02 | MAR/MNAR-aware imputation | Lazar 2016 *JPR* PMID:26906401 |
| 02 | imputed data excluded from reported DE | Webb-Robertson 2015 *JPR* PMID:25855118 |
| 03 | duplicateCorrelation + block | limma User's Guide §9.7; Smyth 2005 |
| 03 | robust eBayes | Phipson 2016 PMID:27008012 |
| 03 | 2×2 factorial / interaction contrasts | Law 2020 *F1000Res* PMID:32346458 |
| 03 | Xiao π-score (`P^|logFC|`) | Xiao 2014 PMID:24478644 |
| 03 | BH-FDR | Benjamini & Hochberg 1995 *JRSS-B* |
| 04 | fgsea, moderated-t ranking | Korotkevich 2021 *bioRxiv* 10.1101/060012 |
| 04 | ORA universe = all tested proteins | Timmons 2015 PMID:26450340; Reimand 2019 PMID:30664679 |
| 04 | EnrichmentMap Jaccard/overlap dedup | Merico 2010 PMID:21085593; Reimand 2019 PMID:30664679 |
| 04 | WGCNA soft-power / preservation | Langfelder 2008 PMID:19114008; 2011 PMID:21901115 |
| 04 | GSVA / singscore | Hänzelmann 2013 PMID:23323831; Foroutan 2018 PMID:30400809 |

---

## Hard stops (refuse without explicit approval)
- Changing a statistical threshold, model formula, contrast, or DE/normalization/imputation
  method.
- Deleting any output a downstream stage still reads.
- Adding a package dependency without naming why it's canonical.
- Proceeding past the GATE before findings are approved, or past a missing/regressed output.
- Net-adding comments, or renaming things that were already clear.
