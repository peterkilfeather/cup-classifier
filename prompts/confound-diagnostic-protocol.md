# Grilling Prompt — Confound Diagnostic Protocol

**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#2 Confound Diagnostic Protocol](https://github.com/peterkilfeather/cup-classifier/issues/2)

## Context

The dataset has been cleaned — 164 samples across 6 tissues. Full provenance at `docs/data-quality-provenance.md`. Clean metadata at `input/metadata/metadata_cleaned.csv`. The key confounds (Source~Tissue, BCT~Tissue, Year~Tissue, Sex) are documented in `CONTEXT.md` and remain unchanged post-cleaning.

The core hypothesis is that **BCT (blood collection tube) is the dominant batch effect**, but this needs empirical proof. The question this ticket resolves: what protocol to use for that proof.

## What needs deciding

- **Feature modalities to test** — all (methylation, FEM4, fragment length, CNV, end density) or a representative subset?
- **Diagnostic methods** — what combination of PCA (colored by confound variables), PERMANOVA, variance partitioning, classifier-on-source-label negative controls, and other approaches?
- **Thresholds** — what quantitative criterion separates "confound dominates" from "biology dominates"?
- **Scope** — full dataset (164 samples, mixed BCTs) or within-EDTA subset as sensitivity?

## Data summary for the session

| Tissue | Count | Sources | BCTs |
|--------|-------|---------|------|
| healthyblood | 40 | Audubon (27), NIH (13) | EDTA (40) |
| prostate | 38 | Fox Chase (18), Sowalsky (20) | EDTA (18), Streck (20) |
| liver | 24 | Fox Chase (12), Audubon (12) | Citrate (12), EDTA (12) |
| colon | 22 | Fox Chase (22) | ACD (13), EDTA (9) |
| stomach | 21 | Fox Chase (11), Audubon (10) | EDTA (21) |
| pancreas | 19 | Fox Chase (19) | Citrate (19) |

## What this session does NOT decide

- How to fix the confounds once diagnosed (that's #4: Stratification and Batch Correction)
- The multi-source validation protocol (that's #3, a parallel unblocked ticket)
- What analysis to run immediately (#6, blocked by this one and #3)

## Files to load

- `CONTEXT.md` — domain model, confound documentation
- `input/metadata/metadata_cleaned.csv` — the cleaned sample table
- `docs/data-quality-provenance.md` — cleaning decisions
