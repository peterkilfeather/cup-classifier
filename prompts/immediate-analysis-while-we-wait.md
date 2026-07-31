# Grilling Prompt — Immediate Analysis While We Wait

**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#6 Immediate Analysis While We Wait](https://github.com/peterkilfeather/cup-classifier/issues/6)
**Depends on:** #2 Confound Diagnostic Protocol, #3 Multi-Source Tissue Validation Protocol (both resolved)

## Context

This is the last ticket on the map. All upstream decisions are resolved:

| # | Decision | Key output |
|---|----------|-----------|
| 7 | Clean data | 164 samples, 6 tissues, `metadata_cleaned.csv` |
| 2 | Confound diagnostic protocol | 4 methods × 6 modalities, dual thresholds |
| 4 | BCT not dominant, Source is the confound | Source-stratified CV sufficient |
| 5 | 3-phase ablation (sequence it) | Phase 1 (now) → Phase 2a (incoming validation) → Phase 2b (integrate + cross-source FS) |
| 3 | Cross-source FS method designed | Per-tissue decomposition: coefficient comparison (stomach/liver), Cohen's d (healthyblood/prostate), deferred (colon/pancreas) |

### What Phase 1 looks like (already decided)

- **Classifier**: L1-penalized logistic regression (scikit-learn, saga solver). No RFECV — redundant on top of L1.
- **CV**: Source-stratified CV via StratifiedKFold stratified by Tissue, with post-hoc source coverage verification.
- **Best modalities** (from #4): FEM4 > probe-averaged methylation > fragment length.
- **BCT**: Primary on all 164 samples; EDTA-only sensitivity (96 samples) as positive control.
- **Success criteria for Phase 2a**: frozen model macro-F1 ≥ 0.7 × Phase 1 CV macro-F1 AND > 0.17 (chance).
- **Hyperparameters**: pinned in Phase 1 for Phase 2a reproducibility.

### What the original proposal identified

Three work streams:
1. **Confound diagnostics** — PCA, PERMANOVA, classifiers, variance partitioning — **ALREADY EXECUTED** (#2 commit 5b4aaee)
2. **EDTA-only subset sensitivity** — same diagnostics in EDTA subset — **ALREADY EXECUTED**
3. **Pipeline scaffolding** — data loading, CV infrastructure, metric tracking — partially built (scripts exist for diagnostics, but not the modelling pipeline)

## What this ticket decides

### 1. Phase 1 implementation plan — concrete steps

Phase 1 is defined at a high level but needs a detailed implementation plan:
- Data loading and merging (inner join metadata with feature files across all modalities)
- Pipeline code: L1-logreg with Source-stratified CV
- Per-modality evaluation: run FEM4, methylation, fragment length separately and in combination
- Metric tracking: macro-F1 with per-source and per-tissue breakdowns
- EDTA-only sensitivity: run the same pipeline on the 96-sample EDTA subset
- Hyperparameter pinning protocol for Phase 2a reproducibility

What gets built now vs what gets deferred?

### 2. Priority order

Given ~weeks until incoming data, what order should the team work through:
- Single modality first (FEM4 — best performer) then add others?
- Build full multi-modal pipeline from day one?
- Infrastructure (CV, metrics, logging) before any modelling?
- Confound diagnostic code reuse — the `run_diagnostic_parallel.py` infrastructure exists; can the modelling pipeline inherit from it?

### 3. Cross-source feature selection implementation

The method is designed (per-tissue decomposition) but not implemented:
- Stomach/liver: bootstrapped logreg coefficient comparison
- Healthyblood/prostate: Cohen's d distribution comparison
- Colon/pancreas: deferred to Phase 2b

Does implementation start now (so Phase 2b is ready when data arrives), or is it deferred entirely?

### 4. What does "done for Phase 1" mean?

Clear success criteria for Phase 1 completion:
- Pipeline runs end-to-end on all modalities?
- Reproducible outputs in `output/`?
- Hyperparameters and random seeds pinned and documented?
- What deliverables? (figure grid, CSV results table, trained model serialization?)

### 5. Relationship to downstream fog

The original handoff identified downstream decisions (feature modality prioritization, model selection, feature importance, missingness, enriched vs unenriched). Does Phase 1 resolve any of these? For example:
- Running Phase 1 on all modalities will rank them empirically — does that close "feature modality prioritization"?
- Imputation is needed for per-CpG analysis — does Phase 1 define the strategy?
- Enriched vs unenriched methylation: both exist in the data. Does Phase 1 test both?

## What this ticket does NOT decide

- The confound mitigation strategy (that's the whole map — and it's largely resolved)
- The cross-source FS method design (#3, resolved)
- Whether to hold out or integrate (#5, resolved)
- The protocol for testing multi-source generalization (#3, resolved)

## Key files

- `docs/multi-source-validation-protocol.md` — Phase 1 protocol with success criteria
- `docs/confound-diagnostic-protocol.md` — diagnostic protocol
- `docs/data-quality-provenance.md` — cleaning decisions
- `CONTEXT.md` — domain model
- `input/metadata/metadata_cleaned.csv` — cleaned metadata
- `scripts/run_diagnostic_parallel.py` — existing parallel pipeline infrastructure
- `output/diagnostic-protocol/adjudication_*.csv` — modality rankings
