# Grilling Prompt — Multi-Source Tissue Validation Protocol

**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#3 Multi-Source Tissue Validation Protocol](https://github.com/peterkilfeather/cup-classifier/issues/3)
**Updated by:** #5 (Incoming EDTA Data — Hold Out or Integrate?, resolved)

## Context

This ticket started as a narrow question (protocol for testing cross-source generalization in stomach, healthyblood, liver). After #5 resolved with a 3-phase ablation plan, the scope expanded significantly — the cross-source invariant feature selection method that #5's Phase 2b depends on needs to be designed here first.

### What #4 and #5 established

- **BCT is NOT the dominant confound** (max marginal R²=0.04). Source is the actionable confound.
- **Source-stratified CV** is the right mitigation for multi-source tissues.
- **Phase 1** (current 164 samples): develop baseline model.
- **Phase 2a** (incoming 41 samples): evaluate frozen model as external validation.
- **Phase 2b** (integrate ~205 samples): retrain with **cross-source invariant feature selection**.
- Best modalities: FEM4 > probe-averaged methylation > fragment length.

### The design problem

The core of Phase 2b is a feature selection method that selects features whose tissue association is **consistent across sources**. The canonical approach: for each feature, fit a model per source and select features where the tissue-effect coefficients are stable. But the devil is in the detail — especially with small per-source sample sizes.

## What needs deciding

### 1. Cross-source invariant feature selection method

The primary design question. Options to grill through:

- **Bootstrapped logreg per source** — for each feature, fit `feature ~ tissue` within each source (Fox Chase, Audubon, Sowalsky, NIH). Select features where tissue coefficients have the same sign and overlapping confidence intervals across sources. Simple, interpretable. But some sources have very small N per tissue (e.g., Audubon stomach = 10).
- **Interaction test in pooled model** — fit `feature ~ source * tissue`. Select features where the interaction term (source × tissue) is non-significant × where the tissue main effect is significant. Statistically principled, but ANOVA assumptions may not hold for all modalities.
- **Domain-invariant feature learning** — adversarial or regularized approaches (e.g., CORAL, DANN). More complex, may overfit on 164 samples.
- **Ensemble of source-specific classifiers** — train a separate classifier per source, select features important across most source-specific models. Practical, but doesn't directly address consistency.
- **Rank-by-stability** — rank features by tissue association in each source, then select high-rank features that appear in the top-K across multiple sources. Non-parametric, robust to distributional assumptions.

For each candidate: what are the N requirements per source per tissue, and which tissues/sources qualify?

### 2. Phase 2b protocol

- How does cross-source feature selection interact with Source-stratified CV?
- Should feature selection happen inside the CV loop (nested) or outside?
- Metrics: macro-F1? Per-source accuracy? Calibration?
- Success criteria: what delta between Phase 2a (frozen model) and Phase 2b (retrained) constitutes a meaningful improvement?

### 3. Phase 1 relationship

Does the original "test NOW" protocol run:
- **A)** Independently as early signal — run now on current 164 samples using the existing best approach (e.g., source-stratified RF or logreg)
- **B)** Wait for the cross-source feature selection design — Phase 1 becomes a placeholder that inherits the method from this ticket

### 4. Shared metrics across phases

For the 3-phase ablation to be clean, the same evaluation criteria should apply:

| Phase | Train | Test | Question |
|-------|-------|------|----------|
| 1 | Current 164 | Source-stratified CV | Baseline performance |
| 2a | Phase 1 model (frozen) | Incoming 41 | Does performance generalize? |
| 2b | All ~205 (cross-source FS) | Source-stratified CV | Does cross-source selection improve generalization? |

- What is the primary metric? (macro-F1? balanced accuracy? per-source accuracy?)
- What threshold separates "generalizes" from "does not generalize"?
- How do we handle class imbalance across sources?

### 5. Data available for immediate testing

Current multi-source tissues for early validation:

| Tissue | Source A (N) | Source B (N) | BCTs |
|--------|-------------|-------------|------|
| healthyblood | Audubon (27) | NIH (13) | Both EDTA |
| stomach | Fox Chase (11) | Audubon (10) | Both EDTA |
| prostate | Fox Chase (18) | Sowalsky (20) | EDTA vs Streck |
| liver | Fox Chase (12) | Audubon (12) | Citrate vs EDTA |

For stomach and healthyblood: BCT is held constant — cleanest test of cross-source biology.
For liver: BCT changes with Source — confounded, but informative as a stress test.
For prostate: BCT changes with Source — same confound structure as liver.

## What this ticket does NOT decide

- Which ML architecture to use (downstream fog)
- Feature importance methodology for biological interpretation (downstream fog)
- The actual execution of Phase 1/2a/2b (that's #6, blocked by this ticket — the protocol needs to be designed first)

## Key files

- `docs/confound-diagnostic-protocol.md` — diagnostic protocol from #2
- `output/diagnostic-protocol/adjudication_*.csv` — confound adjudication results showing Source is the actionable confound
- `CONTEXT.md` — domain model
- `input/metadata/metadata_cleaned.csv` — cleaned sample table
- `docs/data-quality-provenance.md` — cleaning decisions
