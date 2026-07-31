# Multi-Source Tissue Validation Protocol

**Status:** Resolved via grilling session (2026-07-31)
**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#3 Multi-Source Tissue Validation Protocol](https://github.com/peterkilfeather/cup-classifier/issues/3)
**Updated by:** #5 (Incoming EDTA Data — Hold Out or Integrate?, resolved)
**Blocks:** #6 (Immediate Analysis While We Wait)

---

## Core Question

How should we select features whose tissue-of-origin association is **consistent across sources**, and how does this fit into the 3-phase ablation plan?

---

## Structural constraints (from cleaned metadata)

Only 2 of 4 sources have >1 tissue class, which limits which feature selection methods are viable.

| Source | Tissues | N | Multi-class? |
|--------|---------|---|-------------|
| Fox Chase | colon, liver, pancreas, prostate, stomach | 85 | Yes — 5 classes |
| Audubon | healthyblood, liver, stomach | 46 | Yes — 3 classes |
| Sowalsky | prostate | 20 | No |
| NIH Clinical Center | healthyblood | 13 | No |

Shared tissues between multi-class sources (enabling coefficient comparison):

| Tissue | Source A | Source B | Both multi-class? | BCT |
|--------|----------|----------|-------------------|-----|
| stomach | FC (11, EDTA) | Audubon (10, EDTA) | Yes | Clean |
| liver | FC (15, Citrate) | Audubon (9, EDTA) | Yes | Confounded |
| healthyblood | Audubon (27) | NIH (13) | No — NIH single-class | Clean |
| prostate | FC (18, EDTA) | Sowalsky (20, Streck) | No — Sowalsky single-class | Confounded |
| colon | FC (22) | — | Single-source in Phase 1 | — |
| pancreas | FC (19) | — | Single-source in Phase 1 | — |

Colon and pancreas gain cross-source representation in Phase 2b when incoming data (likely new sources) arrives.

---

## 1. Cross-source invariant feature selection — per-tissue decomposition

No single method applies cleanly to all tissues because the source structure differs per tissue. Instead, use a **per-tissue approach** where the consistency check depends on what each tissue's data supports.

### For tissues shared between multi-class sources (stomach, liver)

**Method: Bootstrapped logreg coefficient comparison**

For each feature independently:
1. In each source (FC, Audubon), fit `feature ~ tissue` as a one-vs-rest univariate logistic regression
2. Get the tissue coefficient $\beta$ and its 95% CI (via bootstrap or closed-form)
3. **Pass** if: coefficients have the same sign AND CIs overlap AND both |β| > 0 (non-zero effect)
4. **Fail** if: opposite signs, or non-overlapping CIs, or one source has β ≈ 0 while the other has |β| > 0.5

**Liver caveat**: BCT is confounded with source (FC=Citrate, Audubon=EDTA). Features that fail the consistency check could fail due to BCT artifact rather than real biological difference. Run the same check on the BCT-controlled EDTA subset (stomach, healthyblood) as a positive control — if a feature passes in stomach but fails in liver, BCT is the likely cause.

**Per-tissue N adequacy**:
- Stomach: FC N=11, Audubon N=10. Univariate logreg with ~10/11 positives vs ~75/36 negatives per source. CIs will be wide → lenient CI overlap criterion.
- Liver: FC N=15, Audubon N=9. Audubon liver N=9 is marginal — bootstrap CIs will be very wide. Most features will "pass" by default. Use coefficient sign + effect size (|β| > 0.2) as the primary criterion, not CI overlap.

### For tissues with single-class sources (healthyblood, prostate)

**Method: Two-sample distribution comparison**

Cannot estimate a tissue coefficient (the source has only one tissue). Instead, check whether the feature's **marginal distribution** is similar across sources.

For each feature:
1. Compute Cohen's d between the two sources' samples of the same tissue
2. **Pass** if |d| < 0.5 (small or negligible difference)
3. **Fail** if |d| ≥ 0.8 (large difference → source-specific artifact)

**What this tests**: Distribution invariance, not tissue-effect consistency. A feature can pass this check (similar values in both sources) but have no tissue-discriminative power. The tissue-discriminative power is separately assessed via the tissue-classification model.

**Healthyblood**: Clean test — BCT constant (both EDTA), sources are Audubon (N=27) vs NIH (N=13). |d| < 0.5 with N=13+27 gives adequate power.

**Prostate**: BCT confounded (FC=EDTA, Sowalsky=Streck). Same Cohen's d check, but flag any feature where |d| ≥ 0.8 as "BCT-susceptible." Cross-reference with the EDTA-subset positive controls.

### For single-source tissues in Phase 1 (colon, pancreas)

No cross-source consistency check possible. These tissues are excluded from the invariant feature selection until Phase 2b when incoming data adds new sources.

### Assembly: combining per-tissue checks into a feature set

Each feature gets a **per-tissue passport** — a list of which tissues it was checked for and whether it passed.

**Selection strategy** (Phase 2b, run inside CV loop):
- **Core invariant set**: Features that pass ALL applicable checks. Highest cross-source confidence.
- **Extended set**: Features that pass checks for ≥50% of applicable tissues. Larger, some risk.
- **Flagged set**: Features that fail one or more checks. Excluded from invariant selection but may still be informative for single-source tissues.

For Phase 2b model training, use the **core invariant set** as the primary feature set. If performance is poor, fall back to the extended set.

---

## 2. Phase 1 protocol — baseline on current 164 samples

### Design decisions (from grilling session)

| Question | Decision |
|----------|----------|
| Feature selection method | L1-penalized logistic regression (implicit feature selection via L1 penalty). No RFECV — redundant on top of L1. |
| Cross-validation | StratifiedKFold(n_splits=5) stratified by Tissue, with post-hoc source coverage verification — each training fold must contain all 4 sources. |
| Modality | Start with FEM4 (256 features, best performer per #4). Probe-averaged methylation as secondary. |
| Relationship to cross-source FS | Independent. Phase 1 runs now without cross-source invariant constraint. The cross-source method is designed for Phase 2b. |

### Pipeline

```
For each modality (FEM4 primary, probe_meth secondary):
  1. Inner join feature matrix with cleaned metadata (164 samples)
  2. StandardScaler per fold
  3. Outer CV: StratifiedKFold(n_splits=5) stratified by Tissue
     - Verify each training fold contains all 4 sources
     - If source coverage fails, reject split and re-shuffle
  4. Per fold:
     a. Scale training data, apply same transform to test
     b. Fit LogisticRegression(
          penalty='l1',
          solver='saga',            # only saga supports L1 + multinomial
          multi_class='multinomial', # native multiclass, not one-vs-rest
          C=<tuned via inner GridSearchCV>,
          class_weight='balanced',
          max_iter=5000,
          random_state=42)
        C grid: np.logspace(-3, 1, 6) → [0.001, 0.01, 0.1, 1, 10]
        Inner CV: StratifiedKFold(n_splits=min(3, min_class_count))
     c. Evaluate on held-out fold
  5. Report: macro-F1 (primary), balanced accuracy (secondary),
     per-source accuracy (diagnostic), chance baseline

Phase 2a reproduction note: The exact C values differ per fold (tuned).
Freeze the PIPELINE (scaling + C-tuning + L1 fit), not any single C.
```

### Outputs

- Frozen model + selected feature set (non-zero L1 coefficients) for Phase 2a replication
- CV performance estimates with 95% CI
- Per-source accuracy breakdown

---

## 3. Phase 2a protocol — external validation on incoming ~41 samples

### Design decisions (from #5)

| Question | Decision |
|----------|----------|
| Model source | Phase 1 model, frozen — no retraining, no peeking |
| Evaluation | Apply frozen Phase 1 model to incoming samples |
| Feature set | Identical to Phase 1 (same pipeline, same selected features) |

### Success criterion

**Signal generalizes** if: macro-F1 on incoming 41 ≥ 0.7 × Phase 1 CV macro-F1 AND > 0.17 (6-class chance).

Per-tissue breakdown: which tissues/classes collapse under distribution shift? Report per-source-if-known (incoming sources may not be known until data arrives).

---

## 4. Phase 2b protocol — integrated ~205 samples with cross-source invariant feature selection

### Design decisions (from grilling session)

| Question | Decision |
|----------|----------|
| Feature selection | Per-tissue decomposition as described in §1. Run INSIDE the CV loop (nested). |
| Cross-validation | GroupKFold(n_splits=4) or LeaveOneGroupOut with Source as group. If source composition of incoming data is known, use that grouping. |
| Modality | Multi-modal: FEM4 + probe-averaged methylation + fragment length, with per-modality invariant feature sets combined. |
| BCT handling | Both routes — primary (all samples, flag BCT-susceptible features) + sensitivity (EDTA-only subset). |

### Success criterion

**Cross-source selection helps** if:
- Phase 2b macro-F1 ≥ Phase 2a macro-F1 + 1 standard error, OR
- Phase 2b closes ≥ 50% of the generalization gap (Phase 1 CV F1 − Phase 2a frozen F1)

If Phase 2a already generalized (≥ 0.7 × Phase 1 CV F1), then Phase 2b must improve absolute macro-F1 by ≥ 0.03.

---

## 5. Shared metrics across phases

| Phase | Train | Test | Primary metric | Secondary metrics |
|-------|-------|------|----------------|-------------------|
| 1 | Current 164 | Source-stratified CV | macro-F1 (95% CI) | Balanced accuracy, per-source accuracy, chance baseline |
| 2a | Phase 1 model (frozen) | Incoming 41 | macro-F1 on incoming | Per-tissue accuracy, per-source breakdown |
| 2b | All ~205 | Source-stratified CV | macro-F1 (95% CI) | Per-source accuracy, BCT-sensitivity comparison (Route B) |

All phases report:
- Macro-averaged F1 (primary — handles class imbalance)
- Balanced accuracy (secondary — interpretable per-class)
- Per-source accuracy (diagnostic — shows source-specific degradation)
- Chance baseline (uniform random classifier on the same class distribution)

---

## 6. BCT handling — both routes (liver and prostate)

### Route A (primary, all samples)

Run the cross-source invariant feature selection on all samples regardless of BCT. Features that pass consistency checks despite different BCTs are **stronger** evidence of biological signal. Features that fail are excluded from the invariant set (ambiguous — could be BCT artifact or real biological difference).

### Route B (sensitivity, EDTA-only subset, N=96)

Run the same cross-source checks on the EDTA-only subset where BCT is constant. This subset supports:
- healthyblood: Audubon (27) vs NIH (13) — both EDTA. Distribution comparison.
- stomach: Fox Chase (11) vs Audubon (10) — both EDTA. Coefficient comparison.
- prostate: FC (17, EDTA) only — single-source. No cross-source check.
- colon: FC (9, EDTA) only — single-source.
- liver: Audubon (9, EDTA) only — single-source.

**Interpretation**: Compare the invariant feature set from Route A vs Route B.

| Route A ∩ Route B | Interpretation |
|-------------------|---------------|
| ≥ 70% overlap | BCT is a minor confound for the invariant features. Consistent with #4 finding (BCT max marginal R² = 0.04). |
| < 50% overlap | BCT is materially contaminating the liver/prostate cross-source check. Report both sets; flag Route A features that fail Route B. |

---

## 7. Open implementation decisions (to be resolved in #6)

- **Incoming data source assignment**: Will incoming 41 samples have known source labels? Needed for Phase 2a per-source breakdown and Phase 2b grouping strategy.
- **Per-CpG and end-density modalities**: Not included in Phase 1 baseline. Should they be added in Phase 2b with PCA-first reduction?
- **Phase 2b CV strategy**: If incoming samples add new sources, `GroupKFold` with Source as group may give different fold sizes than Phase 1's `StratifiedKFold`. To be decided when source labels are known.
- **Frozen model serialization**: Format for saving Phase 1 model (pickle, joblib, ONNX) for Phase 2a loading.

---

## Key files

- `docs/confound-diagnostic-protocol.md` — diagnostic protocol from #2
- `output/diagnostic-protocol/adjudication_*.csv` — confound adjudication results showing Source is the actionable confound
- `CONTEXT.md` — domain model
- `input/metadata/metadata_cleaned.csv` — cleaned sample table
- `docs/data-quality-provenance.md` — cleaning decisions
- `scripts/run_diagnostic_parallel.py` — existing classifier infrastructure (L1-logreg, GridSearchCV)
