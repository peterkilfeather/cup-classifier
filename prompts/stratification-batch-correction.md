# Grilling Prompt — Stratification and Batch Correction: What Can They Achieve?

**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#4 Stratification and Batch Correction — What Can They Achieve?](https://github.com/peterkilfeather/cup-classifier/issues/4)
**Depends on:** #2 Confound Diagnostic Protocol (resolved) — protocol at `docs/confound-diagnostic-protocol.md`

## Concrete inputs from #2 (Confound Diagnostic Protocol)

The diagnostic protocol is settled. Before this session can decide what stratification/batch correction can achieve, it needs to know what the protocol actually says. Here are the key outputs:

### 4 methods × 6 modalities

Every method runs on every modality (probe-averaged methylation 148, per-CpG methylation 32K→PCs, FEM4 256, fragment length 369, CNVkit 38×3 thresholds, end density 31K→PCs):

1. **PCA + confound coloring** — PC1×PC2 colored by Tissue vs BCT vs Source vs Sex. Side-by-side comparison: if BCT separation is cleaner than Tissue, that's the 2-second money slide.
2. **PERMANOVA Type III (marginal)** — `dist_matrix ~ Tissue + Source + BCT + Sex`, each term's R² conditional on all others. Entry-order sensitivity to bracket shared variance (required because BCT~Tissue Cramer's V = 0.756). Year excluded (collinear with Source).
3. **Classifier negative controls** — penalized logistic regression or RF predicting Tissue (positive control), BCT (key test), Source (secondary), shuffled labels (null). Macro-F1 with chance baseline.
4. **Per-feature variance partitioning** — `feature ~ Tissue + Source + BCT + Sex` linear model per feature, decompose via ANOVA SS. Distribution of variance fraction explained by each covariate.

### Dual thresholds for confound dominance

A confound must satisfy BOTH criteria:

| PERMANOVA (conditional R²) | Classifier (F1) | Verdict |
|---|---|---|
| ≥ 0.10 | ≥ 0.8 × Tissue-F1 | **Dominant confound** — must mitigate |
| ≥ 0.10 | < 0.8 × Tissue-F1 | Material — design around it |
| < 0.10 | ≥ 0.8 × Tissue-F1 | Borderline — investigate (feature selection artifact?) |
| < 0.10 | < 2× chance | Weak — note and proceed |

### Scope: Full dataset (164) + EDTA sensitivity (96)

EDTA-only subset holds BCT constant. Interpretation matrix:

| Full | EDTA | Inference |
|------|------|-----------|
| BCT dominates | Tissue strong | Real biology, BCT was noise |
| BCT dominates | Tissue collapses | BCT carried apparent tissue signal |
| Tissue dominates | Tissue strong | BCT is a minor confound |

Pancreas excluded from EDTA subset (all Citrate). Colon (9), liver (9), prostate (17), stomach (21), healthyblood (40).

---

## The hard problem this ticket must resolve

The diagnostic protocol will tell us *which* confounds dominate. But even before those results, the structural problem is clear:

**Colon and pancreas are single-source (Fox Chase).** A classifier can learn "Fox Chase" and correctly classify both tissues — for the wrong reason. No stratification or batch correction can mathematically separate Tissue from Source when they are perfectly correlated.

For multi-source tissues (prostate: Fox Chase + Sowalsky; liver: Fox Chase + Audubon; stomach: Fox Chase + Audubon; healthyblood: Audubon + NIH), the question is whether existing methods can produce trustworthy results.

## What needs deciding

### For multi-source tissues

- **Source-stratified CV** — train on one source, test on another (e.g., Fox Chase stomach → Audubon stomach). Is this sufficient for multi-source classes? What performance gap is acceptable?
- **Batch correction** — ComBat, limma::removeBatchEffect, harmony. Can they work when batch is correlated with biology (Fox Chase liver is Citrate, Audubon liver is EDTA — BCT and Source are entangled)? When does correction help vs remove biological signal?
- **Matching / weighting** — propensity score matching to balance Source within tissue. Feasible with these sample sizes (min class: 9-12)?
- **Multi-task / domain adaptation** — learn source-invariant representations. Worth the complexity?

### For single-source tissues (colon, pancreas — both Fox Chase)

- Are they **unusable** until incoming data arrives?
- Or can certain claims be made from within-source comparisons? (Colon vs pancreas, both Fox Chase: if a classifier can distinguish them, that is real biological signal — they share the same batch. If it *cannot*, that's evidence the features lack tissue-discriminative power.)
- What about within-Fox-Chase comparisons that include a multi-source tissue? (e.g., train on Fox Chase colon + Fox Chase stomach, test on Audubon stomach — does colon benefit from sharing the batch with stomach?)

### Interaction with #2's diagnostic outputs

Once the protocol is executed, the BCT-F1 vs Tissue-F1 comparison directly informs batch correction strategy:
- If BCT-F1 ≪ Tissue-F1: correction may be unnecessary
- If BCT-F1 ≈ Tissue-F1: correction risks removing biological signal
- If Source-F1 (Audubon vs Fox Chase vs Sowalsky) ≈ Tissue-F1: stratification alone is insufficient

## What this ticket does NOT decide

- The diagnostic protocol itself (#2, resolved)
- Whether to hold out incoming data as validation (#5, blocked by this ticket — that decision needs the answer to "can the current data support trustworthy results?")
- What analysis to run immediately (#6, blocked by #2 + #3)

## Key files

- `docs/confound-diagnostic-protocol.md` — the full protocol from #2
- `CONTEXT.md` — domain model, confound documentation
- `input/metadata/metadata_cleaned.csv` — cleaned sample table
- `docs/data-quality-provenance.md` — cleaning decisions
