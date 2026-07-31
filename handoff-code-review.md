# Handoff: Cup Classifier — Confound Diagnostic Protocol

## What was done

Executed the Issue #2 Confound Diagnostic Protocol across 6 feature modalities and 2 scopes (Full 164 samples, EDTA 96 samples), producing adjudicated verdicts for BCT, Source, Sex, and Tissue confounds. The protocol runs 4 converging methods per modality: PCA, PERMANOVA (Type III marginal), classifier negative controls (L1 logistic regression with nested CV), and per-feature variance partitioning.

### Scripts (all under `scripts/`)

- **`run_diagnostic_parallel.py`** — main pipeline (12 modality×scope tasks via ProcessPoolExecutor, ~5 min)
- **`postprocess_permanova.py`** — entry-order sensitivity table and PERMANOVA summary plot
- **`clean_metadata.py`**, **`clean_features.py`** — data cleaning (from issue #7, already resolved)

### Output (`output/diagnostic-protocol/`)

- 85 PNG figures (12 scree + 48 PCA-colored + 12 classifier + 12 varpart + 1 PERMANOVA summary)
- 2 adjudication CSVs (`adjudication_Full_(164).csv`, `adjudication_EDTA_(96).csv`)
- Entry-order sensitivity table (via postprocess script)

### Fixes applied during review

Two bugs were identified and fixed during code review, then the pipeline was re-run to regenerate all outputs:

1. **Variance partitioning silently produced all zeros for EDTA scope.** The formula `feat ~ C(Tissue) + C(Source) + C(BCT) + C(Sex)` included `C(BCT)` with one level (EDTA), crashing `smf.ols` and recording zeros for *all* covariates. Fixed by dynamically filtering to active covariates (`nunique() > 1`), matching the PERMANOVA R script's treatment of single-level factors. EDTA VP_Tissue now correctly shows 0.10--0.22.

2. **Classifier chance baseline was a misleading analytical formula.** `chance_f1` computed F1 under uniform random prediction, which doesn't match class-weighted L1 logistic regression behavior. Replaced with the empirical `Shuffled_Tissue` F1 (already computed as a negative control). The redundant "Chance" horizontal line was removed from classifier plots -- the Shuffled_Tissue bar is the null baseline.

On output count: 5 orphaned PNGs from a previous run (enriched methylation with stale naming convention) were also cleaned up, bringing total from 90 to 85.

### Supporting files

- **`GLOSSARY.md`** — canonical terms for confound analysis
- **`MISSION.md`** — learning mission (understanding why BCT is not dominant)
- **`RESOURCES.md`** — curated sources on PERMANOVA, batch effects
- **`NOTES.md`** — teaching session log
- **`learning-records/`** — 2 records of concepts learned/confirmed
- **`lessons/`** — 3 HTML lessons (BCT aliasing, variance partitioning, Source as real confound)
- **`assets/lesson.css`** — shared stylesheet for lessons

## Core findings

**BCT is NOT a dominant confound.** Unique marginal R² never exceeds 0.04 (end density). Classifier ability to predict BCT (F1 up to 0.73) is aliasing through Tissue correlation (Cramer's V structure: ACD→colon, Streck→prostate, Citrate→liver+pancreas). Entry-order sensitivity proves it: BCT jumps from R²=0.018 to 0.155 when entered first (FEM4). EDTA subset confirms: tissue signal strengthens when BCT is removed.

**Source is the actionable confound.** Marginal R²=0.014-0.052, persists in EDTA (independent of tube type), classifier F1 often exceeds Tissue F1. Source-stratified CV (train on one source, test on another) is the correct mitigation. Batch correction would harm — Source is aliased with Tissue.

**Best modalities:** FEM4 (Tissue R²=0.150, F1=0.705) > probe-averaged methylation (0.131, 0.494) > fragment length (0.141, 0.459).

## What the code review should focus on

The main deliverable is `scripts/run_diagnostic_parallel.py`. Key things to check:

1. **PERMANOVA**: R `vegan::adonis2` call with `by="margin"` — single-level factors (BCT in EDTA) dynamically dropped from formula. Entry-order sensitivity runs three orderings.
2. **Classifier**: L1-penalized logistic regression with nested CV. 5 outer folds, 3 inner folds (GridSearchCV over 6 C values). Stratified folds. Class-balanced. Falls back to train/test split when classes are too small.
3. **Variance partitioning**: `smf.ols` formula API for proper categorical encoding. ANOVA Type II SS. Covariates grouped by name (summing across dummy columns).
4. **PCA for high-dim**: `elbow_k()` with min_k=5 floor. Randomized SVD. NaN imputation before scaling.
5. **Parallelism**: `ProcessPoolExecutor(max_workers=4)` across 12 tasks. Each worker loads its own data files independently.
6. **EDTA subset**: Filters to BCT=='EDTA' for the scope. PERMANOVA auto-drops single-level BCT factor.
7. **Adjudication**: Dual threshold logic (R² ≥ 0.10 AND F1 ≥ 0.8×Tissue-F1 for "Dominant"). Prints table and saves CSV.

## How to run the review

```bash
# In a fresh session referencing this file:
# 1. Check syntax
python3 -m py_compile scripts/run_diagnostic_parallel.py

# 2. Run a quick smoke test (one modality, one scope, 99 permutations)
#    (modify the R script to use permutations=99 for speed)
python3 scripts/run_diagnostic_parallel.py

# 3. Check output integrity
ls output/diagnostic-protocol/*.png | wc -l
cat output/diagnostic-protocol/adjudication_Full_\(164\).csv
```

## Previous session context (tl;dr)

- Issue #4 (Stratification and Batch Correction) was **closed** with findings comment
- The user spent 3 months believing BCT was the dominant confound — the protocol corrected this
- Three teaching lessons were created to explain the findings
- The user wants the work reviewed before it's considered complete
