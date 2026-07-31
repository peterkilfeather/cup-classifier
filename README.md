# Cup Classifier — Confound Diagnostic Protocol

**Issue:** [#2 Confound Diagnostic Protocol](https://github.com/peterkilfeather/cup-classifier/issues/2)
**Unblocks:** [#4 Stratification and Batch Correction](https://github.com/peterkilfeather/cup-classifier/issues/4)

Quantifies batch/confound effects across 6 feature modalities for the
cfDNA tissue-of-origin classifier.

## Results

Full adjudication tables and 85 figures in `output/diagnostic-protocol/`.

**Core finding:** BCT (blood collection tube) is *not* a dominant confound.
Unique BCT variance never exceeds R²=0.021 across any modality. The apparent
classifier ability to predict BCT is aliasing through Tissue (Cramer's V=0.756).
Source is the actionable confound (R²=0.01-0.05). Source-stratified CV is the
recommended mitigation.

| Modality | Tissue R² | BCT R² | Source R² | Tissue F1 | BCT F1 |
|---|---|---|---|---|---|
| Methylation (probe-avg) | 0.131 | 0.005 | 0.029 | 0.494 | 0.517 |
| FEM4 (256) | 0.150 | 0.018 | 0.047 | 0.705 | 0.734 |
| Fragment length (369) | 0.141 | 0.004 | 0.052 | 0.459 | 0.511 |
| CNVkit thr0.10 | 0.054 | 0.011 | 0.017 | 0.463 | 0.422 |
| Per-CpG (32K→5 PCs) | 0.068 | 0.015 | 0.025 | 0.349 | 0.294 |
| End density (31K→5 PCs) | 0.062 | 0.021 | 0.024 | 0.357 | 0.476 |

## How to re-run

```bash
# Prerequisites: Python 3.12+ with numpy, pandas, scikit-learn, statsmodels,
# matplotlib, seaborn, and R 4.x with vegan package.

# 1. Data cleaning (skip if input/metadata/metadata_cleaned.csv exists)
python3 scripts/clean_metadata.py
python3 scripts/clean_features.py

# 2. Full diagnostic protocol (4 methods × 6 modalities × 2 scopes, ~5 min)
python3 scripts/run_diagnostic_parallel.py

# 3. Entry-order sensitivity table
python3 scripts/postprocess_permanova.py
```

Output goes to `output/diagnostic-protocol/`.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/clean_metadata.py` | Raw metadata → 164 clean samples |
| `scripts/clean_features.py` | Drop all-NaN features per modality |
| `scripts/run_diagnostic_parallel.py` | Main pipeline: PCA, PERMANOVA, classifier, variance partitioning |
| `scripts/postprocess_permanova.py` | PERMANOVA entry-order sensitivity table + plot |

## Methods

Four converging methods applied to each modality:

1. **PCA + confound coloring** — PC1×PC2 colored by Tissue/BCT/Source/Sex
2. **PERMANOVA Type III (marginal)** — `dist ~ Tissue + Source + BCT + Sex`
   with entry-order sensitivity (BCT-first, sequential)
3. **Classifier negative controls** — L1-penalized logistic regression with
   nested CV predicting Tissue/BCT/Source/shuffled labels
4. **Per-feature variance partitioning** — OLS per feature, ANOVA Type II SS,
   grouped by covariate

High-dim modalities (per-CpG, end density) are PCA-reduced (K=5 minimum via
elbow) before PERMANOVA and classifier; variance partitioning uses raw features
with 5000-feature subsampling for speed.
