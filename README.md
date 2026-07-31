# Cup Classifier

Tissue-of-origin classifier from cfDNA methylation and fragmentomic features.
Cancer of Unknown Primary (CUP) application.

## Phase 1: Baseline Tissue Classifier

**Issue:** [#6 Immediate Analysis While We Wait](https://github.com/peterkilfeather/cup-classifier/issues/6)
**Blocks:** Phase 2a (external validation on incoming data)

L1-penalized logistic regression with Source-stratified CV on 164 current samples.

### Results

| Scope | Modality | macro-F1 (CV) | Features | Median C |
|-------|----------|---------------|----------|----------|
| Full (164) | **FEM4 (256)** | **0.705 ± 0.123** | 253 / 256 | 10.0 |
| Full (164) | Methylation (probe-avg) | 0.515 ± 0.086 | 147 / 148 | 10.0 |
| Full (164) | Fragment length (369) | 0.459 ± 0.102 | 199 / 368 | 1.58 |
| Full (164) | End density (31K→20 PCs) | 0.415 ± 0.100 | 20 PCs / 31K | 0.25 |
| Full (164) | Per-CpG methylation (32K→20 PCs) | 0.399 ± 0.120 | 20 PCs / 32K | 0.25 |
| EDTA (96) | **FEM4 (256)** | **0.715 ± 0.137** | 227 / 256 | 10.0 |
| EDTA (96) | Methylation (probe-avg) | 0.453 ± 0.089 | 145 / 148 | 10.0 |
| EDTA (96) | Fragment length (369) | 0.527 ± 0.105 | 177 / 368 | 1.58 |
| EDTA (96) | End density (31K→20 PCs) | 0.472 ± 0.090 | 20 PCs / 31K | 0.25 |
| EDTA (96) | Per-CpG methylation (32K→20 PCs) | 0.351 ± 0.030 | 20 PCs / 32K | 0.25 |

FEM4 is the strongest single modality (macro-F1 0.705). All 5 modalities well above
chance (0.167). EDTA sensitivity results consistent with full dataset — BCT is not
driving tissue classification. Results, CV tables, figures, and serialized models
in `output/phase1/`.

### Run

```bash
python3 scripts/run_phase1_pipeline.py
# Options: --modalities fem4 probe_meth, --skip-edta, --skip-plots
```

### Outputs

`output/phase1/`:
- `phase1_summary.csv` — combined results across all modalities/scopes
- `{modality}_cv_metrics.csv` — per-fold metrics (macro-F1, balanced accuracy, per-source accuracy)
- `{modality}_hyperparameters.json` — pinned hyperparams for Phase 2a reproduction
- `models/{modality}_full_model.joblib` — frozen model (scaler + L1-logreg) for Phase 2a
- `figures/{modality}_cv_results.png` — per-fold performance + per-source accuracy
- `figures/{modality}_coefficients.png` — L1 coefficient heatmap

---

## Confound Diagnostic Protocol

**Issue:** [#2 Confound Diagnostic Protocol](https://github.com/peterkilfeather/cup-classifier/issues/2)

Quantifies batch/confound effects across 6 feature modalities.

**Core finding:** BCT is *not* a dominant confound. Unique BCT variance never
exceeds R²=0.021 across any modality. Source is the actionable confound
(R²=0.01-0.05). Source-stratified CV is the recommended mitigation.

| Modality | Tissue R² | BCT R² | Source R² | Tissue F1 | BCT F1 |
|---|---|---|---|---|---|
| Methylation (probe-avg) | 0.131 | 0.005 | 0.029 | 0.494 | 0.517 |
| FEM4 (256) | 0.150 | 0.018 | 0.047 | 0.705 | 0.734 |
| Fragment length (369) | 0.141 | 0.004 | 0.052 | 0.459 | 0.511 |
| CNVkit thr0.10 | 0.054 | 0.011 | 0.017 | 0.463 | 0.422 |
| Per-CpG (32K→5 PCs) | 0.068 | 0.015 | 0.025 | 0.349 | 0.294 |
| End density (31K→5 PCs) | 0.062 | 0.021 | 0.024 | 0.357 | 0.476 |

### Run

```bash
# Data cleaning (skip if metadata_cleaned.csv exists)
python3 scripts/clean_metadata.py
python3 scripts/clean_features.py

# Diagnostic protocol (4 methods × 6 modalities × 2 scopes, ~5 min)
python3 scripts/run_diagnostic_parallel.py
python3 scripts/postprocess_permanova.py
```

Output goes to `output/diagnostic-protocol/`.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/data_loading.py` | Shared data loading (used by diagnostic + Phase 1) |
| `scripts/run_phase1_pipeline.py` | Phase 1: L1-logreg with Source-stratified CV |
| `scripts/run_diagnostic_parallel.py` | Confound diagnostics: PCA, PERMANOVA, classifier, variance partitioning |
| `scripts/postprocess_permanova.py` | PERMANOVA entry-order sensitivity |
| `scripts/clean_metadata.py` | Raw metadata → 164 clean samples |
| `scripts/clean_features.py` | Drop all-NaN features per modality |

## Methods (confound diagnostics)

Four converging methods applied per modality:
1. **PCA + confound coloring** — PC1×PC2 colored by Tissue/BCT/Source/Sex
2. **PERMANOVA Type III (marginal)** — `dist ~ Tissue + Source + BCT + Sex` with entry-order sensitivity
3. **Classifier negative controls** — L1-logreg predicting Tissue/BCT/Source/shuffled labels
4. **Per-feature variance partitioning** — OLS + ANOVA Type II SS

High-dim modalities (per-CpG, end density) are PCA-reduced (K=5 minimum via
elbow) before PERMANOVA and classifier; variance partitioning uses raw features
with 5000-feature subsampling for speed.
