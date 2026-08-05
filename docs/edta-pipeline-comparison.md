# EDTA Pipeline Comparison: `cup-classifier` vs `classifier`

> Generated 2026-07-31. Compares the EDTA-only classification pipelines between
> `/xscratch/farney/cup-classifier` (Phase 1) and `/xscratch/elnitski/farney/classifier`
> (the de novo BCT pipeline + EDTA followup).

## TL;DR The "better scores" question

The most comparable EDTA-only results across both repos are:

| Condition | Repo | Classes | Samples | F1 macro | Std |
|---|---|---|---|---|---|
| **EDTA FEM4 (256)** | `cup-classifier` | 5 | 96 | **0.7153** | 0.1365 |
| edta 3-family (660 feats) | `classifier` | 5 | 96 | 0.6599 | 0.0229 |
| edta 3-family + CpG all148 | `classifier` | 5 | 96 | 0.6661 | 0.0424 |
| edta 3-family + CpG tt39 | `classifier` | 5 | 96 | 0.6566 | 0.0430 |
| edta FEM4 alone | `classifier` | **4** | 56 | **0.7596** | 0.0368 |
| edta 3-family | `classifier` | **4** | 56 | **0.7395** | 0.0345 |
| edta 3-family + CpG tt39 | `classifier` | **4** | 56 | **0.7711** | 0.0306 |

**On the 5-class problem (colon/liver/prostate/stomach/healthyblood), `cup-classifier` scores are HIGHER.**

The "better scores" from the other repo are on the **4-class cancer-only** EDTA set (56 samples, no healthyblood). Excluding healthyblood — the largest and most separable class — makes the problem fundamentally easier. The score gap (~0.04–0.06) between the other repo's 4-class and this repo's 5-class is expected from class-count reduction alone, regardless of algorithm.

**If you were comparing 4-class (other) vs 5-class (current), class count explains most of the difference. If you were comparing the other repo's 5-class results (0.6599) against this repo's (0.7153), this repo is actually ahead.**

---

## Pipeline comparison by dimension

### 1. Classifier and hyperparameters

| | `cup-classifier` (Phase 1) | `classifier` (EDTA followup) |
|---|---|---|
| **Solver** | `saga` (supports L1) | `lbfgs` (L2 only) |
| **Penalty** | L1 (`l1_ratio=1`, elasticnet) | L2 |
| **C selection** | Inner 3-fold GridSearchCV, 6 values: `[0.001, 0.006, 0.04, 0.25, 1.58, 10]` | `LogisticRegressionCV(cv=3)` — internal CV over sklearn-default C grid |
| **Feature selection** | **Embedded** — L1 sparsity zeroes coefficients during optimization. No separate selection step. | **Explicit** — per-fold RFECV (`rfecv_infold_filter`, step=0.2, min_features=10) selects subset, THEN final L2 classifier trains on selected features |
| **Imputation** | Per-fold **mean** (custom `_impute_with_mean`), then `StandardScaler` | `SimpleImputer(strategy="median")` in Pipeline, then `StandardScaler` |
| **Max iterations** | 5000 | 5000 |

**Impact**: L1 with saga on 256 FEM4 features / 96 samples will zero coefficients during optimization. The C grid may not find the optimal sparsity level. The two-stage RFECV → L2 approach searches a fundamentally different space. **This is the single largest algorithmic difference.**

**Evidence**:
- `cup-classifier`: `scripts/run_phase1_pipeline.py` lines 224-254 (`LogisticRegression(l1_ratio=1, solver='saga')` with `GridSearchCV({'C': C_GRID})`)
- `classifier`: `scripts/denovo/bct_followup3_edta.py` lines 69-77 (`LogisticRegressionCV(penalty="l2", solver="lbfgs", cv=3, scoring="f1_macro", balanced)`) + RFECV at `scripts/denovo/phase5c_cohensd_selection.py`

### 2. Cross-validation

| | `cup-classifier` | `classifier` |
|---|---|---|
| **Folds × repeats** | 1 × 5-fold | 20 × 5-fold |
| **Stratification** | By `Tissue` | By `Tissue_Source` (e.g. `colon_Fox Chase`) |
| **Source-coverage check** | Yes — `source_covering_split()` rejects folds where any training fold lacks all sources. Retries up to 100× with incremented seed. | Not needed — `Tissue_Source` stratification inherently encodes source in the stratification key |
| **Seed** | `random_state=42` | `seed + rep` (42…61) |
| **Total test predictions** | 5 per modality | 100 per condition |
| **Metric aggregation** | Mean of 5 fold F1 values | Mean of 20 per-repeat F1 values (collect OOF predictions per repeat, compute macro-F1, then mean across repeats) |

**Impact**: The other repo's 20×5 CV gives dramatically more stable estimates (100 test predictions vs 5). The per-repeat F1 aggregation differs from simple fold-mean. Stratification by `Tissue_Source` vs `Tissue` changes which samples land in each fold, which matters for small-source classes (Audubon liver=9, Audubon stomach=10).

**Evidence**:
- `cup-classifier`: `scripts/run_phase1_pipeline.py` lines 66-97 (`source_covering_split`)
- `classifier`: `scripts/denovo/bct_phase_de_classification.py` lines 497-503 (`StratifiedKFold(n_splits=N_SPLITS, random_state=seed + rep)` with `strat_key = merged["Tissue"] + "_" + merged["Source"]`)

### 3. Data sources and sample composition

| | `cup-classifier` | `classifier` |
|---|---|---|
| **Raw metadata** | `input/metadata/metadata_cleaned.csv` — 164 samples after exclusions | `metadata_riskscores_all.csv` + Phase A BCT workbook — 186 confirmed BCT |
| **EDTA total** | **96** (healthy=40, stomach=21, prostate=17, colon=9, liver=9) | **111** raw (healthy=52, stomach=21, prostate=18, colon=9, liver=11) → **96** after 3-family inner join |
| **BCT types** | EDTA, Citrate, **ACD**, Streck | EDTA, Citrate, Streck (ACD excluded from `CONFIRMED_BCT`) |
| **Join key** | `TWIST_ID` | `Sample_ID` |
| **Healthy label** | `healthyblood` | `healthy` (normalized) |

**Key metadata difference**: `cup-classifier` has a more restrictive data-cleaning pipeline — 220 raw samples → 164 clean (ovary removed, gDNA excluded, metadata-only samples dropped, missing BCT excluded, duplicates collapsed). `classifier` keeps 186 "confirmed BCT" samples (less aggressive cleaning, relies on BCT filter + feature family join to drop samples).

The 15 extra EDTA samples in the other repo (12 healthy, 1 prostate, 2 liver) mostly drop during the 3-family inner join, landing both repos at 96 for 5-class. But the specific 96 survivors may differ due to different join keys (`TWIST_ID` vs `Sample_ID`) and different per-modality sample coverage.

**Evidence**:
- `cup-classifier`: `input/metadata/metadata_cleaned.csv` (164 rows); `CONTEXT.md` "Data quality" section
- `classifier`: `scripts/denovo/shared_utils.py` `filter_confirmed_bct()`; `scripts/denovo/bct_phase_de_classification.py` `assemble_noncpg_matrix()` inner join

### 4. Feature construction

| | `cup-classifier` | `classifier` |
|---|---|---|
| **FEM4 file** | Identical (`md5sum` matches) | Same file |
| **FEM4 loading** | Drops `['tissue', 'N_motifs']`, sets index on `sample` | Renames `sample`→`Sample_ID`, drops `EXCLUDE_TISSUES` rows (includes `healthy` by default), deduplicates on `Sample_ID` |
| **Fragment length** | `fragment_length_features_qc.csv` (QC-filtered — `tri_450_510` pre-removed) | `fragment_length_features.csv` (unfiltered; `tri_450_510` removed downstream by `clean_input_dataframe`) |
| **Feature set used in EDTA** | Individual modalities (FEM4 alone reported in summary) | 3-family: FEM4 + Fragment + CNVkit (±CpG variants in some conditions) |

**Evidence**:
- `cup-classifier`: `scripts/data_loading.py` `MODALITY_CONFIGS` — fragment: `*_qc.csv`, FEM4 drops `['tissue','N_motifs']`
- `classifier`: `scripts/denovo/shared_utils.py` `load_fem4_features()`, `load_fragment_features()`

### 5. EDTA filtering logic

| | `cup-classifier` | `classifier` |
|---|---|---|
| **Where** | `meta[meta['BCT'] == 'EDTA']` before loading features | `merged_all_h[merged_all_h['BCT'] == 'EDTA']` after feature family assembly |
| **When** | Metadata subset → load features for those samples | Load all features → inner-join families → filter to EDTA |
| **BCT source** | Direct column in `metadata_cleaned.csv` | Phase A `bct_metadata.csv` (derived from BCT XLSX) |
| **Pancreas** | 0 EDTA samples → automatically excluded | Same |

**Evidence**:
- `cup-classifier`: `scripts/run_phase1_pipeline.py` lines 965-970
- `classifier`: `scripts/denovo/bct_followup3_edta.py` lines 138-145

### 6. Metric calculation

| | `cup-classifier` | `classifier` |
|---|---|---|
| **Primary** | `f1_score(average='macro')` per fold → mean across 5 folds | Per-repeat F1: collect all OOF predictions in a repeat → `f1_score(average='macro')` → mean ± std across 20 repeats. Also reports "pooled F1" (all OOF together). |
| **Std calculation** | Across 5 folds (high variance with n=5) | Across 20 repeats (more stable) |
| **95% CI** | `mean ± 1.96 × std` (5-fold std) | Not reported |
| **OOF saved** | No — only per-fold metric table | Yes — full per-sample records in `oof_predictions/` |

**Evidence**:
- `cup-classifier`: `scripts/run_phase1_pipeline.py` lines 321-370 (per-fold metrics → mean)
- `classifier`: `scripts/denovo/bct_phase_de_classification.py` lines 550-570 (per-repeat F1 loop)

### 7. What is NOT different (EDTA-only)

These differences exist between the repos' wider pipelines but are NOT active in EDTA-only comparisons:

- **BCT correction** (`MeanShiftCorrector`): Not applied — EDTA is a single tube type. Confirmed: `bct_followup3_edta.py` conditions all have `corrected=False`.
- **CpG site condensation**: Only used in the `all148`/`tt39` conditions, not in `raw3` (3-family baseline).
- **Source correction** (`SourceCorrector`): Only in `benchmark_dimensionality_reduction.py`, not in EDTA followup.
- **PCA/PLS/Autoencoder**: Only in `benchmark_dimensionality_reduction.py`, not in EDTA followup.

---

## Root causes of F1 differences (ranked)

1. **Classifier: L1 (saga) vs L2 (lbfgs) + RFECV** — The core algorithmic divergence. L1 with saga on 256 features × 96 samples zeroes coefficients during optimization. The two-stage RFECV → L2 searches a fundamentally different space. The C grids also differ (6 values vs 3).

2. **CV stability: 1×5 vs 20×5** — 100 test estimates vs 5 changes both the mean and variance. The single 5-fold is a point estimate with high variance (std=0.1365 for the current repo's EDTA FEM4). The 20-repeat estimate has std≈0.03.

3. **Stratification axis: Tissue vs Tissue_Source** — Changes fold composition for small-source classes. `cup-classifier` has an explicit source-coverage retry loop; `classifier` relies on `Tissue_Source` key stratification.

4. **Feature set: 256 FEM4-only vs 660 3-family** — Adding fragment and CNVkit features (which have weaker signal in the EDTA subset, per feature ablation: Fragment F1=0.54, CNVkit F1=0.37) can dilute the stronger FEM4 signal.

5. **Data cleaning differences** — Though both end at 96 samples for 5-class, the specific survivors may differ due to different join keys (`TWIST_ID` vs `Sample_ID`) and different per-modality coverage.

6. **Fragment length source file** — Current repo uses QC-filtered version; other repo uses raw version. The `tri_450_510` column is dropped either way but other differences may exist.

---

## Unresolved questions

1. **Which specific EDTA results were you comparing?** The other repo's best 5-class EDTA (3-family) is 0.6599 — lower than this repo's 0.7153 (FEM4-only). The higher scores (0.7395–0.7711) are 4-class cancer-only. If you were comparing the 4-class results against this repo's 5-class results, class count is the dominant factor.

2. **Phase 1 summary truncated**: `phase1_summary.csv` only shows 2 rows (Full/EDTA FEM4) but all 10 modality×scope runs exist in the output directory. The summary was likely generated by a partial run (specific `--modalities` flag). I can regenerate it.

3. **Combined EDTA model**: Neither repo ran an EDTA-filtered combined model (FEM4 + probe_meth + fragment_length) that would be directly comparable. The pipeline supports it (`--combine` flag in Phase 1, `edta_4c_raw3` in the other repo) but the 5-class combined model isn't in the summary.
