# Confound Diagnostic Protocol

**Status:** Resolved via grilling session (2026-07-31)
**Map:** [#1 Confound Mitigation Strategy](https://github.com/peterkilfeather/cup-classifier/issues/1)
**Ticket:** [#2 Confound Diagnostic Protocol](https://github.com/peterkilfeather/cup-classifier/issues/2)
**Blocks:** #4 (Stratification and Batch Correction), #6 (Immediate Analysis While We Wait)

## Core Question

Which confounds (BCT, Source, Year, Sex) drive feature variation, and is BCT the dominant batch effect as hypothesized?

---

## 1. Feature Modalities

### Primary (well-conditioned, p << n or p ~ n)

| Modality | # features | Notes |
|----------|-----------|-------|
| Methylation (probe-averaged, enriched) | 148 | Biological signal, matches THEMIS MFR concept |
| FEM4 | 256 | THEMIS's strongest modality (weight 0.58) |
| Fragment length | 369 | Known tube-type sensitivity in literature |
| CNVkit | 38 × 3 thresholds | Low-dim, fast, THEMIS CAFF analog |

### Secondary (high-dim, p >> n, PCA-reduced)

| Modality | # raw features | Reduction |
|----------|---------------|-----------|
| Methylation (per-CpG, unenriched) | 32,084 | PCA → top K PCs (elbow-determined) |
| End density (100kb bins) | ~31,000 | PCA → top K PCs (elbow-determined) |

**Rationale:** All modalities run through the same pipeline for cross-modality comparison. High-dim modalities use PCA-first to enable distance-based and classifier methods without p >> n instability. Per-CpG methylation preserves within-probe CpG variation that probe-averaging dilutes — PCA naturally captures this.

### Data preparation

- **Inner join** all feature files with cleaned metadata (`input/metadata/metadata_cleaned.csv`, 164 samples)
- **Drop all-NaN features** before analysis (70 per-CpG sites confirmed all-NaN; zero in other modalities)
- **Zero-variance features** retained (biologically valid, see data-quality-provenance.md)
- **Per-CpG missingness**: ~6% per-sample on average; handle via feature-wise complete-case for variance partitioning, or impute (mean) for PCA — to be decided at implementation

---

## 2. Diagnostic Methods

All four methods run on every modality, producing converging evidence.

### 2.1 PCA + Confound Coloring

- **Input:** Feature matrix (samples × features), centered and scaled
- **Output:** PC1 vs PC2, PC2 vs PC3 scatter plots
- **Colorings (separate plots):**
  - Tissue (hue) — does biological signal appear?
  - BCT (hue) — does tube type separate samples?
  - Source (hue) — does clinical site drive clustering?
  - Sex (shape/hue) — is there a sex effect?
  - Year_Drawn (continuous color gradient) — temporal trends only, not a formal test
- **PPT value:** Side-by-side PC1×PC2 colored by Tissue vs colored by BCT. If BCT separation is cleaner, the slide makes the case in 2 seconds.

### 2.2 PERMANOVA (adonis) — Type III (Marginal)

- **Model:** `dist_matrix ~ Tissue + Source + BCT + Sex`
- **Distance:** Euclidean on centered/scaled features (PCA-reduced space for high-dim modalities; full feature space for low-dim)
- **Type:** Marginal (Type III) — each term's R² computed conditional on all others, enabling fair cross-term comparison
- **Entry-order sensitivity:** Also run each term first and last to bracket shared variance (required because BCT~Tissue Cramer's V = 0.756)
- **Output:** Stacked bar chart — R² per variable (+ residuals), one bar per modality, with error bars or ranges from entry-order sensitivity
- **PPT value:** Grid of stacked bars; if BCT's conditional R² consistently ≥ 0.10 across modalities, the case is quantitative

**Year_Drawn excluded** from the model: Cramer's V with Source = 0.904 (near-perfect collinearity), making it a proxy for Source/BCT, not an independent axis. Used as PCA color gradient instead.

### 2.3 Classifier Negative Controls

- **Targets:** Tissue (positive control), BCT (key test), Source (secondary test), shuffled labels (null)
- **Classifier:** Penalized logistic regression (L1 or elastic net) with nested cross-validation, OR random forest — simpler is better for interpretability
- **Evaluation:** Macro-averaged F1 (handles class imbalance)
- **Chance baseline:** Macro-F1 from uniform random classifier on the same class distribution
- **Output:** Grouped bar chart — macro-F1 per target per modality, with chance line
- **PPT value:** If BCT-prediction F1 rivals Tissue-prediction F1, that's the money slide

### 2.4 Per-Feature Variance Partitioning

- **Method:** For each feature (probe, CpG, motif, bin, fragment-length bin, CNV segment), fit linear model: `feature ~ Tissue + Source + BCT + Sex`. Decompose variance via ANOVA sum-of-squares.
- **Computational cost:** Negligible — 32K linear models on 164 samples runs in seconds
- **Applied to all modalities** — not PCA-reduced. The per-feature decomposition is consistent and interpretable across low-dim and high-dim modalities.
- **Output:** Violin/box plot — distribution of variance fraction explained by each covariate across all features. Heatmap of top features most driven by each confound.
- **PPT value:** "The typical feature has X% of its variance explained by BCT vs Y% by Tissue" — a single compelling number.

---

## 3. Thresholds

A confound is **confirmed dominant** only if both PERMANOVA AND classifier criteria are met.

### PERMANOVA (Type III conditional R²)

| Threshold | Interpretation |
|-----------|---------------|
| R²(confound \| others) < 0.05 | Weak — unlikely to distort downstream analysis |
| 0.05 ≤ R² < 0.10 | Detectable — note in supplement, may need mitigation |
| 0.10 ≤ R² < 0.25 | **Material** — should be addressed in stratification/batch correction |
| R² ≥ 0.25 | **Dominant** — confound explains more unique variance than any biological term |

R² ≥ 0.10 is the primary action threshold (domain-standard for batch-effect materiality in omics).

### Classifier Negative Control

| Threshold | Interpretation |
|-----------|---------------|
| F1 < 2× chance | Confound not learnable from these features |
| 2× chance ≤ F1 < 0.8 × Tissue-F1 | Confound is learnable but weaker than biological signal |
| F1 ≥ 0.8 × Tissue-F1 | **Confound signal rivals biological signal** — strong evidence of dominance |

### Final adjudication

| PERMANOVA | Classifier | Verdict |
|-----------|-----------|---------|
| ≥ 0.10 | F1 ≥ 0.8×Tissue-F1 | **Dominant confound** — must mitigate |
| ≥ 0.10 | < 0.8×Tissue-F1 | Material confound — design around it |
| < 0.10 | F1 ≥ 0.8×Tissue-F1 | Borderline — investigate further (feature selection artifact?) |
| < 0.10 | < 2× chance | Weak confound — note and proceed |

---

## 4. Scope

### Primary: Full dataset (164 samples)

All 6 tissues, all BCT types. Establishes overall confound landscape.

### Sensitivity: EDTA-only subset (96 samples)

BCT held constant. Asks: "does biological tissue signal persist when tube type is removed?" Tissues represented: healthyblood (40), stomach (21), prostate (17), colon (9), liver (9). Pancreas excluded (all Citrate).

**Interpretation of sensitivity comparison:**

| Full dataset | EDTA subset | Inference |
|-------------|-------------|-----------|
| BCT dominates | Tissue signal strong | Real biology, BCT was noise |
| BCT dominates | Tissue signal collapses | BCT was carrying apparent tissue signal |
| Tissue dominates | Tissue signal strong | BCT is a minor confound |
| Tissue dominates | Tissue signal collapses | Something weird — investigate |

---

## 5. Implementation Notes

### Software
- **PCA**: scikit-learn `PCA` with scaling
- **PERMANOVA**: `skbio.stats.distance.permanova` or `vegan::adonis2` via R
- **Classifier**: scikit-learn `LogisticRegression(penalty='l1', solver='saga')` with `GridSearchCV`
- **Variance partitioning**: custom ANOVA decomposition or `statsmodels` OLS
- **Plots**: `matplotlib` + `seaborn` for publication-quality figures

### Modality-specific preprocessing
- **High-dim (per-CpG, end density):** PCA fit on full feature matrix (centered, scaled) → retain top K PCs. **K selection rule:** elbow of scree plot (derivative-based), capped at max(20, floor(n_samples/5)). Elbow prioritised over cap; if elbow > cap, use cap and report cumulative variance at chosen K. K will typically be 15–25 for 164 samples. Use these PCs as the "modality representation" for PERMANOVA and classifier.
- **Per-feature variance partitioning** uses raw features directly (not PCA-reduced), consistently across all modalities. High-dim modality PCs are only used for PERMANOVA and classifier methods.
- **Fragment length:** 369 features, already well-conditioned. No PCA needed.
- **Low-dim (probe_meth, CNV):** Use raw features directly.
- **CNVkit thresholds:** Run all three (0.05, 0.10, 0.20) separately; report the one with strongest confound signal (conservative — threshold = parameter choice, not p-hacking).

### Per-CpG missingness
~6% per-sample missingness. Three options, to be decided at implementation:
1. Mean-impute per CpG (fast, conservative — shrinks variance)
2. Feature-wise complete-case (drops some CpGs if missing in any sample)
3. kNN imputation (more accurate but slower)
Recommendation: start with (1) for PCA, (2) for variance partitioning.

---

## 6. Deliverables

Per modality, per scope:

- **PCA grid** (4 panels: Tissue, BCT, Source, Sex) + scree plot
- **PERMANOVA stacked bar** with entry-order sensitivity bounds
- **Classifier bar chart** (Tissue vs BCT vs Source F1, chance line, per modality)
- **Variance partition violin plot** per covariate
- **Adjudication table** — verdict per (modality, confound) pair
- **Summary slide** — "BCT is/is not the dominant batch effect" with supporting panel

All figures designed for PowerPoint insertion (no academic figure cruft; clean labels, readable legends).
