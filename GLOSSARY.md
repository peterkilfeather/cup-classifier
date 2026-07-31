# Confound Analysis Glossary

Terms used in the cup-classifier confound diagnostic protocol.

## PERMANOVA
**Permutational Multivariate Analysis of Variance**: a distance-based method that partitions variance among covariates and tests significance by permuting rows of the distance matrix. No normality assumption. Implemented in R as `vegan::adonis2`.

## Marginal (Type III) test
Each term's contribution evaluated *last* — after all other terms are already in the model. Answers: "what does this variable explain that nothing else already explains?" This is the strictest test of a term's unique contribution.

## Sequential (Type I) test
Terms evaluated in the order they enter the model. Earlier terms get credit for variance shared with later terms. Answers: "how much does this variable explain if it gets first pick?" Sensitive to entry order.

## Shared variance
Variance in the features that could be explained by either of two correlated covariates. Example: ACD tubes are only used for colon samples — variance from "ACD vs EDTA" is indistinguishable from "colon vs other tissues."

## Unique variance
Variance that can only be explained by one covariate, not any other. The quantity measured by the marginal (Type III) test.

## R²
Fraction of total variance in the distance matrix (PERMANOVA) or a single feature (variance partitioning) that a covariate explains. Ranges 0–1.

## F1 (Macro)
Harmonic mean of precision and recall, averaged across classes. Used for the classifier negative control because it handles class imbalance better than accuracy.

## Variance partitioning
Fitting a linear model per feature (`feature ~ covariates`), then decomposing each feature's variance via ANOVA sums of squares into fractions attributable to each covariate.

## Elbow K
The number of principal components chosen by finding the largest drop in singular values (second derivative maximum). Used to reduce high-dim modalities before PERMANOVA and classifier.

## Batch correction
Methods (ComBat, limma::removeBatchEffect, harmony) that adjust feature values to remove batch-associated variation. **Not needed** when the batch effect's unique variance is below threshold and the batch is aliased with biology.

## Entry-order sensitivity
Running PERMANOVA with covariates in different orders to bracket shared variance. The range of R² values (marginal, sequential, reversed-order) shows how much variance is shared between covariates.

## CRM (Cramer's V)
A measure of association between categorical variables, ranging 0 (independent) to 1 (perfect correlation). Chi-square-based, analogous to Pearson correlation for categories.

---

# Phase 1 Pipeline Glossary

## Source-stratified CV
A cross-validation strategy that splits data by tissue class (StratifiedKFold) and then verifies that every training fold contains all sources present in the dataset. If a fold is missing a source, the split is rejected and re-shuffled (up to 100 retries). Prevents the classifier from learning source-specific artifacts rather than biological tissue signal.

## L1-penalized logistic regression (multinomial)
A multiclass logistic regression with L1 (lasso) regularization, which drives irrelevant feature coefficients to zero — performing implicit feature selection. Used with `solver='saga'` (the only solver supporting L1 + multinomial in scikit-learn). The `C` parameter controls regularization strength (lower = stronger). Tuned via inner GridSearchCV per fold.

## Per-fold PCA
Principal Component Analysis fit on training data only, then applied to transform both training and test sets. Prevents data leakage that would occur if PCA were fit on the full dataset before cross-validation. Used for high-dimensional modalities (per-CpG methylation with ~32K features, end density with ~31K bins). K=20 fixed components (not elbow-selected).

## Imputation inside CV
Replacing missing values (NaNs) using column means computed from the training set only, then applying those same means to the test set. Contrast with pre-CV imputation, which would leak test-set statistics into training. Used for low-dimensional modalities (FEM4, probe-averaged methylation, fragment length) that may have sporadic missing values.
