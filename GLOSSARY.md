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
