# Complementarity of PERMANOVA and variance partitioning; Source as primary confound

Two lessons delivered back-to-back:

**Lesson 2 (Variance Partitioning):** User learned that VP fits a linear model per feature (`feature ~ Tissue + Source + BCT + Sex`) and decomposes variance via ANOVA Type II SS. The key insight: VP and PERMANOVA use completely different math (per-feature OLS vs distance-based permutation test) but converge on the same ranking: Tissue > Source > BCT. VP's BCT mean is slightly inflated by outlier features; PERMANOVA is the more conservative estimate. The correct interpretation is converging evidence across methods.

**Lesson 3 (Source as real confound):** User learned that Source differs from BCT in three critical ways: (1) higher unique R² (0.014–0.052 across modalities, ~3-10× BCT), (2) persists in EDTA subset (BCT held constant), proving it's independent of tube type, (3) classifier F1 for Source often exceeds Tissue F1. The practical implication: Source-stratified CV (train on one Source, test on another) is the correct mitigation. Batch correction would harm because Source is aliased with Tissue. Single-source tissues (colon, pancreas) are limited to within-batch comparisons.

**Evidence of understanding:** User requested these specific lessons after the first, indicating they understand the conceptual framework and want to drill into the next layers.

**Implications:** User now has the full confound picture. Next sessions could cover: (a) implementing Source-stratified CV in practice, (b) how to interpret the performance gap between standard and stratified CV, (c) designing a write-up of the confound analysis for collaborators.
