# Confound Analysis Resources

## Knowledge

- [Article: "What is PERMANOVA?" — Michael J. Anderson (2001)](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x)
  The original PERMANOVA paper. Use for: understanding the permutation-based distance test.
- [Vignette: `vegan::adonis2` documentation — Jari Oksanen](https://rdrr.io/github/vegandevs/vegan/man/adonis2.html)
  R function reference for the actual implementation used here. Use for: the `by=` argument (marginal vs terms) and how permutations work.
- [Article: "What is a marginal (Type III) test?" — UCLA IDRE Stats](https://stats.oarc.ucla.edu/other/mult-pkg/whatstat/)
  Use for: the difference between Type I, II, and III sums of squares.
- [Book: *Model Selection and Multimodel Inference* — Burnham & Anderson (2002)](https://doi.org/10.1007/b97636)
  Use for: the conceptual framework of variance partitioning and evidence ratios.
- [Article: "Confounding and batch effects in omics" — Leek et al. (2010)](https://doi.org/10.1038/nrg2825)
  Use for: the difference between a confound that matters (unique variance) and one that's merely correlated.

## Wisdom (Communities)

- [r/bioinformatics](https://reddit.com/r/bioinformatics)
  Use for: practical omics confound questions, batch correction strategies.
- [Bioconductor Support Forum](https://support.bioconductor.org/)
  Use for: PERMANOVA, variance partitioning, and batch-correction tool questions (tag `vegan`, `adonis`).
