# Mission: Understand how confound analysis works — and why BCT is not a dominant confound

## Why
I've spent 3 months assuming BCT (blood collection tube type) is the dominant batch effect in my cfDNA tissue-classifier dataset. The diagnostic protocol just told me BCT's unique contribution is negligible — and I need to understand *why* the statistics say this, so I can trust the result and explain it to collaborators. Deeper goal: build reliable intuition for interpreting PERMANOVA, variance partitioning, and confound diagnostics in future omics work.

## Success looks like
- I can explain the difference between marginal (Type III) and sequential (Type I) PERMANOVA to a colleague
- I can articulate why a classifier can "predict" BCT even when BCT has near-zero unique variance
- I know when to trust a PERMANOVA result and when to be skeptical
- I can design my own confound diagnostics for future datasets

## Constraints
- Learning happens in the cup-classifier repo workspace — use our actual data
- Lessons should be short, focused, with interactive elements where possible

## Out of scope
- Mathematical derivations of PERMANOVA or ANOVA formulas
- General R/Python programming (we use these as tools, not subjects)
