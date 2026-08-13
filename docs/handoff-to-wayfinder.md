# Handoff: cup-classifier domain model

## Repo

`github.com/peterkilfeather/cup-classifier` at `/xscratch/farney/cup-classifier`

## What we did

Domain modeling session. Built a comprehensive `CONTEXT.md` documenting the project's language, batch effects, data quality, and leakage risks.

## Project goal

Predict tissue-of-origin (colon, liver, pancreas, prostate, stomach, ovary, healthyblood) from cfDNA methylation and fragmentomic features in plasma from Cancer of Unknown Primary (CUP) patients.

Two multiclass classifiers:
- **Full Classifier** — 7-class (6 cancer tissues + healthyblood)
- **TOO Classifier** — 6-class (cancer tissues only)

Key aim: identify which feature sets and specific features drive classification performance.

## Data structure

```
input/
├── metadata/                       # metadata, BCT, source, tissue labels
├── methylation/
│   ├── probe_meth/                 # probe-averaged methylation (filtered, enriched, unfiltered)
│   └── probe_cpg/                  # per-CpG methylation (enriched + unenriched filtered)
├── fragmentomic/
│   ├── end_density/                # full 100kb genome-wide fragment-end CPM matrix (~31K bins)
│   ├── ALL_fem4_features.tsv       # 256 4-mer frequencies
│   ├── fragment_length_features.csv # 371 length features (ratios, bin counts)
│   └── fragle_stats_all.xlsx       # ctDNA burden scores (post-processing output)
├── cnvkit/                         # CNV features (3 thresholds) + raw .cns files
└── archived/                       # leakage-flagged and outdated files
```

## Key facts for wayfinding

### Confounds (all documented in CONTEXT.md)

The dataset has multiple interdependent confounds that make naive classification untrustworthy:

- **Source ~ Tissue**: colon/ovary/pancreas = Fox Chase only; healthyblood = Audubon + NIH; prostate = Fox Chase + Sowalsky; stomach = Fox Chase + Audubon; liver = Fox Chase + Audobon
- **Year ~ Tissue**: all healthyblood from 2023 or 2025; all liver from 2024
- **BCT ~ Tissue**: ACD (citrate) = colon only; Streck = prostate only (Sowalsky); Citrate = liver + pancreas; EDTA spans most tissues
- **Sex**: Sowalsky 100% male; dataset skews male
- **Audobon/Audubon**: same site, inconsistent spelling

### Missingness
220 metadata samples; 196-198 present in each feature modality. 22 metadata-only samples have no methylation/FEM4 data.

### Data leakage
Pre-selected 100kb bin files archived due to label-informed feature selection. Full genome-wide matrix is safe for feature selection inside CV.

### Related work
THEMIS paper (Bie et al. 2023, Nat Commun) uses the same EM-seq chemistry with similar modalities. Documented in `docs/related-work.md`.

### Feature counts
- Methylation probes: 148
- FEM4: 256 features
- Fragment length: 371 features
- CNVkit: 39 features
- End density: ~31K genome-wide bins

### tt39 panel
39 probes from previous group work, identified as informative for tissue discrimination. IDs in `docs/references/tt39-probes.txt`.

### Risk scores
Outdated predictions from an earlier model. Archived.

## Decisions captured at CONTEXT.md level

- Classifier architecture (2 multiclass classifiers)
- Feature format semantics (probe_meth vs probe_cpg, enriched vs unenriched, QC filtering)
- Probe annotations from 450K array
- THEMIS paper as relevant prior art
- All confounds documented but not yet mitigated

## Starting decisions for wayfinder (not exhaustive)

The decisions below were surfaced during domain modeling. The wayfinder may discover additional decisions and should not be limited to this list. Options listed are illustrative examples only — not a closed set. Research and brainstorming are expected.

1. **Confound mitigation strategy** — how to handle Source ~ Tissue, BCT, Year confounding
2. **Feature modality prioritization** — which features to start with and in what order
3. **Model selection** — which ML approach
4. **Evaluation protocol** — what constitutes a credible result given the confounds
5. **Feature importance methodology** — how to determine which features drive performance
6. **Missingness handling** — imputation strategy or sample exclusion
7. **Enriched vs unenriched methylation** — which capture method to use
