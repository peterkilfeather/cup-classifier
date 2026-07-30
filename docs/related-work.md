# Related Work

## THEMIS (THorough Epigenetic Marker Integration Solution)

Bie et al., *Nature Communications*, 2023. https://doi.org/10.1038/s41467-023-41774-w

**Platform:** Enzyme-based (TET2/APOBEC) whole-methylome sequencing of cfDNA at shallow depth (~2X, 60M reads). Same EM-seq chemistry used in this project.

**Cohort:** 497 healthy + 780 cancer patients (breast, colorectal, esophageal, gastric, liver, lung, pancreatic) across all stages.

**Feature modalities extracted from shallow WMS:**

| Modality | Description | Granularity |
|---|---|---|
| MFR | Fraction of fully methylated fragments per window | 1,846 windows of 1 Mb |
| FSI | Short/long fragment ratio (100-166 bp / 169-240 bp) | 502 windows of 5 Mb |
| CAFF | Copy number from size-selected short (<151 bp) and long (>220 bp) fragments | Chromosome arms |
| FEM | 256 × 4-mer fragment end motif frequencies | Per sample |

**Ensemble classifier:** Regularized logistic regression integrating all 4 modalities. AUC 0.966 for cancer detection. Modality contribution weights: FEM (0.58) > FSI (0.34) > MFR (0.33) > CAFF (0.06).

**Cancer Signal Origin localization:** Used TCGA ATAC-seq chromatin accessibility peaks (18 tissue-specific clusters) to profile methylation + fragment coverage. Random forest classifier achieved 54-65% accuracy across 7 cancer types.

**Relevance to this project:**
- Same chemistry (EM-seq) and same feature types (methylation, fragmentation, FEM, CNV)
- FEM was their strongest modality — supports the value of the FEM4 data in this project
- Their CSO approach is conceptually analogous to the 148 targeted probes used here
- 100-kb bins used for CNA and fragmentation — same bin size as the unreviewed end-density profiles
- Demonstrated split-by-hospital cross-validation as a robustness check (directly relevant to the Source confound documented in CONTEXT.md)

**Key differences from this project:**
- THEMIS uses shallow whole-genome WMS; this project uses targeted probes (148 regions)
- THEMIS cohort is much larger (1,277 vs ~220) and focuses on screening (earlier stage)
- This project targets CUP (confirmed metastatic cancer); THEMIS targets general screening
- THEMIS developed its own feature extraction pipeline; this project uses CNVkit for CNV
