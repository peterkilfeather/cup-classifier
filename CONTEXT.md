# Cup Classifier

A classifier that predicts cancer tissue-of-origin from cell-free DNA (cfDNA) methylation and fragmentomic features in plasma from patients with Cancer of Unknown Primary (CUP). A key aim is to identify which feature sets — and which specific features within them — drive classification performance.

## Language

**Cancer of Unknown Primary (CUP)**:
A metastatic cancer whose original tissue-of-origin cannot be determined by standard clinical workup. The classifier's clinical use case.
_Avoid_: Cancer of unknown primary site, occult primary

**Tissue-of-Origin**:
The anatomical site or tissue type where the cancer originated (e.g. colon, lung, liver, pancreas, stomach, prostate, ovary). One of the classifier's prediction targets.
_Avoid_: Tumor type, cancer type, primary site (though used in data files — these are synonyms)

**Class Label**:
One of the values in the `Tissue` column of the metadata: `colon`, `liver`, `pancreas`, `prostate`, `stomach`, `ovary`, `healthyblood`. The prediction target for both classifiers.
_Avoid_: Tissue type (ambiguous with Tissue-of-Origin)

**Full Classifier**:
The multiclass classifier that predicts across all 7 class labels: {healthyblood, colon, liver, pancreas, prostate, stomach, ovary}. Distinguishes both cancer-vs-healthy and tissue-of-origin among cancers.
_Avoid_: The full model, classifier v1

**Tissue-of-Origin (TOO) Classifier**:
The multiclass classifier that predicts among the 6 cancer tissue class labels: {colon, liver, pancreas, prostate, stomach, ovary}. Applied to samples already known to be cancer (e.g. clinically confirmed CUP).
_Avoid_: TOO-only model, reduced classifier

**Healthy Control**:
A plasma sample from an individual with no cancer diagnosis. The class label in the data is `healthyblood`. Functions both as a negative class (Full Classifier) and as a baseline normalization reference.
_Avoid_: Normal, non-cancer

**Cell-Free DNA (cfDNA)**:
DNA fragments shed into the bloodstream by dying cells, including both normal and tumor cells. The classifier's input material.
_Avoid_: ctDNA (circulating tumor DNA is a subset, not a synonym)

**Methylation Feature Format**:
Whether the methylation data is averaged per probe or kept per individual CpG site. Two variants:
- **Probe-averaged** (`probe_meth`): methylation values are collapsed across all CpG sites within each 122bp probe region. The feature is the probe-level beta fraction.
- **Per-CpG** (`probe_cpg`): each individual CpG site within a probe is a separate feature.
_Avoid_: Long format, wide format (these describe file shape, not feature semantics)

**Capture Method**:
Whether a targeted enrichment was performed during library preparation.
- **Enriched** (`_enriched`): the cfDNA library is hybridized to biotinylated RNA/DNA baits matching the 148 probe regions, then physically pulled down with streptavidin beads before sequencing. Higher coverage on target regions.
- **Unenriched** (`_unenriched`, no suffix): no pull-down step. Lower expected coverage on the 148 probe regions.
_Avoid_: Targeted, untargeted (ambiguous with other enrichment methods)

**Bisulfite Conversion QC (mNonCpGlt4)**:
A QC filter that excludes samples where the non-CpG methylation fraction (`nonCpG_frac`) is >= 0.04. In somatic tissues, CpH methylation should be essentially zero; elevated values indicate incomplete bisulfite conversion, which would produce unreliable CpG methylation measurements. Files tagged `mNonCpGlt4` have passed this filter; files tagged `unfiltered` have not.
_Avoid_: (none — this is the standard term)

**Probe (Methylation)**:
One of 148 targeted genomic regions (122 bp each) designed to capture specific CpG sites relevant to tissue-of-origin discrimination. Probes were chosen from the Illumina 450K array and are assay-independent (usable across WGBS, EM-seq, and TWIST panels). Annotations are in `probe_annotations_450k.csv`.
_Avoid_: Panel, amplicon, region (too generic)

**Methylation Beta Fraction**:
The fraction of methylated molecules at a given probe or CpG site, ranging 0–1. Calculated as `CpG_meth / CpG_total`. The primary methylation feature.
_Avoid_: Beta score, methylation level, CpG_frac (though used in data files — this is a synonym)

**Fragmentomic Feature**:
A feature derived from cfDNA fragment length distribution, including short/long fragment ratios (e.g. short_long_ratio_165), modal length, and bin counts across size windows.
_Avoid_: Fragmentomics feature, fragment size feature

**Fragment End Motif (FEM4)**:
The frequency of each 4-mer sequence at the 5' end of cfDNA fragments. A 256-dimensional feature vector (4^4 = 256 possible tetramers) per sample. Stored in `input/fragmentomic/ALL_fem4_features.tsv`.
_Avoid_: End motif, fragment start motif, 4-mer

**Copy Number Feature (CNV)**:
A feature derived from cfDNA read-depth/coverage changes, reflecting chromosomal gains and losses in tumor-derived cfDNA. Processed via cnvkit.
_Avoid_: CNA (copy number alteration — correct in the paper, but CNV is broader)

**100kb Fragment End-Density Profile** (on hold):
A feature set derived from 5-prime 1bp fragment-end coverage at 100kb bins genome-wide, normalized to counts per million (CPM). Bins were selected by comparing each cancer tissue to healthy controls using a significance + effect-size filter (FDR < 0.05, |delta mean CPM| >= 10), yielding 24 candidate bins (17 tissue-unique, 7 shared). Present in `input/unreviewed/` as `all_samples_selected_100kb_bin_features.tsv` and `all_samples_tissue_unique_100kb_bin_features.tsv`.
_Avoid_: (on hold pending clarification — the bin selection may have used the same samples intended for classifier training, which would cause data leakage)

**Risk Score**:
A per-sample score in `metadata_riskscores_all.csv`. Origin is not yet documented.

## Clinical Site
One of the hospitals or tissue banks that contributed plasma samples: Fox Chase, Audubon, Sowalsky, NIH Clinical Center. Relevant for batch effect analysis.
_Avoid_: Source, hospital, cohort

## Batch Effects and Confounding

The dataset has potential confounding between tissue-of-origin and non-biological variables. These confounds must be addressed in experimental design and model validation.

### Source ~ Tissue confound

Most cancer tissue types come from a single clinical site:

| Tissue | Sources |
|--------|---------|
| colon | Fox Chase only |
| ovary | Fox Chase only |
| pancreas | Fox Chase only |
| prostate | Fox Chase + Sowalsky (25 of 47 from Sowalsky) |
| liver | Fox Chase (15) + Audobon (12) |
| stomach | Fox Chase (12) + Audubon (12) |
| healthyblood | Audubon (34) + NIH Clinical Center (26) |

This means a classifier might learn site-specific technical artifacts (library prep, sequencer, storage) rather than biological tissue-of-origin signal. Samples from Fox Chase dominate most cancer classes, while healthy controls come exclusively from other sites.

### Year Drawn ~ Tissue confound

- All 2023 healthyblood (34, Audubon)
- All 2025 healthyblood (26, NIH)
- All 2024 liver (12, Audobon)
- Most prostate from 2018-2021

Temporal effects (reagent lot, sequencing run, protocol drift) could produce spurious correlations.

### BCT (Blood Collection Tube) ~ Tissue confound

- ACD tubes: only colon (13 samples)
- Citrate tubes: dominantly liver (15) + pancreas (25)
- Streck tubes: only prostate (21, all Sowalsky)
- EDTA tubes: most diverse (healthyblood, stomach, prostate, liver, colon)

### Sex confound

Sowalsky samples (all prostate) are 100% male. Other sources skew male (Fox Chase 69%, Audubon 63%). The dataset overall is male-dominated.

### Ethnicity coding inconsistency

Same concept recorded differently across sites: "White" vs "Caucasian" vs "white"; "Ukraine" vs "Slavic"; Sowalsky/Origene use "NA". This makes ethnicity analysis unreliable without normalization.

### Audobon / Audubon

Liver samples are coded as Source="Audobon", stomach and healthyblood as "Audubon" — likely the same institution with inconsistent spelling in the metadata.

### Implications

- High classification accuracy does not confirm biology — the model may learn site or year. Performance on confound-free held-out samples (e.g., pancreas from Fox Chase only, no external validation) is unreliable.
- A truly rigorous evaluation requires at minimum cross-validation stratified by Source, or better, external validation from a source unseen during training.
- Feature importance analysis must distinguish between biologically meaningful features and features that correlate with batch (e.g., read depth, GC bias, fragment length distributions that differ by library prep).

## Data Leakage

Data leakage occurs when information from outside the training set influences model development, producing an optimistically biased estimate of generalization performance. In this project, leakage can arise from several sources:

### Feature selection using the full dataset

If features are selected (e.g., differentially methylated probes, informative fragment-end bins) by comparing training labels across the entire dataset before cross-validation, the feature set is contaminated with label information. This is the most common leakage pattern in genomic classifier development. The 100kb fragment end-density profiles (see above) are flagged for this concern pending clarification.

### Sample overlap between feature sets

Some samples in the feature files may have no matching entry in the methylation data (22 metadata samples are missing from the probe-meth and FEM4 feature tables). If analyses are run on different sample subsets without tracking the overlap, leakage can occur when the same sample appears in both a discovery and validation set.

### Preprocessing informed by all samples

Normalization steps applied across all samples together (rather than within each cross-validation fold) leak information between train and test splits. Examples include:
- Global mean centering or scaling
- PCA or dimensionality reduction fitted on all samples
- Imputation of missing values using statistics from the full dataset

### Avoiding leakage

- Any feature selection or dimensionality reduction must be performed inside the cross-validation loop, on training folds only.
- Normalization parameters (means, scales, PCA loadings) must be estimated on training folds and applied to held-out folds.
- Sample overlap between feature sets must be tracked and documented explicitly.
- External validation on an independent cohort (different source, different collection period) is the strongest guard against both leakage and batch confounding.
