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

### Unknown origin (to be resolved)

**100kb End Bin Features**:
A set of features derived from fragment-end coverage at 100kb genomic bins. Present in `input/unreviewed/`. Origin and exact computation are not yet documented.

**Risk Score**:
A per-sample score in `metadata_riskscores_all.csv`. Origin is not yet documented.

## Clinical Site
One of the hospitals or tissue banks that contributed plasma samples: Fox Chase, Audubon, Sowalsky, NIH Clinical Center. Relevant for batch effect analysis.
_Avoid_: Source, hospital, cohort
