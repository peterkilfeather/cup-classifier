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
One of the values in the `Tissue` column of the cleaned metadata: `colon`, `liver`, `pancreas`, `prostate`, `stomach`, `healthyblood`. The prediction target for both classifiers. (Ovary was eliminated during data cleaning — no BCT annotation on any sample in that class.)
_Avoid_: Tissue type (ambiguous with Tissue-of-Origin)

**Full Classifier**:
The multiclass classifier that predicts across all 6 class labels: {healthyblood, colon, liver, pancreas, prostate, stomach}. Distinguishes both cancer-vs-healthy and tissue-of-origin among cancers.
_Avoid_: The full model, classifier v1

**Tissue-of-Origin (TOO) Classifier**:
The multiclass classifier that predicts among the 5 cancer tissue class labels: {colon, liver, pancreas, prostate, stomach}. Applied to samples already known to be cancer (e.g. clinically confirmed CUP).
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

**Probe Aggregation**:
Collapsing per-CpG beta values to probe level as an unweighted mean of observed site betas within a probe (sites with missing values skipped). Computed downstream from the per-CpG wide matrices. **Not the same as Probe-Averaged (probe_meth)**: probe_meth betas are read-weighted (`CpG_meth / CpG_total` over all reads in the 122 bp region), while aggregation gives every measured site equal weight and spans flanking-inclusive site sets (89/148 probes have site spans beyond 122 bp). Used as the same-data granularity baseline: per-CpG vs per-CpG-aggregated isolates the granularity question holding the measurement pipeline fixed.
_Avoid_: Probe-averaged (distinct term), probe-level mean

**Per-Sample Missingness**:
The fraction of a sample's features that are missing (NaN) in a feature matrix. Reported per version (e.g. per-CpG unenriched mean 6.0%, enriched 37.1%, post-join 164 samples). The control metric for every cross-version comparison.
_Avoid_: Missingness alone (ambiguous with per-feature missingness)

**Per-Feature Missingness**:
The fraction of samples where a feature is missing (NaN). Its extremes are the all-NaN features (NaN in every sample). Long-format methylation files express this structurally: an absent sample×probe (or sample×CpG) row is a missing value.
_Avoid_: Missingness alone (ambiguous with per-sample missingness)

**All-NaN Feature**:
A feature (e.g. a CpG site) with zero observed values across all samples. Dropped before any analysis; **un-imputable** — no observed values exist to impute from, and a mean-imputed all-NaN feature is a constant with zero variance, carrying no information. Counts: 70 (per-CpG unenriched), 1,156 (per-CpG enriched), 0 (all probe-averaged versions).
_Avoid_: "70 instances of NaNs" (they are features, not cells), "impute them"

**Coverage (Methylation)**:
The number of sequencing reads observed at a probe region or individual CpG site — the `n_reads` / `CpG_total` / `total_count` columns from which beta fractions are computed. Coverage drives missingness: sites below effective coverage are absent (NaN) from the feature matrix.
_Avoid_: Depth (colloquial; file columns use read counts)

**Fragmentomic Feature**:
A feature derived from cfDNA fragment properties. Includes:
- **Fragment length distribution**: short/long fragment ratios, modal length, bin counts across size windows
- **Fragment End-Density profile**: 5-prime 1bp fragment-end coverage at 100kb bins genome-wide, normalized to CPM
_Avoid_: Fragmentomics feature, fragment size feature

**Fragment End Motif (FEM4)**:
The frequency of each 4-mer sequence at the 5' end of cfDNA fragments. A 256-dimensional feature vector (4^4 = 256 possible tetramers) per sample. Stored in `input/fragmentomic/ALL_fem4_features.tsv`.
_Avoid_: End motif, fragment start motif, 4-mer

**Copy Number Feature (CNV)**:
A feature derived from cfDNA read-depth/coverage changes, reflecting chromosomal gains and losses in tumor-derived cfDNA. Processed via cnvkit.
_Avoid_: CNA (copy number alteration — correct in the paper, but CNV is broader)

**100kb Fragment End-Density Profile**:
A feature set derived from 5-prime 1bp fragment-end coverage at 100kb bins genome-wide, normalized to counts per million (CPM). The full genome-wide matrix (`input/fragmentomic/end_density/all_samples_ends_100kb_CPM_matrix.tsv`, ~31K bins) contains all bins with no label-informed selection — usable for feature selection inside cross-validation without leakage.

Two pre-selected subsets are archived at `input/archived/`:
- `all_samples_selected_100kb_bin_features.tsv`: 24 bins selected by comparing cancer tissues to healthy controls (FDR < 0.05, |delta mean CPM| >= 10)
- `all_samples_tissue_unique_100kb_bin_features.tsv`: 17 tissue-unique bins from that selection
- These pre-selected subsets are **flagged for potential data leakage** — the selection used the same samples, and any analysis using these files must account for the leaked label information.
_Avoid_: (none)

**tumortype39 / tt39 Panel**:
A set of 39 methylation probes (450K array IDs) provided for consideration in this project; no provenance record exists in this repo (selection basis unknown). Treated as leakage-adjacent pending evidence otherwise — any comparison including it must note this. Probe IDs are listed in `docs/references/tt39-probes.txt`. The original annotation file is archived at `input/archived/tumortype39_annotated_seq.txt`.

**Risk Score**:
A per-sample prediction probability from an earlier classifier model, stored in `input/archived/metadata_riskscores_all.csv`. Outdated — not relevant for current model development.
_Avoid_: (do not use for training or evaluation)

## Meeting Records

The meeting pipeline (current example `meetings/260731/`; pattern from `meetings/260611`–`260702` in the elnitski/classifier repo) produces a chain of artifacts with distinct roles. Transcription runs on the cherwell RTX 3090 from `/tmp` scratch while the local GPUs are serving-loaded; scratch is wiped after outputs are rsync'd back.

**Recording:**
The raw audio of a client meeting, `meetings/<YYMMDD>/recording_<YYMMDD>.m4a`. Renamed from the device's export name to the meeting-date convention; source of truth for everything downstream.
_Avoid_: keeping the device filename (e.g. "New Recording 69.m4a")

**Raw Transcript:**
The whisperx output: `out/<YYMMDD>/recording_<YYMMDD>.json` plus `.srt/.tsv/.txt/.vtt` renderings. Per-segment text, word timestamps, anonymous speaker labels (`SPEAKER_00`…). Never edited; the reference for all downstream artifacts.
_Avoid_: "transcript" bare (ambiguous with the corrected version)

**Speaker Turns:**
Consecutive same-speaker whisperx segments merged into turns (`<YYMMDD>.speaker_turns.txt`) by `post_process_<YYMMDD>.py`, which also prints QC (coverage %, segments/turns/words, per-speaker talk time). Still anonymous.

**Corrected Transcript:**
Speaker turns after `apply_corrections_<YYMMDD>.py` applies (a) a provisional speaker-name map and (b) domain-term ASR corrections (`<YYMMDD>.corrected.txt` + `corrections_report.md`). PROVISIONAL until the speaker map is confirmed — diarization gives anonymous labels; identities are inferred from content.
_Avoid_: treating the corrected transcript as authoritative before speaker confirmation

**Meeting Summary:**
The structured deliverable (`<YYMMDD>_summary.md`): Attendees, TL;DR, Decisions, Action Items, Analysis-Relevant Technical Points, Key Numbers, Open Questions, Pivotal Quotes (timestamped), Logistics. A synthesis written from the corrected transcript in the project's canonical language.
_Avoid_: copying transcript text verbatim into the summary (it is a synthesis)

**ASR Glossary (initial_prompt):**
The domain-term passage passed to whisperx `--initial_prompt` to bias transcription toward project vocabulary. Derived from this CONTEXT.md + GLOSSARY.md at the time of the meeting; reflects SPOKEN vocabulary (e.g. "source", "collection tube") even where written artifacts use canonical terms (Clinical Site, BCT). Kept byte-identical in `transcribe_<YYMMDD>.sh` to what actually ran.

## Clinical Site
One of the hospitals or tissue banks that contributed plasma samples: Fox Chase, Audubon (includes samples originally coded as "Audobon"), Sowalsky, NIH Clinical Center. Relevant for batch effect analysis.
_Avoid_: Source, hospital, cohort

## Batch Effects and Confounding

The dataset has potential confounding between tissue-of-origin and non-biological variables. These confounds must be addressed in experimental design and model validation.

### Source ~ Tissue confound

Most cancer tissue types come from a single clinical site:

| Tissue | Sources |
|--------|---------|
| colon | Fox Chase only (22) |
| pancreas | Fox Chase only (19) |
| prostate | Fox Chase (18) + Sowalsky (20) |
| liver | Fox Chase (15) + Audubon (9) |
| stomach | Fox Chase (11) + Audubon (10) |
| healthyblood | Audubon (27) + NIH Clinical Center (13) |

This means a classifier might learn site-specific technical artifacts (library prep, sequencer, storage) rather than biological tissue-of-origin signal. Samples from Fox Chase dominate most cancer classes, while healthy controls come exclusively from other sites.

### Year Drawn ~ Tissue confound

- All healthyblood from 2023 (27, Audubon) or 2025 (13, NIH)
- All 2024 liver (9, Audubon)
- Most prostate from 2018-2021

Temporal effects (reagent lot, sequencing run, protocol drift) could produce spurious correlations.

### BCT (Blood Collection Tube) ~ Tissue confound

- ACD tubes (a citrate formulation): only colon (13 samples)
- Citrate tubes: dominantly liver (15) + pancreas (19) + 1 prostate
- Streck tubes: only prostate (20, all Sowalsky)
- EDTA tubes: most diverse (healthyblood 40, stomach 21, prostate 17, colon 9, liver 9)

### Sex confound

Sowalsky samples (all prostate, 20 samples) are 100% male. Other sources skew male (Fox Chase 80%, Audubon 65%, NIH 54%). The dataset overall is male-dominated (125 of 164).

## Order Dependence

**Order-Independent Attribution (Type II / Type III / Marginal)**:
Variance attribution that does not depend on covariate order. Variance partitioning uses Type II sums of squares (each term's SS after all *other* terms, no interactions); PERMANOVA headline results use marginal (Type III) tests. Both answer "what does this covariate uniquely explain given everything else" — the answer does not change with covariate ordering.
_Avoid_: "in order" when describing these methods (slide wording pitfall — the covariate *list* has an order, the attribution does not)

**Entry-Order Sensitivity**:
The robustness check that *does* depend on order: sequential (Type I) runs with each covariate entering first/last, used to bracket shared variance between covariates. Numbers change across orderings (e.g. FEM4 BCT R² 0.155 entering first vs 0.018 marginal); conclusions are invariant (Tissue retains the largest unique component in every ordering tested).
_Avoid_: treating entry-order sensitivity as the headline result, or claiming "order doesn't matter" without the marginal-vs-sequential distinction

## Data Quality

Data cleaning was performed in issue #7 (see `docs/data-quality-provenance.md` for full decision record). The raw metadata (220 samples) was reduced to 164 clean samples.

### Exclusions

| Reason | Samples dropped |
|--------|----------------|
| gDNA samples (not cfDNA) | 2 |
| Not sent for sequencing | 4 |
| Ovary class (no BCT on any sample) | 11 |
| Metadata-only (no feature data) | 17 |
| Missing BCT (have feature data) | 21 |
| Duplicate metadata row | 1 |

### Cleaned feature modalities

| Modality | Samples present | Notes |
|----------|----------------|-------|
| Methylation (probe_meth) | 164 (after filtering) | Zero all-NaN probes across all variants |
| FEM4 | 164 (after filtering) | Zero all-NaN features |
| Fragment length | 196 (pre-join) | `tri_450_510` removed (all-NaN); cleaned file at `*_qc.csv` |
| CNVkit | 164 (after filtering) | Zero all-NaN features across all 3 thresholds |
| End density | 164 (after filtering) | Zero all-NaN bins |

Inner join with `metadata_cleaned.csv` is the expected filtering step.

### Standardizations applied

- **Source**: "Audobon" → "Audubon" (same institution, correct spelling)
- **Race**: White/Caucasian/white → White; values mapped to White/Black/Asian/Other/NA
- **Ethnicity**: Ukraine/Slavic → Slavic; Non-Spanish variants → Non-Hispanic; values mapped to Non-Hispanic/Slavic/Hispanic/Unknown/NA

### Per-feature missingness

Methylation data is in long format (one row per sample per probe), and not every sample-probe pair may be present if coverage was insufficient. This per-feature missingness must be handled downstream (imputation, min-coverage filtering).

## Data Leakage

Data leakage occurs when information from outside the training set influences model development, producing an optimistically biased estimate of generalization performance. In this project, leakage can arise from several sources:

### Feature selection using the full dataset

If features are selected (e.g., differentially methylated probes, informative fragment-end bins) by comparing training labels across the entire dataset before cross-validation, the feature set is contaminated with label information. This is the most common leakage pattern in genomic classifier development. The pre-selected 100kb bin files in `input/archived/` are flagged for this concern. The full genome-wide end-density matrix in the same directory has no pre-selection and is safe to use with feature selection inside cross-validation.

### Sample overlap between feature sets

Some samples in the feature files may have no matching entry in the methylation data (17 metadata-only samples remain in the cleaned set after other exclusions). If analyses are run on different sample subsets without tracking the overlap, leakage can occur when the same sample appears in both a discovery and validation set.

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

## Implications

- High classification accuracy does not confirm biology — the model may learn site or year. Performance on confound-free held-out samples (e.g., pancreas from Fox Chase only, no external validation) is unreliable.
- A truly rigorous evaluation requires at minimum cross-validation stratified by Source, or better, external validation from a source unseen during training.
- Feature importance analysis must distinguish between biologically meaningful features and features that correlate with batch (e.g., read depth, GC bias, fragment length distributions that differ by library prep).
