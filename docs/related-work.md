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
- 100-kb bins used for CNA and fragmentation — same bin size as the end-density profiles
- Demonstrated split-by-hospital cross-validation as a robustness check (directly relevant to the Source confound documented in CONTEXT.md)

**Key differences from this project:**
- THEMIS uses shallow whole-genome WMS; this project uses targeted probes (148 regions)
- THEMIS cohort is much larger (1,277 vs ~220) and focuses on screening (earlier stage)
- This project targets CUP (confirmed metastatic cancer); THEMIS targets general screening
- THEMIS developed its own feature extraction pipeline; this project uses CNVkit for CNV

---

## Methylation feature granularity: probe-averaged vs per-CpG (issue #8, step 2a)

**Sources:**
- Liu MC et al., *Ann Oncol* 2020;31:745–759 (GRAIL CCGA substudy 2). DOI: 10.1016/j.annonc.2020.02.011
- Klein EA et al., *Ann Oncol* 2021;32:1167–1177 (CCGA substudy 3, clinical validation). DOI: 10.1016/j.annonc.2021.05.806
- Conway AM et al. (CUPiD), *Nat Commun* 2024;15:3292. DOI: 10.1038/s41467-024-47195-7
- Stackpole ML et al. (cfMethyl-Seq), *Nat Commun* 2022;13:5566. DOI: 10.1038/s41467-022-32995-6
- Moran S et al. (EPICUP), *Lancet Oncol* 2016;17:1386–1395. DOI: 10.1016/S1470-2045(16)30297-2
- Shen SY et al., *Nature* 2018;563:579–583. DOI: 10.1038/s41586-018-0703-0

**Key findings:**
The clinically validated targeted cfDNA TOO classifiers all operate at region or fragment level, not per-CpG beta: GRAIL's panel targets 103,456 regions / 1.12M CpGs but models per-fragment methylation states within regions (CSO accuracy 88.7%, 87.0–90.2%, in true positives; TOO accuracy 93% in substudy 2); CUPiD aggregates 450K-array probes into 22,179 DMRs and reaches 96.8% TOO accuracy on cfDNA; cfMethyl-Seq uses ~117 bp region markers and reaches 89.1% TOO accuracy. Per-CpG-level features appear mainly in whole-methylome work (Shen 2018, cfMeDIP-seq, per-CpG WGBS-like resolution) and in tissue-array classifiers (EPICUP: 450K single-CpG probes, 97.7% sensitivity on tissue), where data are near-complete — cfDNA per-CpG data are not. No published head-to-head compares per-CpG vs probe-averaged features on identical targeted cfDNA data.

**Implication for this project:** Region-level (probe-level) aggregation is the validated norm in targeted cfDNA methylation; the per-CpG version screen directly measures whether per-CpG granularity adds signal on our 148 probes, a question the literature leaves open.

## Enriched (hybrid capture) vs unenriched cfDNA methylation (issue #8, step 2b)

**Sources:**
- Liu MC et al. 2020 (above): hybrid-capture (Twist) targeted methylation panel, median 139x unique on-target depth
- Buckley DN et al., *NAR Genom Bioinform* 2022;4:lqac099. DOI: 10.1093/nargab/lqac099
- Stackpole ML et al. 2022 (above): cfMethyl-Seq CpG-island enrichment vs WGBS
- Shen SY et al. 2018 (above): cfMeDIP-seq antibody enrichment on small cfDNA inputs

**Key findings:**
Hybrid capture concentrates sequencing on informative regions: GRAIL reaches 139x median on-target depth over 17.2 Mb/1.12M CpGs, a depth WGBS cannot afford genome-wide. Capture is reproducible (between-capture beta R² = 0.92, >90% target recovery, 172 plasma samples; Buckley) and agrees with WGBS gold-standard beta values (R² = 0.79). Unenriched alternatives exist: cfMethyl-Seq achieves >12x CpG-island enrichment over WGBS without capture, and cfMeDIP-seq uses immunoprecipitation on small inputs. The stated caveat for small targeted panels is fragment yield: they capture only a fraction of tumor cfDNA fragments, risking false negatives at low tumor fraction — GRAIL's panel is large partly for this reason. Cost tradeoff: capture adds bait/protocol cost but cuts sequencing cost per informative base; unenriched loses depth per target at equal read budget, which shows up as per-site missingness (observed here: 6% unenriched vs 37% enriched per-CpG — driven by the enriched site set including sparsely covered sites, not by sample failure).

**Implication for this project:** Enriched vs unenriched is a cost/depth tradeoff, not a missingness panacea; our two capture conditions differ mainly in which CpG sites are measurable and how completely, and this must be reported alongside any version-screen comparison.

## Missingness and imputation practice for per-CpG features (~6% here) (issue #8, step 2c)

**Sources:**
- Nardini C et al., *BMC Bioinformatics* 2020;21:268. DOI: 10.1186/s12859-020-03592-5
- Stackpole ML et al. 2022 (above): 10x coverage threshold for reliable beta
- Kang S et al. (CancerLocator), *Genome Biol* 2017;18:53. DOI: 10.1186/s13059-017-1191-5

**Key findings:**
Standard practice has two layers. (1) Minimum-coverage filtering: cfMethyl-Seq validates beta values against RRBS at ≥10x coverage (Pearson r = 0.987), and 10x is the commonly used per-CpG threshold. (2) Imputation: a 7-method benchmark (mean, kNN, etc.) found beta values impute better than M-values, mid-range betas are hardest to impute, and MAR/MNAR patterns are harder than MCAR — coverage-driven missingness in sequencing data is MNAR, so naive mean imputation is a ceiling, not a gold standard. Model-based alternatives handle low coverage explicitly: CancerLocator's per-CpG probabilistic model works on low-coverage plasma data without imputation. Our observed 6.0% (unenriched per-CpG) is modest and matches the issue's expectation; 37% (enriched) exceeds what mean imputation can faithfully repair.

**Implication for this project:** For per-CpG features, prefer coverage-threshold filtering and/or aggregation over imputation; at 37% missingness, per-fold mean imputation is questionable and missingness-matched comparisons (already mandated by the protocol) must report it.

## Per-CpG → probe-aggregation benchmarks (issue #8, step 2d)

**Sources:**
- Conway AM et al. 2024 (above): 450K probe values → 22,179 DMR features
- Liu MC et al. 2020; Klein EA et al. 2021 (above): region-level fragment-methylation features
- Stackpole ML et al. 2022 (above): ~117 bp region markers
- Moran S et al. 2016 (above): single-CpG (450K) tissue classifier

**Key findings:**
There is no published head-to-head benchmark of per-CpG vs probe-aggregated features for targeted cfDNA methylation classifiers — this is an explicit gap. The closest evidence is indirect and consistent: successful cfDNA systems aggregate (GRAIL regions, CUPiD DMRs, cfMethyl-Seq ~117 bp markers), while per-CpG features succeed on near-complete data (450K arrays; EPICUP) or in WGBS-scale work with aggressive feature selection (Shen 2018). Our own inventory provides the missing empirical anchor for this dataset: aggregating observed per-CpG values to probe level leaves 0 all-NaN probes per sample even at 37% missingness (enriched), and probe-level missingness is ≤0.7% in all 5 versions — so probe aggregation is a lossless (with respect to missingness) dimensionality reduction here.

**Implication for this project:** Probe-aggregation of per-CpG features is a defensible dimensionality-reduction baseline; whether per-CpG granularity adds signal beyond it must be settled by our own version screen, since the literature offers no direct benchmark.
