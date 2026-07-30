# Data Quality Provenance Log

**Issue**: #7 (Data Quality and Cleaning)
**Map parent**: #1 (Confound Mitigation Strategy)
**Date**: 2026-07-30
**Resolver**: wayfinder session

This document records every data quality issue identified, the decision made, and the rationale. All downstream tickets inherit this clean state.

---

## Summary

| Metric | Value |
|--------|-------|
| Raw metadata samples | 220 |
| Samples dropped | 56 |
| Clean metadata samples | 164 |
| Tissue classes remaining | 6 (ovary eliminated) |
| Feature-quality issues found | 1 (`tri_450_510` all-NaN) |

---

## Decisions

### 1. gDNA samples — DROP

**Samples**: V10_S3 (colon_gDNA), V2_S1 (prostate_gDNA)
**Reason**: These are genomic DNA from tumor tissue (gDNA_TumorCOAD/gDNA_TumorPRAD), not cfDNA. Processed by Laura's lab with an expired library prep kit. Not valid for cfDNA methylation/fragmentomic analysis.
**Reference**: Sheet3 in Metadata_ALL_Sep2025_with_BCT_20260103.xlsx
**Confirmed by Peter**: Yes

### 2. "Not sent for sequencing" — DROP

**Samples**: V14_S8 (prostate, Fox Chase), V37_S11 (prostate, Sowalsky), V24_S5 (liver, Audubon), V24_S6 (liver, Audubon)
**Reason**: Annotated "Sample not sent for sequencing" (V14_S8, V37_S11) or "Not enough sample to send for sequencing" (V24_S5, V24_S6). Confirmed zero feature data across all modalities.
**Confirmed by Peter**: Yes

### 3. V36_S7 duplicate — DROP

**Samples**: V36_S7 (healthyblood, NIH Clinical Center) — both metadata rows
**Reason**: Duplicate entry in metadata (two identical rows for NIH_21). No feature data in any modality.
**Confirmed by Peter**: Yes

### 4. Ovary samples (all) — DROP

**Samples**: V39_S2 through V39_S12 (11 ovary samples)
**Reason**: None of the 11 ovary samples have BCT annotation. Since BCT is a critical confound variable and cannot be reliably imputed, the entire ovary class is excluded.
**Note**: All 11 samples had complete feature data across all 4 modalities — this is a material loss of data.
**Confirmed by Peter**: Yes

### 5. Metadata-only samples (no feature data) — DROP

**Samples** (17):
- V12_S1, V12_S2, V13_S3, V17_S1, V17_S2, V17_S3, V17_S5 (healthyblood, Audubon, Ukraine 2023)
- V34_S10, V36_S13, V37_S6 (healthyblood, NIH Clinical Center, 2025)
- V35_S2 (healthyblood, Sereti lab — also missing Source and BCT)
- V15_S9, V16_S7, V16_S10, V16_S12, V17_S8, V28_S6 (pancreas, Fox Chase)

**Reason**: Present in metadata file but absent from all 4 feature modalities (methylation, FEM4, fragment length, CNVkit). No feature data available for any analysis.
**Confirmed by Peter**: Yes

### 6. Missing BCT with feature data — DROP

**Samples** (21):
| Count | Samples | Tissue | Source |
|-------|---------|--------|--------|
| 3 | V14_S1, V14_S3, V14_S4 | prostate | Fox Chase |
| 2 | V19_S3, V19_S6 | stomach | Audubon |
| 1 | V19_S12 | stomach | Fox Chase |
| 1 | V24_S3 | colon | Fox Chase |
| 1 | V24_S4 | liver | Audubon |
| 1 | V28_S7 | pancreas | Fox Chase |
| 4 | V37_S2, V37_S3, V37_S4, V37_S9 | prostate | Sowalsky |
| 8 | V34_S1, V34_S3, V34_S4, V34_S7, V34_S9, V36_S8, V36_S14, V36_S16 | healthyblood | NIH Clinical Center |

**Reason**: BCT is a critical confound variable directly confounded with tissue type. Samples with unknown BCT cannot be used in confound-mitigated analysis. Imputation was considered (e.g., NIH protocol uses EDTA, Sowalsky uses Streck) but rejected in favor of dropping.
**Confirmed by Peter**: Yes

### 7. Source naming: Audobon → Audubon — STANDARDIZE

**Scope**: 12 liver samples relabelled from "Audobon" to "Audubon"
**Reason**: Same institution (Ukraine hospital network). "Audubon" is the correct spelling and the majority usage (46 samples vs 12).
**Confirmed by Peter**: Yes

### 8. Ethnicity/race normalization — STANDARDIZE

**Race normalization mapping**:
| Raw | Normalized |
|-----|-----------|
| White, white, Caucasian | White |
| Black | Black |
| Asian | Asian |
| Other | Other |
| NaN/NA (Sowalsky, Origene, NIH) | NA |

**Ethnicity normalization mapping**:
| Raw | Normalized |
|-----|-----------|
| Non-Spanish, Non-His, Non-Spanish | Non-Hispanic |
| Ukraine, Slavic | Slavic |
| Spanish, NOS; Hispan, Puerto Rican | Hispanic |
| Unknown Whether Span, Unknown | Unknown |
| NaN/NA (Sowalsky, Origene, NIH) | NA |

**Confirmed by Peter**: Yes

### 9. V14_S1 whitespace — FIX

**Issue**: TWIST ID `V14_S1 ` (trailing space) in raw metadata did not match `V14_S1` in feature files.
**Fix**: Trim whitespace on all TWIST ID values. (Sample is dropped per decision #6 regardless.)
**Applied**: Yes

### 10. Fragment length feature `tri_450_510` — DROP

**Issue**: Feature is all-NaN (and zero-variance) across all samples.
**Fix**: Removed from fragment-length feature set. Cleaned file at `input/fragmentomic/fragment_length_features_qc.csv`.
**Status**: Done.

---

## Post-Cleaning Dataset Summary

**164 samples across 6 tissue classes:**

| Tissue | Count | Sources |
|--------|-------|---------|
| healthyblood | 40 | Audubon (27), NIH Clinical Center (13) |
| prostate | 38 | Fox Chase (18), Sowalsky (20) |
| liver | 24 | Fox Chase (12), Audubon (12) |
| colon | 22 | Fox Chase (22) |
| stomach | 21 | Fox Chase (11), Audubon (10) |
| pancreas | 19 | Fox Chase (19) |

**BCT coverage**: 100% — all remaining samples have BCT annotated.

**Confound notes (unchanged from CONTEXT.md)**:
- Source ~ Tissue: colon/ovary(removed)/pancreas = Fox Chase only; healthyblood = Audubon + NIH; prostate = Fox Chase + Sowalsky; stomach/liver = Fox Chase + Audubon
- Year ~ Tissue: all healthyblood from 2023 or 2025; all Audubon liver from 2024
- BCT ~ Tissue: ACD = colon only; Streck = prostate only (Sowalsky); Citrate = liver + pancreas; EDTA spans most tissues

---

## Feature-level QC results

All feature files were scanned for all-NaN and zero-variance features.

| Modality | File | All-NaN | Zero-variance | Action |
|----------|------|---------|---------------|--------|
| Methylation (probe_meth) | Enriched filtered wide | 0 | 0 | None needed |
| Methylation (probe_meth) | Enriched filtered long | 0 | 0 | None needed |
| Methylation (probe_meth) | Filtered (unenriched) long | 0 | 0 | None needed |
| Methylation (probe_meth) | Unfiltered long | 0 | 0 | None needed |
| Methylation (probe_cpg) | Unenriched filtered wide | 0 | ~2.6% sampled | Kept (biologically valid) |
| FEM4 | ALL_fem4_features.tsv | 0 | 0 | None needed |
| Fragment length | fragment_length_features.csv | **1** (`tri_450_510`) | 1 | Removed → `*_qc.csv` |
| CNVkit | thr 0.05 / 0.10 / 0.20 | 0 | 0 | None needed |
| End density | 100kb CPM matrix | 0 | ~7% sampled | Kept (biologically valid) |

**Dropped samples in feature files**: All files still contain rows for dropped samples (34–52 per file). These should be filtered by inner join with the cleaned metadata in downstream analysis — no separate cleaned copies produced per Peter's instruction.

**Zero-variance features**: Retained. These represent biologically constant signal (e.g., CpG sites with invariant methylation, genomic bins with no fragment-end coverage), not data errors. Confirmed by Peter.

## Files

- **Cleaning script**: `scripts/clean_metadata.py` — run to reproduce the cleaned CSV from raw input
- **Cleaned metadata**: `input/metadata/metadata_cleaned.csv` (164 samples, 14 columns)
- **Cleaned fragment-length features**: `input/fragmentomic/fragment_length_features_qc.csv` (tri_450_510 removed)
- **Raw metadata**: `input/metadata/Metadata_ALL_Sep2025_with_BCT_20260103.xlsx` (unchanged)
- **This log**: `docs/data-quality-provenance.md`
