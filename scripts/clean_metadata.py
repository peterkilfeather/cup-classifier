#!/usr/bin/env python3
"""
clean_metadata.py — Produce cleaned metadata CSV from raw input.

Usage:
    python scripts/clean_metadata.py

Reads:  input/metadata/Metadata_ALL_Sep2025_with_BCT_20260103.xlsx  (sheet 'Metadata')
Writes: input/metadata/metadata_cleaned.csv

Every exclusion and standardization is enumerated inline.
Decisions made via issue #7 (Data Quality and Cleaning),
confirmed by Peter in grilling conversation.
"""

import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "input/metadata/Metadata_ALL_Sep2025_with_BCT_20260103.xlsx"
OUT = REPO / "input/metadata/metadata_cleaned.csv"

# ── 1. Load raw metadata ──────────────────────────────────────────────────────
meta = pd.read_excel(RAW, sheet_name="Metadata", header=0)
# Keep only the 19 substantive columns; drop unnamed trailing cols
meta = meta.iloc[:, :19]

# Clean TWIST ID (strip whitespace)
meta["TWIST_ID"] = meta["TWIST ID"].str.strip()
meta = meta.drop(columns=["TWIST ID"])

# ── 2. Define drop sets (each tagged with reason) ─────────────────────────────
# Each entry: sample_id -> reason string
# Order matches the provenance log (docs/data-quality-provenance.md)

DROP = {}

# (a) gDNA samples — genomic DNA, not cfDNA
for sid in ("V10_S3", "V2_S1"):
    DROP[sid] = "gDNA sample (not cfDNA)"

# (b) "Not sent for sequencing" — confirmed zero feature data
for sid in ("V14_S8", "V37_S11", "V24_S5", "V24_S6"):
    DROP[sid] = "Sample not sent for sequencing"

# (c) V36_S7 duplicate — both rows, no feature data
DROP["V36_S7"] = "Duplicate metadata row; no feature data"

# (d) All ovary samples — no BCT annotation on any V39 sample
for i in range(2, 13):
    DROP[f"V39_S{i}"] = "Ovary sample with no BCT annotation"

# (e) Metadata-only samples — present in metadata, absent from all 4 feature modalities
META_ONLY = [
    "V12_S1", "V12_S2", "V13_S3",
    "V15_S9", "V16_S7", "V16_S10", "V16_S12",
    "V17_S1", "V17_S2", "V17_S3", "V17_S5", "V17_S8",
    "V28_S6",
    "V34_S10", "V35_S2", "V36_S13", "V37_S6",
]
for sid in META_ONLY:
    DROP[sid] = "Metadata-only sample: no feature data in any modality"

# (f) Missing BCT — have feature data but BCT is critical confound
BCT_MISSING = [
    "V14_S1", "V14_S3", "V14_S4",
    "V19_S3", "V19_S6", "V19_S12",
    "V24_S3", "V24_S4", "V28_S7",
    "V34_S1", "V34_S3", "V34_S4", "V34_S7", "V34_S9",
    "V36_S8", "V36_S14", "V36_S16",
    "V37_S2", "V37_S3", "V37_S4", "V37_S9",
]
for sid in BCT_MISSING:
    DROP[sid] = "Missing BCT annotation"

# ── 3. Apply drops ───────────────────────────────────────────────────────────
n_before = len(meta)
meta = meta[~meta["TWIST_ID"].isin(DROP)]
n_dropped = n_before - len(meta)

# Verify no dropped sample remains (belt-and-suspenders)
remaining_dropped = set(DROP) & set(meta["TWIST_ID"])
assert not remaining_dropped, f"Still present: {remaining_dropped}"

# ── 4. Standardizations ──────────────────────────────────────────────────────

# 4a. Source: Audobon → Audubon (same institution, correct spelling)
meta["Source"] = meta["Source"].str.replace("Audobon", "Audubon", regex=False)

# 4b. Race normalization
RACE_MAP = {
    "White": "White", "white": "White", "Caucasian": "White",
    "Black": "Black",
    "Asian": "Asian",
    "Other": "Other",
}
meta["Race_norm"] = meta["Race"].map(RACE_MAP).fillna("NA")

# 4c. Ethnicity normalization
ETH_MAP = {
    "Non-Spanish, Non-His": "Non-Hispanic",
    "Non-Spanish": "Non-Hispanic",
    "Spanish, NOS; Hispan": "Hispanic",
    "Puerto Rican": "Hispanic",
    "Ukraine": "Slavic",
    "Slavic": "Slavic",
    "Unknown Whether Span": "Unknown",
    "Unknown": "Unknown",
}
meta["Ethnicity_norm"] = meta["Ethnicity"].map(ETH_MAP).fillna("NA")

# ── 5. Select and order columns for downstream use ────────────────────────────
COLS = [
    "TWIST_ID", "Barcode", "Tissue", "Source", "BCT",
    "Sex", "Age", "Year Drawn",
    "Race_norm", "Ethnicity_norm",
    "Status", "Histology", "Tumor_site",
    "Pathologic Stage", "Clinical Stage", "Stage_type",
]
out = meta[COLS].rename(columns={"Year Drawn": "Year_Drawn"})
out = out.sort_values("TWIST_ID").reset_index(drop=True)

# ── 6. Write ──────────────────────────────────────────────────────────────────
out.to_csv(OUT, index=False)

# ── 7. Summary ────────────────────────────────────────────────────────────────
print(f"Wrote {len(out)} samples to {OUT}")
print(f"Dropped {n_dropped} samples ({len(DROP)} unique IDs in exclusion list)")
print()
print("Tissue distribution:")
for tissue, count in out["Tissue"].value_counts().sort_index().items():
    print(f"  {tissue}: {count}")
print()
print("BCT coverage: "
      f"{out['BCT'].notna().sum()} / {len(out)} "
      f"({out['BCT'].notna().mean() * 100:.0f}%)")
print()
print("Source values:", sorted(out["Source"].unique()))
