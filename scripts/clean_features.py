#!/usr/bin/env python3
"""
clean_features.py — Remove known-bad features from feature matrices.

Usage:
    python scripts/clean_features.py

Reads raw feature files, drops features flagged by QC (all-NaN, etc.),
writes cleaned copies. Currently handles:

  - fragment_length_features.csv → removes tri_450_510 (all-NaN)

Extend the BAD_FEATURES dict as new issues are discovered.
"""

import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# ── Registry of bad features ──────────────────────────────────────────────────
# Each entry: relative_input_path -> { 'output_suffix': str, 'drop_cols': list }
# Reasoning documented in docs/data-quality-provenance.md

BAD_FEATURES = {
    "input/fragmentomic/fragment_length_features.csv": {
        "output_suffix": "_qc",
        "drop_cols": ["tri_450_510"],
        "reason": "All-NaN (and zero-variance) across all samples",
    },
}


def main():
    for rel_path, spec in BAD_FEATURES.items():
        src = REPO / rel_path
        if not src.exists():
            print(f"[SKIP] {rel_path} — not found")
            continue

        df = pd.read_csv(src)
        print(f"\n[READ] {rel_path}  ({df.shape[0]} rows × {df.shape[1]} cols)")

        before = df.shape[1]
        missing = [c for c in spec["drop_cols"] if c not in df.columns]
        if missing:
            print(f"  [WARN] Columns not found in file: {missing}")
        present = [c for c in spec["drop_cols"] if c in df.columns]
        if not present:
            print(f"  [SKIP] Nothing to drop")
            continue

        df = df.drop(columns=present)
        after = df.shape[1]
        print(f"  [DROP] {present} ({spec['reason']})")

        # Build output path
        stem = src.stem  # e.g. "fragment_length_features"
        suffix = src.suffix  # e.g. ".csv"
        out_name = f"{stem}{spec['output_suffix']}{suffix}"
        out_path = src.parent / out_name
        df.to_csv(out_path, index=False)
        print(f"  [WRITE] {out_path.relative_to(REPO)}  ({df.shape[0]} rows × {after} cols)")

    print("\nDone.")


if __name__ == "__main__":
    main()
