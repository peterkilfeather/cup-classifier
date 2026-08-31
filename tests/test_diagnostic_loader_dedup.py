"""Regression test: duplicate sample-row dedup in run_diagnostic_parallel.load_modality.

Ties to issue #11 (FEM4 n=169 bug) and #17. The FEM4 feature file
(input/fragmentomic/ALL_fem4_features.tsv) contained duplicate sample rows
(V17_S16, V7_S10, V7_S13, V7_S4, V7_S9 each twice — 198 unique IDs over 203
rows). The diagnostic loader lacked the dedup step that data_loading.py has,
so 164 unique samples + 5 duplicated rows = 169 passed the metadata mask,
inflating every FEM4 PERMANOVA/VP/classifier run. The fix lives in
run_diagnostic_parallel.py::load_modality (df.index.duplicated(keep='first')).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

import data_loading
import run_diagnostic_parallel as rdp

# The 5 sample IDs duplicated in the real ALL_fem4_features.tsv (issue #11).
DUP_SAMPLES = ['V17_S16', 'V7_S10', 'V7_S13', 'V7_S4', 'V7_S9']


def _write_fixture(path, rows):
    """Write a feature TSV (sample col + 2 feature cols) with duplicated sample rows."""
    lines = ['sample\tF1\tF2']
    for s, f1, f2 in rows:
        lines.append(f'{s}\t{f1}\t{f2}')
    path.write_text('\n'.join(lines) + '\n')


def _cfg(path):
    return {'file': str(path), 'sep': '\t', 'sample_col': 'sample',
            'drop_cols': [], 'high_dim': False}


def _dup_rows():
    """The 5 real duplicated IDs, each twice; copies carry distinct values
    so the kept row is identifiable."""
    rows = []
    for i, s in enumerate(DUP_SAMPLES):
        rows.append((s, 0.1 * (i + 1), 1.0 * (i + 1)))
        rows.append((s, 99.0, 99.0))  # duplicate copy: must NOT be the kept row
    return rows


def test_dedup_removes_duplicate_sample_rows(tmp_path):
    # Mirror the real bug shape: the 5 duplicated IDs each appear twice.
    f = tmp_path / 'features.tsv'
    _write_fixture(f, _dup_rows())
    meta = DUP_SAMPLES + ['V1_S1']  # V1_S1 in metadata but not in the file

    X, ids, feats = rdp.load_modality(_cfg(f), np.array(meta))

    assert len(ids) == 5                  # deduped unique count, not 10
    assert len(set(ids)) == 5             # no duplicates survive
    assert set(ids) == set(DUP_SAMPLES)
    assert X.shape == (5, 2)
    # keep='first': the kept row is the first occurrence, not the 99.0 copy
    kept = X[np.where(ids == 'V7_S4')[0][0], 0]
    assert kept == pytest.approx(0.4)


def test_dedup_then_intersection_counts_duplicates_once(tmp_path):
    # A duplicated sample present in metadata must count once; a sample not in
    # metadata must be excluded.
    rows = [('A', 0.1, 1.0), ('B', 0.2, 2.0), ('B', 0.2, 2.0),
            ('X', 0.9, 9.0)]  # X not in metadata
    f = tmp_path / 'features.tsv'
    _write_fixture(f, rows)
    meta = ['A', 'B']

    X, ids, feats = rdp.load_modality(_cfg(f), np.array(meta))

    assert sorted(ids.tolist()) == ['A', 'B']
    assert X.shape == (2, 2)


def test_data_loading_cross_check(tmp_path):
    # Cross-check only (ticket #17 wording): the shared loader dedups the same
    # way, so the two loaders agree on the deduped count. Not the primary
    # target — data_loading.py already had the dedup block and was never the
    # bug; a test there alone would give false assurance.
    f = tmp_path / 'features.tsv'
    _write_fixture(f, _dup_rows())

    X, ids, feats = data_loading.load_modality(_cfg(f), np.array(DUP_SAMPLES))

    assert len(ids) == 5
    assert len(set(ids)) == 5
    assert X.shape == (5, 2)
