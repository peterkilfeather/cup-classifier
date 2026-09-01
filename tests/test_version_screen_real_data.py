"""Real-data smoke for the 13 version-screen rows (issue #16 acceptance).

Asserts every screen config loads on the real feature files with the full
164-sample metadata set. Skips when the enriched per-CpG files are absent
(they are gitignored, so a fresh clone cannot run these).

Feature-count expectations:
- probe-avg rows: exactly 148 probes (file-structural).
- tt39 probe-avg rows: exactly 39 probes (panel-structural).
- aggregation rows: exactly 148 probes (manifest-structural, unless a probe
  is all-NaN across all 164 samples — currently none are).
- per-CpG rows: site counts vary with all-NaN column drops over the sample
  set, so assert protocol-consistent bounds instead of exact values.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

import data_loading

ENRICHED_WIDE = data_loading.PROBE_CPG_DIR / (
    'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv')

SCREEN_ROWS = [
    'probe_meth', 'probe_meth_unenriched', 'probe_meth_unfiltered_qc',
    'probe_cpg', 'probe_cpg_enriched',
    'probe_cpg_agg_unenriched', 'probe_cpg_agg_enriched',
    'probe_meth_tt39_enriched', 'probe_meth_tt39_unenriched',
    'probe_cpg_tt39_unenriched', 'probe_cpg_tt39_unenriched_lasso',
    'probe_cpg_tt39_enriched', 'probe_cpg_tt39_enriched_lasso',
]

EXACT_FEATURES = {
    'probe_meth': 148, 'probe_meth_unenriched': 148,
    'probe_meth_unfiltered_qc': 148,
    'probe_cpg_agg_unenriched': 148, 'probe_cpg_agg_enriched': 148,
    'probe_meth_tt39_enriched': 39, 'probe_meth_tt39_unenriched': 39,
}

# (lo, hi) bounds on per-CpG site counts after all-NaN drops (protocol's
# "~5-8K tt39 sites" estimate is stale; actuals are 0.7-1.4K).
BOUNDS_FEATURES = {
    'probe_cpg': (30000, 33000),
    'probe_cpg_enriched': (50000, 57000),
    'probe_cpg_tt39_unenriched': (700, 1200),
    'probe_cpg_tt39_unenriched_lasso': (700, 1200),
    'probe_cpg_tt39_enriched': (1200, 2000),
    'probe_cpg_tt39_enriched_lasso': (1200, 2000),
}


@pytest.mark.skipif(not ENRICHED_WIDE.exists(),
                    reason='enriched per-CpG files not present (gitignored)')
def test_all_13_screen_rows_load_on_real_data():
    meta = data_loading.load_metadata()
    ids = meta['TWIST_ID'].values

    for name in SCREEN_ROWS:
        cfg = data_loading.MODALITY_CONFIGS[name]
        X, sids, feats = data_loading.load_modality(cfg, ids, impute=False)
        assert len(sids) == len(meta), f'{name}: sample count {len(sids)} != {len(meta)}'
        assert X.shape == (len(sids), len(feats)), name
        if name in EXACT_FEATURES:
            assert len(feats) == EXACT_FEATURES[name], name
        else:
            lo, hi = BOUNDS_FEATURES[name]
            assert lo <= len(feats) <= hi, f'{name}: {len(feats)} features'


def test_scope_counts_on_real_metadata():
    meta = data_loading.load_metadata()
    too = meta[meta['Tissue'] != 'healthyblood']
    too_edta = meta[(meta['BCT'] == 'EDTA') & (meta['Tissue'] != 'healthyblood')]

    assert len(too) == 124
    assert set(too['Tissue'].unique()) == {'colon', 'liver', 'pancreas',
                                           'prostate', 'stomach'}
    assert len(too_edta) == 56
    assert 'pancreas' not in set(too_edta['Tissue'].unique())
