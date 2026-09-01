"""Version screen (issue #16): loader mechanics + pipeline scope changes.

Covers the ticket's TDD list:
- new MODALITY_CONFIGS entries (13 screen rows registered, correct flags)
- long_format pivot + duplicate sample x probe ValueError
- tt39 probe-level / manifest site-level subsetting
- per-CpG -> probe aggregation (mean of observed sites)
- TOO scope inference in build_summary (TOO prefix checked before EDTA)
- dr='lasso' raw-LASSO path (skips PCA, reports raw feature totals)

Real-data load smoke runs manually per ticket acceptance (the enriched
per-CpG files are gitignored, so committed tests must be hermetic).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import joblib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

import data_loading
import run_phase1_pipeline as pipeline


# The 13 version-screen row names (protocol run matrix, rows 1-13).
SCREEN_ROWS = [
    'probe_meth', 'probe_meth_unenriched', 'probe_meth_unfiltered_qc',
    'probe_cpg', 'probe_cpg_enriched',
    'probe_cpg_agg_unenriched', 'probe_cpg_agg_enriched',
    'probe_meth_tt39_enriched', 'probe_meth_tt39_unenriched',
    'probe_cpg_tt39_unenriched', 'probe_cpg_tt39_unenriched_lasso',
    'probe_cpg_tt39_enriched', 'probe_cpg_tt39_enriched_lasso',
]


# ── fixtures ──────────────────────────────────────────

def _write_long(path, rows):
    """LONG probe_meth-style file (sample, probe_id, CpG_frac + extras)."""
    lines = ['sample\tprobe_id\tCpG_meth\tCpG_total\tCpG_frac']
    for s, p, frac in rows:
        lines.append(f'{s}\t{p}\t10\t20\t{frac}')
    path.write_text('\n'.join(lines) + '\n')


def _write_manifest(path, pairs):
    """Manifest: probe_id + feature_id rows (one per measured site)."""
    lines = ['probe_id\tchr\tpos\tcpg_idx\tfeature_id']
    for probe_id, fid in pairs:
        lines.append(f'{probe_id}\tchr1\t1\t1\t{fid}')
    path.write_text('\n'.join(lines) + '\n')


def _write_wide(path, feats, rows):
    """Wide feature file: sample column + one column per feature."""
    lines = ['\t'.join(['sample'] + list(feats))]
    for s, vals in rows:
        lines.append('\t'.join([s] + ['' if v is None else str(v) for v in vals]))
    path.write_text('\n'.join(lines) + '\n')


def _cfg(**extra):
    base = {'file': 'x', 'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
            'high_dim': False}
    base.update(extra)
    return base


# ── long_format ──────────────────────────────────────

def test_long_format_pivots_to_wide(tmp_path):
    f = tmp_path / 'long.tsv'
    _write_long(f, [('A', 'cg2', 0.25), ('A', 'cg1', 0.5),
                    ('B', 'cg2', 0.1), ('B', 'cg1', 0.75)])
    cfg = _cfg(file=str(f), long_format=True)

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A', 'B']))

    assert set(ids.tolist()) == {'A', 'B'}
    assert sorted(feats.tolist()) == ['cg1', 'cg2']
    a = np.where(ids == 'A')[0][0]
    vals = {p: X[a, i] for i, p in enumerate(feats)}
    assert vals['cg1'] == pytest.approx(0.5)
    assert vals['cg2'] == pytest.approx(0.25)


def test_long_format_raises_on_duplicate_sample_probe(tmp_path):
    f = tmp_path / 'long.tsv'
    _write_long(f, [('A', 'cg1', 0.5), ('A', 'cg1', 0.9), ('B', 'cg1', 0.1)])
    cfg = _cfg(file=str(f), long_format=True)

    with pytest.raises(ValueError, match='[Dd]uplicate'):
        data_loading.load_modality(cfg, np.array(['A', 'B']))


# ── probe-level tt39 subsetting ──────────────────────

def test_probe_subset_restricts_columns(tmp_path):
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['cg1', 'cg2', 'cg3', 'cg4'], [('A', [0.1, 0.2, 0.3, 0.4])])
    cfg = _cfg(file=str(f), probe_subset=['cg4', 'cg2'])

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A']))

    assert X.shape == (1, 2)
    assert sorted(feats.tolist()) == ['cg2', 'cg4']


def test_probe_subset_missing_probe_raises(tmp_path):
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['cg1', 'cg2'], [('A', [0.1, 0.2])])
    cfg = _cfg(file=str(f), probe_subset=['cg1', 'cgX'])

    with pytest.raises(ValueError, match='cgX'):
        data_loading.load_modality(cfg, np.array(['A']))


def test_probe_subset_applies_after_long_pivot(tmp_path):
    f = tmp_path / 'long.tsv'
    _write_long(f, [('A', 'cg1', 0.5), ('A', 'cg2', 0.25),
                    ('B', 'cg1', 0.75), ('B', 'cg2', 0.1)])
    cfg = _cfg(file=str(f), long_format=True, probe_subset=['cg2'])

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A', 'B']))

    assert X.shape == (2, 1)
    assert feats.tolist() == ['cg2']


# ── per-CpG manifest site subsetting ─────────────────

def test_manifest_subset_restricts_sites(tmp_path):
    man = tmp_path / 'manifest.tsv'
    _write_manifest(man, [('P1', 'f1'), ('P1', 'f2'), ('P2', 'f3'), ('P3', 'f4')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3', 'f4'],
                [('A', [0.1, 0.2, 0.3, 0.4]), ('B', [0.5, 0.6, 0.7, 0.8])])
    cfg = _cfg(file=str(f), manifest=str(man), probe_subset=['P1', 'P2'])

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A', 'B']))

    assert X.shape == (2, 3)
    assert sorted(feats.tolist()) == ['f1', 'f2', 'f3']


def test_manifest_subset_missing_probe_raises(tmp_path):
    man = tmp_path / 'manifest.tsv'
    _write_manifest(man, [('P1', 'f1')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1'], [('A', [0.1])])
    cfg = _cfg(file=str(f), manifest=str(man), probe_subset=['P1', 'P2'])

    with pytest.raises(ValueError, match='P2'):
        data_loading.load_modality(cfg, np.array(['A']))


# ── per-CpG -> probe aggregation ─────────────────────

def test_agg_by_probe_means_observed_sites(tmp_path):
    man = tmp_path / 'manifest.tsv'
    _write_manifest(man, [('P1', 'f1'), ('P1', 'f2'), ('P2', 'f3')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3'],
                [('A', [0.4, 0.8, None]), ('B', [None, None, None]),
                 ('C', [0.2, 0.6, None])])
    cfg = _cfg(file=str(f), manifest=str(man))

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A', 'B', 'C']),
                                               impute=False)

    # P2 is all-NaN across samples -> dropped by the all-NaN feature step
    assert sorted(feats.tolist()) == ['P1']
    a = np.where(ids == 'A')[0][0]
    assert X[a, 0] == pytest.approx(0.6)  # mean of observed sites (0.4+0.8)/2
    b = np.where(ids == 'B')[0][0]
    assert np.isnan(X[b, 0])


def test_agg_by_probe_ignores_unobserved_sites_per_sample(tmp_path):
    man = tmp_path / 'manifest.tsv'
    _write_manifest(man, [('P1', 'f1'), ('P1', 'f2'), ('P2', 'f3')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3'],
                [('A', [0.4, 0.8, 0.2]), ('B', [None, 0.6, 0.5])])
    cfg = _cfg(file=str(f), manifest=str(man))

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A', 'B']))

    b = np.where(ids == 'B')[0][0]
    vals = {p: X[b, i] for i, p in enumerate(feats)}
    assert vals['P1'] == pytest.approx(0.6)  # single observed site, not 0.3
    assert vals['P2'] == pytest.approx(0.5)


# ── config registry (13 screen rows) ─────────────────

def test_screen_rows_are_registered():
    for name in SCREEN_ROWS:
        assert name in data_loading.MODALITY_CONFIGS, name
    assert list(pipeline.PHASE1_MODALITIES) == SCREEN_ROWS


def test_screen_config_flags():
    tt39 = data_loading.TT39_PROBES
    assert len(tt39) == 39

    tt39_rows = ('probe_meth_tt39_enriched', 'probe_meth_tt39_unenriched',
                 'probe_cpg_tt39_unenriched', 'probe_cpg_tt39_enriched',
                 'probe_cpg_tt39_unenriched_lasso', 'probe_cpg_tt39_enriched_lasso')
    for name in tt39_rows:
        assert data_loading.MODALITY_CONFIGS[name]['probe_subset'] == tt39

    long_rows = ('probe_meth_unenriched', 'probe_meth_unfiltered_qc',
                 'probe_meth_tt39_unenriched')
    for name in long_rows:
        assert data_loading.MODALITY_CONFIGS[name]['long_format'] is True

    for name in ('probe_cpg_tt39_unenriched_lasso', 'probe_cpg_tt39_enriched_lasso'):
        cfg = data_loading.MODALITY_CONFIGS[name]
        assert cfg['dr'] == 'lasso'
        assert not pipeline.should_use_pca(cfg)

    for name in ('probe_cpg_tt39_unenriched', 'probe_cpg_tt39_enriched',
                 'probe_cpg', 'probe_cpg_enriched'):
        assert pipeline.should_use_pca(data_loading.MODALITY_CONFIGS[name])

    for name in ('probe_cpg_agg_unenriched', 'probe_cpg_agg_enriched'):
        cfg = data_loading.MODALITY_CONFIGS[name]
        assert cfg['manifest'] is not None
        assert cfg['high_dim'] is False


def test_should_use_pca():
    assert pipeline.should_use_pca({'high_dim': True}) is True
    assert pipeline.should_use_pca({'high_dim': True, 'dr': 'lasso'}) is False
    assert pipeline.should_use_pca({'high_dim': False, 'dr': 'lasso'}) is False
    assert pipeline.should_use_pca({'high_dim': False}) is False


# ── scope inference ──────────────────────────────────

def _fake_result(label, **over):
    r = {'label': label, 'macro_f1_mean': 0.5, 'macro_f1_std': 0.1,
         'median_C': 1.0, 'n_nonzero_full': 5, 'n_pcs': None,
         'used_pca': False, 'n_total_features': 148}
    r.update(over)
    return r


def test_build_summary_too_prefix_wins_over_edta():
    results = [
        _fake_result('Methylation (probe-avg)'),
        _fake_result('TOO Methylation (probe-avg)'),
        _fake_result('TOO EDTA Methylation (probe-avg)'),
    ]

    df = pipeline.build_summary(results)

    scopes = dict(zip(df['Modality'], df['Scope']))
    assert scopes['Methylation (probe-avg)'] == 'Full'
    assert scopes['TOO Methylation (probe-avg)'] == 'TOO'
    assert scopes['TOO EDTA Methylation (probe-avg)'] == 'TOO EDTA'


def test_build_summary_totals_for_pca_vs_lasso_rows():
    pca = _fake_result('TOO Per-CpG methylation (enriched)', used_pca=True,
                       n_pcs=20, n_total_features=54300)
    lasso = _fake_result('TOO Per-CpG methylation (tt39, raw LASSO)',
                         used_pca=False, n_pcs=None, n_total_features=6112)

    df = pipeline.build_summary([pca, lasso])

    totals = dict(zip(df['Modality'], df['total_features']))
    assert totals['TOO Per-CpG methylation (enriched)'] == 20
    assert totals['TOO Per-CpG methylation (tt39, raw LASSO)'] == 6112
    assert df.loc[df['Modality'] == 'TOO Per-CpG methylation (enriched)',
                  'n_original_features'].iloc[0] == 54300
    assert pd.isna(df.loc[df['Modality'] == 'TOO Per-CpG methylation (tt39, raw LASSO)',
                          'n_original_features'].iloc[0])


# ── scopes (get_scope) ───────────────────────────────

def test_get_scope_subsets():
    meta = pd.DataFrame({
        'TWIST_ID': ['s1', 's2', 's3', 's4', 's5', 's6'],
        'Tissue': ['colon', 'colon', 'healthyblood', 'pancreas', 'pancreas', 'pancreas'],
        'BCT': ['EDTA', 'Citrate', 'EDTA', 'EDTA', 'Streck', 'Streck'],
    })

    full, prefix = pipeline.get_scope('full', meta)
    assert prefix == ''
    assert len(full) == 6

    too, prefix = pipeline.get_scope('too', meta)
    assert prefix == 'TOO '
    assert set(too['TWIST_ID']) == {'s1', 's2', 's4', 's5', 's6'}

    too_edta, prefix = pipeline.get_scope('too-edta', meta)
    assert prefix == 'TOO EDTA '
    assert set(too_edta['TWIST_ID']) == {'s1', 's4'}

    with pytest.raises(ValueError):
        pipeline.get_scope('bogus', meta)


# ── raw-LASSO pipeline path ──────────────────────────

def test_lasso_dr_runs_without_pca(tmp_path, monkeypatch):
    out = tmp_path / 'out'
    figs = out / 'figures'
    models = out / 'models'
    monkeypatch.setattr(pipeline, 'OUT', out)
    monkeypatch.setattr(pipeline, 'FIGS', figs)
    monkeypatch.setattr(pipeline, 'MODELS', models)
    figs.mkdir(parents=True)
    models.mkdir(parents=True)

    f = tmp_path / 'wide.tsv'
    feats = [f'f{i}' for i in range(10)]
    rows = []
    for i in range(12):
        vals = ['' if (j % 3 == 0 and i % 2 == 0) else f'{(i * 10 + j) % 7 / 10:.2f}'
                for j in range(10)]
        rows.append((f's{i + 1}', vals))
    _write_wide(f, feats, rows)

    cfg = _cfg(file=str(f), high_dim=True, dr='lasso', label='Fake (raw LASSO)')
    monkeypatch.setitem(data_loading.MODALITY_CONFIGS, 'fake_lasso', cfg)

    meta = pd.DataFrame({
        'TWIST_ID': [f's{i + 1}' for i in range(12)],
        'Tissue': ['colon'] * 6 + ['liver'] * 6,
        'Source': ['A', 'B'] * 6,
    })

    r = pipeline.run_modality_pipeline('fake_lasso', meta)

    assert r is not None
    assert r['used_pca'] is False
    assert r['n_pcs'] is None
    assert r['n_total_features'] == 10

    art = joblib.load(models / 'fake_lasso_full_model.joblib')
    assert 'pca' not in art
    assert 'scaler' in art
    assert list(art['feature_names']) == feats
