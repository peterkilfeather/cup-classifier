"""Version screen addendum (issue #18): in-probe-only rows (mechanics).

Covers the ticket's TDD list:
- chr normalization in the manifest filter (manifest 'chr1' vs annotations '1')
- window edge: sites at exactly ±half_width included (boundary inclusive)
- zero-site probe drop (no in-window sites -> no column), parameterized
  over the per-CpG and aggregation paths
- half_width config key honored at runtime (default 61)
- in-probe aggregation means only in-window sites
- raw-LASSO pipeline path for the in-probe per-CpG rows
- tt39 in-window guard: all 39 tt39 probes sit fully in-window under the
  default window (real tracked manifests + annotations; hermetic)
- config registry flags for the 4 addendum rows

Skip-guarded real-data smoke (enriched per-CpG files are gitignored):
- the 4 addendum configs load with expected feature counts
- under the default window cg14861089 has zero in-probe sites and drops
  from the aggregation columns
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

ENRICHED_WIDE = data_loading.PROBE_CPG_DIR / (
    'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv')

# The 4 addendum row names (protocol run matrix, rows 14-17).
ADDENDUM_ROWS = [
    'probe_cpg_inprobe_unenriched', 'probe_cpg_inprobe_enriched',
    'probe_cpg_agg_inprobe_unenriched', 'probe_cpg_agg_inprobe_enriched',
]


# ── fixtures ──────────────────────────────────────────

def _write_annotations(path, rows):
    """probe_annotations_450k-style CSV (probe_id,chr,mapinfo + extras).

    chr uses the annotations naming (no 'chr' prefix); mapinfo is the 450K
    target position of the probe.
    """
    lines = ['probe_id,chr,mapinfo,ucsc_refgene_name,ucsc_refgene_group,relation_to_island']
    for probe_id, chr_, mapinfo in rows:
        lines.append(f'{probe_id},{chr_},{mapinfo},,,')
    path.write_text('\n'.join(lines) + '\n')


def _write_manifest(path, rows):
    """Manifest: probe_id + chr + pos + feature_id (chr uses manifest
    naming — 'chr1' style)."""
    lines = ['probe_id\tchr\tpos\tcpg_idx\tfeature_id']
    for probe_id, chr_, pos, fid in rows:
        lines.append(f'{probe_id}\t{chr_}\t{pos}\t1\t{fid}')
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


def _in_probe_env(tmp_path, monkeypatch, ann_rows, man_rows):
    """Write synthetic annotations + manifest + wide file; point the module
    at the synthetic annotations. Returns the manifest path."""
    ann = tmp_path / 'annotations.csv'
    _write_annotations(ann, ann_rows)
    monkeypatch.setattr(data_loading, 'PROBE_ANNOTATIONS', ann)
    man = tmp_path / 'manifest.tsv'
    _write_manifest(man, man_rows)
    return man


# ── in-probe site filter (manifest x annotations join) ─

def test_in_probe_filter_normalizes_chr_naming(tmp_path, monkeypatch):
    # manifest chr 'chr1' must match annotations chr '1'; without
    # normalization the join produces no mapinfo and the loader raises.
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P1', 'chr1', 100500, 'f2')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2'], [('A', [0.1, 0.9])])
    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True)

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A']))

    assert X.shape == (1, 1)
    assert feats.tolist() == ['f1']


def test_in_probe_window_boundary_inclusive(tmp_path, monkeypatch):
    # Sites exactly at ±half_width are in; one base pair beyond is out.
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000)],
                        man_rows=[('P1', 'chr1', 100000 - 61, 'f1'),
                                  ('P1', 'chr1', 100000, 'f2'),
                                  ('P1', 'chr1', 100000 + 61, 'f3'),
                                  ('P1', 'chr1', 100000 + 62, 'f4')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3', 'f4'], [('A', [0.1, 0.2, 0.3, 0.4])])
    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True)

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A']))

    assert sorted(feats.tolist()) == ['f1', 'f2', 'f3']


def test_half_width_config_key_controls_membership(tmp_path, monkeypatch):
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P1', 'chr1', 100020, 'f2'),
                                  ('P1', 'chr1', 100100, 'f3')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3'], [('A', [0.1, 0.2, 0.3])])

    _, _, tight = data_loading.load_modality(
        _cfg(file=str(f), manifest=str(man), in_probe_only=True,
             half_width=10), np.array(['A']))
    _, _, loose = data_loading.load_modality(
        _cfg(file=str(f), manifest=str(man), in_probe_only=True,
             half_width=30), np.array(['A']))

    assert tight.tolist() == ['f1']
    assert loose.tolist() == ['f1', 'f2']


def test_in_probe_manifest_probe_missing_from_annotations_raises(tmp_path, monkeypatch):
    # Only P1 is annotated; the manifest also carries P2 -> the join cannot
    # compute membership and must fail loudly rather than silently drop.
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P2', 'chr1', 500000, 'f2')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2'], [('A', [0.1, 0.2])])
    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True)

    with pytest.raises(ValueError, match='P2'):
        data_loading.load_modality(cfg, np.array(['A']))


# ── zero-site probe drop (parameterized over both paths) ─

@pytest.mark.parametrize('agg', [False, True],
                         ids=['per-cpg', 'aggregated'])
def test_in_probe_zero_site_probe_drops(tmp_path, monkeypatch, agg):
    # P2's only site lies 1000 bp from its mapinfo: no in-window sites.
    # Per-CpG path: no feature column survives the usecols filter.
    # Aggregation path: the probe never enters probe_cols -> drops as a
    # consequence of the filter, with no hardcoded probe names.
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000), ('P2', 1, 900000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P2', 'chr1', 901000, 'f2')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2'], [('A', [0.2, 0.8])])
    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True,
               high_dim=agg is False, dr='lasso' if not agg else None,
               aggregate=True if agg else None)

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A']))

    # Per-CpG path keeps site columns; aggregation path emits probe columns.
    # Either way P2 (no in-window sites) yields nothing.
    expected = ['P1'] if agg else ['f1']
    assert X.shape == (1, 1)
    assert feats.tolist() == expected


# ── in-probe aggregation ─────────────────────────────

def test_in_probe_agg_means_only_window_sites(tmp_path, monkeypatch):
    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P1', 'chr1', 100061, 'f2'),
                                  ('P1', 'chr1', 100200, 'f3')])
    f = tmp_path / 'wide.tsv'
    _write_wide(f, ['f1', 'f2', 'f3'], [('A', [0.2, 0.8, 0.9])])
    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True,
               aggregate=True)

    X, ids, feats = data_loading.load_modality(cfg, np.array(['A']))

    assert feats.tolist() == ['P1']
    a = np.where(ids == 'A')[0][0]
    assert X[a, 0] == pytest.approx(0.5)  # (0.2+0.8)/2; flanking 0.9 excluded


# ── config registry (4 addendum rows) ────────────────

def test_addendum_rows_are_registered():
    for name in ADDENDUM_ROWS:
        assert name in data_loading.MODALITY_CONFIGS, name
    assert data_loading.DEFAULT_HALF_WIDTH == 61


def test_addendum_config_flags():
    for name in ('probe_cpg_inprobe_unenriched', 'probe_cpg_inprobe_enriched'):
        cfg = data_loading.MODALITY_CONFIGS[name]
        assert cfg['in_probe_only'] is True
        assert cfg['high_dim'] is True
        assert cfg['dr'] == 'lasso'
        assert cfg['manifest'] is not None
        assert not pipeline.should_use_pca(cfg)

    for name in ('probe_cpg_agg_inprobe_unenriched', 'probe_cpg_agg_inprobe_enriched'):
        cfg = data_loading.MODALITY_CONFIGS[name]
        assert cfg['in_probe_only'] is True
        assert cfg['high_dim'] is False
        assert cfg['aggregate'] is True
        assert cfg['manifest'] is not None
        assert 'dr' not in cfg


# ── raw-LASSO pipeline path on an in-probe row ───────

def test_in_probe_lasso_runs_without_pca(tmp_path, monkeypatch):
    out = tmp_path / 'out'
    figs = out / 'figures'
    models = out / 'models'
    monkeypatch.setattr(pipeline, 'OUT', out)
    monkeypatch.setattr(pipeline, 'FIGS', figs)
    monkeypatch.setattr(pipeline, 'MODELS', models)
    figs.mkdir(parents=True)
    models.mkdir(parents=True)

    man = _in_probe_env(tmp_path, monkeypatch,
                        ann_rows=[('P1', 1, 100000), ('P2', 1, 900000)],
                        man_rows=[('P1', 'chr1', 100000, 'f1'),
                                  ('P1', 'chr1', 100040, 'f2'),
                                  ('P2', 'chr1', 901000, 'f3')])
    f = tmp_path / 'wide.tsv'
    feats = ['f1', 'f2', 'f3']
    rows = []
    for i in range(12):
        vals = [f'{(i * 10 + j) % 7 / 10:.2f}' for j in range(3)]
        rows.append((f's{i + 1}', vals))
    _write_wide(f, feats, rows)

    cfg = _cfg(file=str(f), manifest=str(man), in_probe_only=True,
               high_dim=True, dr='lasso', label='Fake (in-probe raw LASSO)')
    monkeypatch.setitem(data_loading.MODALITY_CONFIGS, 'fake_inprobe_lasso', cfg)

    meta = pd.DataFrame({
        'TWIST_ID': [f's{i + 1}' for i in range(12)],
        'Tissue': ['colon'] * 6 + ['liver'] * 6,
        'Source': ['A', 'B'] * 6,
    })

    r = pipeline.run_modality_pipeline('fake_inprobe_lasso', meta)

    assert r is not None
    assert r['used_pca'] is False
    assert r['n_pcs'] is None
    assert r['n_total_features'] == 2  # in-probe sites only (f1, f2)

    art = joblib.load(models / 'fake_inprobe_lasso_full_model.joblib')
    assert 'pca' not in art
    assert list(art['feature_names']) == ['f1', 'f2']


# ── tt39 in-window guard (real tracked inputs, hermetic) ─

def test_tt39_probes_fully_in_window():
    """Guard: all 39 tt39 probes sit fully in-probe under the default window
    (verified 2026-09-08: max site distance 61 bp, both captures). Re-check
    under any future half_width change."""
    half_width = data_loading.DEFAULT_HALF_WIDTH
    for cap in ('unenriched', 'enriched'):
        man = pd.read_csv(data_loading.PROBE_CPG_DIR
                          / f'all_samples.probe_cpg_{cap}_filtered.mNonCpGlt4.manifest.tsv',
                          sep='\t')
        in_feats = set(data_loading._in_probe_feature_ids(man, half_width))
        for probe in data_loading.TT39_PROBES:
            sites = man.loc[man['probe_id'] == probe, 'feature_id']
            assert len(sites) > 0, f'{cap}: {probe} missing from manifest'
            assert set(sites).issubset(in_feats), f'{cap}: {probe} has flanking sites'


# ── real-data smoke (skip-guarded) ───────────────────

@pytest.mark.skipif(not ENRICHED_WIDE.exists(),
                    reason='enriched per-CpG files not present (gitignored)')
def test_addendum_rows_load_on_real_data():
    """4 addendum configs load on the real feature files with expected
    feature counts (bounds absorb all-NaN feature drops over the 164-sample
    universe; exact memberships are window-dependent, derived at runtime)."""
    meta = data_loading.load_metadata()
    ids = meta['TWIST_ID'].values
    expected = {
        'probe_cpg_inprobe_unenriched': (3500, 3560),
        'probe_cpg_inprobe_enriched': (5700, 5900),
        'probe_cpg_agg_inprobe_unenriched': (147, 147),
        'probe_cpg_agg_inprobe_enriched': (147, 147),
    }
    for name, (lo, hi) in expected.items():
        cfg = data_loading.MODALITY_CONFIGS[name]
        X, sids, feats = data_loading.load_modality(cfg, ids, impute=False)
        assert len(sids) == len(meta), f'{name}: sample count {len(sids)} != {len(meta)}'
        assert X.shape == (len(sids), len(feats)), name
        assert lo <= len(feats) <= hi, f'{name}: {len(feats)} features'


@pytest.mark.skipif(not ENRICHED_WIDE.exists(),
                    reason='enriched per-CpG files not present (gitignored)')
def test_cg14861089_zero_in_probe_drops_on_real_data():
    """Under the default window cg14861089 has zero in-probe sites (its 29/49
    measured sites are all flanking) and drops from the aggregation columns."""
    meta = data_loading.load_metadata()
    ids = meta['TWIST_ID'].values
    for name in ('probe_cpg_agg_inprobe_unenriched', 'probe_cpg_agg_inprobe_enriched'):
        X, sids, feats = data_loading.load_modality(
            data_loading.MODALITY_CONFIGS[name], ids, impute=False)
        assert 'cg14861089' not in feats.tolist(), name
