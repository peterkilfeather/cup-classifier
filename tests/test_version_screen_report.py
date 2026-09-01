"""Version screen report builder (ticket #12): tests for the reporting step.

The reporting script reads run outputs from output/phase1/ (summary CSV,
per-run cv_metrics CSVs, hyperparameters JSON) and computes the hand-back
table: per-scope macro-F1 ± std, n samples, n features, median C, missingness,
plus per-fold paired deltas vs the same-capture probe-avg row (TOO scopes).

Hermetic: no real data; synthetic outputs via tmp_path, loader stubbed.
"""
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))

import report_version_screen as rep
import run_phase1_pipeline as pipeline


# ── fixtures ──────────────────────────────────────────

def _write_summary(path, rows):
    """Synthetic phase1_summary.csv (columns as built by build_summary)."""
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


def _write_cv_metrics(path, folds, values):
    pd.DataFrame({'fold': folds, 'macro_f1': values,
                  'balanced_accuracy': values,
                  'src_acc_SrcA': [0.5] * len(folds),
                  'src_acc_SrcB': [0.6] * len(folds)}).to_csv(path, index=False)


def _write_hp(path, n_samples, n_features, median_c):
    (path).write_text(json.dumps({'n_samples': n_samples,
                                  'n_features': n_features,
                                  'n_pcs': 20 if n_features == 20 else None,
                                  'median_C': median_c}))


def _synthetic_out(tmp_path):
    """Output dir with all 13 rows x 2 scopes of synthetic run artifacts."""
    out = tmp_path / 'phase1'
    out.mkdir()
    labels = {
        'probe_meth': 'Methylation (probe-avg)',
        'probe_meth_unenriched': 'Methylation (probe-avg, unenriched)',
        'probe_meth_unfiltered_qc': 'Methylation (probe-avg, unfiltered QC)',
        'probe_cpg': 'Per-CpG methylation (32K->PCs)',
        'probe_cpg_enriched': 'Per-CpG methylation (enriched)',
        'probe_cpg_agg_unenriched':
            'Per-CpG methylation (aggregated to probes, unenriched)',
        'probe_cpg_agg_enriched':
            'Per-CpG methylation (aggregated to probes, enriched)',
        'probe_meth_tt39_enriched': 'Methylation (probe-avg, tt39, enriched)',
        'probe_meth_tt39_unenriched': 'Methylation (probe-avg, tt39, unenriched)',
        'probe_cpg_tt39_unenriched': 'Per-CpG methylation (tt39, unenriched)',
        'probe_cpg_tt39_unenriched_lasso':
            'Per-CpG methylation (tt39, unenriched, raw LASSO)',
        'probe_cpg_tt39_enriched': 'Per-CpG methylation (tt39, enriched)',
        'probe_cpg_tt39_enriched_lasso':
            'Per-CpG methylation (tt39, enriched, raw LASSO)',
    }
    n_feats = {'probe_meth': 148, 'probe_meth_unenriched': 148,
               'probe_meth_unfiltered_qc': 148, 'probe_cpg': 20,
               'probe_cpg_enriched': 20, 'probe_cpg_agg_unenriched': 148,
               'probe_cpg_agg_enriched': 148, 'probe_meth_tt39_enriched': 39,
               'probe_meth_tt39_unenriched': 39, 'probe_cpg_tt39_unenriched': 20,
               'probe_cpg_tt39_unenriched_lasso': 900,
               'probe_cpg_tt39_enriched': 20, 'probe_cpg_tt39_enriched_lasso': 1300}
    # Distinct base macro-F1 per modality so paired deltas are non-trivial;
    # probe_cpg (0.70) vs its anchor probe_meth_unenriched (0.60) = +0.1.
    base_f1 = {m: 0.55 + 0.01 * i for i, m in enumerate(pipeline.PHASE1_MODALITIES)}
    base_f1['probe_meth'] = 0.60
    base_f1['probe_meth_unenriched'] = 0.60
    base_f1['probe_cpg'] = 0.70
    summary_rows = []
    for scope, prefix in [('Full', ''), ('TOO', 'TOO ')]:
        for m in pipeline.PHASE1_MODALITIES:
            tag = 'TOO_' + m if scope == 'TOO' else m
            f1 = base_f1[m]
            _write_cv_metrics(out / f'{tag}_cv_metrics.csv',
                              [1, 2, 3, 4, 5], [f1 + i * 0.01 for i in range(5)])
            _write_hp(out / f'{tag}_hyperparameters.json',
                      124 if scope == 'TOO' else 164, n_feats[m], 1.0)
            summary_rows.append({
                'Scope': scope, 'Modality': prefix + labels[m],
                'macro_F1_mean': round(f1, 4), 'macro_F1_std': 0.05,
                'macro_F1_CI95_lower': round(f1 - 0.1, 4),
                'macro_F1_CI95_upper': round(f1 + 0.1, 4),
                'median_C': 1.0,
                'selected_features': 140, 'total_features': n_feats[m],
            })
    _write_summary(out / 'phase1_summary.csv', summary_rows)
    return out


# ── anchor mapping (protocol: pair same captures) ─────

def test_same_capture_anchor_map():
    assert rep.anchor_for('probe_cpg') == 'probe_meth_unenriched'
    assert rep.anchor_for('probe_cpg_enriched') == 'probe_meth'
    assert rep.anchor_for('probe_cpg_agg_unenriched') == 'probe_meth_unenriched'
    assert rep.anchor_for('probe_cpg_agg_enriched') == 'probe_meth'
    assert rep.anchor_for('probe_cpg_tt39_unenriched') == 'probe_meth_tt39_unenriched'
    assert rep.anchor_for('probe_cpg_tt39_unenriched_lasso') == 'probe_meth_tt39_unenriched'
    assert rep.anchor_for('probe_cpg_tt39_enriched') == 'probe_meth_tt39_enriched'
    assert rep.anchor_for('probe_cpg_tt39_enriched_lasso') == 'probe_meth_tt39_enriched'
    # probe-avg rows are anchors themselves, not compared to themselves
    for m in ['probe_meth', 'probe_meth_unenriched', 'probe_meth_unfiltered_qc',
              'probe_meth_tt39_enriched', 'probe_meth_tt39_unenriched']:
        assert rep.anchor_for(m) is None


def test_anchor_map_covers_all_screen_rows():
    for m in pipeline.PHASE1_MODALITIES:
        if m in rep.SAME_CAPTURE_ANCHOR:
            assert rep.anchor_for(m) in pipeline.PHASE1_MODALITIES


def test_aggregated_anchor_same_data_only():
    # Protocol's clean granularity test pairs per-CpG vs its OWN aggregation.
    # tt39 rows (39-probe subsets) have no aggregated counterpart in the run
    # matrix — pairing them against the full 148-probe agg rows would mix
    # probe-set with granularity.
    assert rep.AGGREGATED_ANCHOR == {
        'probe_cpg': 'probe_cpg_agg_unenriched',
        'probe_cpg_enriched': 'probe_cpg_agg_enriched',
    }
    for m in pipeline.PHASE1_MODALITIES:
        if m not in rep.AGGREGATED_ANCHOR:
            continue
        anchor = rep.AGGREGATED_ANCHOR[m]
        # same capture: unenriched pairs unenriched, enriched pairs enriched
        assert ('_enriched' in m) == ('_enriched' in anchor)


# ── tag convention (mirrors pipeline) ─────────────────

def test_tag_for_matches_pipeline_prefixes():
    assert rep.tag_for('probe_cpg', 'full') == 'probe_cpg'
    assert rep.tag_for('probe_cpg', 'too') == 'TOO_probe_cpg'
    assert rep.tag_for('probe_cpg', 'too-edta') == 'TOO_EDTA_probe_cpg'


# ── paired per-fold deltas ────────────────────────────

def test_paired_fold_deltas():
    mod = pd.DataFrame({'fold': [1, 2, 3], 'macro_f1': [0.5, 0.6, 0.55]})
    anchor = pd.DataFrame({'fold': [1, 2, 3], 'macro_f1': [0.4, 0.45, 0.6]})
    d = rep.paired_fold_deltas(mod, anchor)
    assert d['delta'].tolist() == pytest.approx([0.1, 0.15, -0.05])
    assert d['mean_delta'].iloc[0] == pytest.approx((0.1 + 0.15 - 0.05) / 3)
    assert d['n_folds_positive'].iloc[0] == 2


def test_paired_fold_deltas_joins_on_fold():
    mod = pd.DataFrame({'fold': [2, 3], 'macro_f1': [0.6, 0.55]})
    anchor = pd.DataFrame({'fold': [1, 2, 3], 'macro_f1': [0.4, 0.45, 0.6]})
    d = rep.paired_fold_deltas(mod, anchor)
    assert d['fold'].tolist() == [2, 3]


# ── missingness stats ─────────────────────────────────

def test_missingness_stats():
    X = np.array([[np.nan, 1.0, np.nan],
                  [2.0, np.nan, 3.0],
                  [4.0, 5.0, 6.0]])
    cell_frac, sample_frac = rep.missingness_stats(X)
    assert cell_frac == pytest.approx(3 / 9)
    assert sample_frac == pytest.approx(2 / 3)


def test_missingness_stats_no_missing():
    X = np.ones((4, 5))
    assert rep.missingness_stats(X) == (0.0, 0.0)


# ── build_report end-to-end (synthetic outputs) ───────

def test_build_report_from_outputs(tmp_path):
    out = _synthetic_out(tmp_path)

    def fake_missing(mod_name, scope):
        return (0.02, 0.10) if mod_name == 'probe_cpg' else (0.0, 0.0)

    summary, deltas = rep.build_report(out, missingness_fn=fake_missing)

    # Per-scope table: 2 scopes x 13 rows
    assert len(summary) == 26
    too = summary[summary['scope'] == 'TOO']
    assert len(too) == 13
    row = too[too['modality'] == 'probe_cpg'].iloc[0]
    assert row['macro_F1_mean'] == pytest.approx(0.7)
    assert row['macro_F1_std'] == pytest.approx(0.05)
    assert row['macro_F1_CI95_lower'] == pytest.approx(0.6)
    assert row['macro_F1_CI95_upper'] == pytest.approx(0.8)
    assert row['balanced_accuracy'] == pytest.approx(0.72)  # mean(0.70..0.74)
    assert row['n_samples'] == 124
    assert row['n_features'] == 20
    assert row['median_C'] == 1.0
    assert row['missing_cell_frac'] == pytest.approx(0.02)
    assert row['src_acc_SrcA'] == pytest.approx(0.5)
    assert row['src_acc_SrcB'] == pytest.approx(0.6)

    # Deltas: TOO scope only; 8 pairings x probe-averaged anchor (all per-CpG
    # rows incl. aggregated) + 2 same-data probe-aggregation pairings, one
    # row per fold (5). Full scope never compared. tt39 rows pair ONLY
    # against their probe-averaged tt39 anchor (no same-data agg row exists).
    assert len(deltas) == 50
    assert (deltas['scope'] == 'TOO').all()
    assert (deltas['anchor_type'].value_counts().to_dict()
            == {'probe-averaged': 40, 'probe-aggregation': 10})
    d = deltas[(deltas['modality'] == 'probe_cpg')
               & (deltas['anchor_type'] == 'probe-averaged')].iloc[0]
    assert d['anchor'] == 'probe_meth_unenriched'
    assert d['n_folds_positive'] == 5
    assert d['mean_delta'] == pytest.approx(0.1)


def test_build_report_missing_anchor_skips_deltas(tmp_path):
    out = _synthetic_out(tmp_path)
    # Drop the anchor's cv metrics: that pairing is skipped, not crashed.
    # probe_meth anchors the enriched rows (probe_cpg_enriched and
    # probe_cpg_agg_enriched, probe-averaged pairings) -> 2 pairings dropped.
    (out / 'TOO_probe_meth_cv_metrics.csv').unlink()
    summary, deltas = rep.build_report(out, missingness_fn=lambda m, s: (0.0, 0.0))
    # probe_meth's own TOO row is skipped too (no cv_metrics)
    assert len(summary) == 25
    assert len(deltas) == 40
    assert not ((deltas['modality'] == 'probe_cpg_enriched')
                & (deltas['anchor_type'] == 'probe-averaged')).any()
    # probe_cpg's own probe-averaged anchor (unenriched) is untouched
    assert ((deltas['modality'] == 'probe_cpg')
            & (deltas['anchor_type'] == 'probe-averaged')).any()


def test_build_report_missing_own_metrics_skips(tmp_path):
    out = _synthetic_out(tmp_path)
    # Missing the modality's OWN cv_metrics (but not its anchor's): skip the
    # row, do not raise
    (out / 'TOO_probe_cpg_cv_metrics.csv').unlink()
    summary, deltas = rep.build_report(out, missingness_fn=lambda m, s: (0.0, 0.0))
    assert len(summary) == 25
    assert 'probe_cpg' not in summary[summary['scope'] == 'TOO']['modality'].values
    assert not ((deltas['modality'] == 'probe_cpg')
                & (deltas['anchor_type'] == 'probe-averaged')).any()


def test_build_report_requires_summary_columns(tmp_path):
    out = _synthetic_out(tmp_path)
    summary_path = out / 'phase1_summary.csv'
    df = pd.read_csv(summary_path).drop(columns=['median_C'])
    df.to_csv(summary_path, index=False)
    with pytest.raises(ValueError, match='median_C'):
        rep.build_report(out, missingness_fn=lambda m, s: (0.0, 0.0))
