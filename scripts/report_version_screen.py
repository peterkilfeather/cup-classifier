"""Version screen hand-back report (ticket #12).

Reads the 39-run outputs written by scripts/run_phase1_pipeline.py into
output/phase1/ and produces the hand-back table:

- Per-scope table: row, macro-F1 mean ± std (with 95% CI), balanced
  accuracy, n samples, n features, median C, per-source accuracy, and
  missingness (fraction of missing cells + fraction of samples with any
  missing — recomputed from the raw feature files, since the pipeline
  does not emit missingness).
- Per-fold paired deltas vs the same-capture probe-averaged row (TOO
  scopes): per-fold macro-F1 delta, mean delta, fold-consistency.
  Anchors follow docs/version-screen-protocol.md ("pair same captures").
  Also pairs each per-CpG row with its per-CpG aggregation row (same-data
  granularity test — restricted to same-data pairs: no aggregated row
  exists for the tt39 probe subset).

No pipeline changes. Run: python3 scripts/report_version_screen.py
Outputs: output/phase1/version_screen_report.csv, version_screen_paired_deltas.csv
"""
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'scripts'))

import data_loading
import run_phase1_pipeline as pipeline

OUT = pipeline.OUT
log = pipeline.log

# Scope names as accepted by pipeline.get_scope, with summary labels
SCOPES = ['full', 'too', 'too-edta']
SCOPE_TO_SUMMARY = {'full': 'Full', 'too': 'TOO', 'too-edta': 'TOO EDTA'}

# Same-capture probe-averaged anchor per per-CpG row (protocol run matrix).
# probe-averaged rows are anchors themselves and are excluded from the table.
SAME_CAPTURE_ANCHOR = {
    'probe_cpg': 'probe_meth_unenriched',
    'probe_cpg_enriched': 'probe_meth',
    'probe_cpg_agg_unenriched': 'probe_meth_unenriched',
    'probe_cpg_agg_enriched': 'probe_meth',
    'probe_cpg_tt39_unenriched': 'probe_meth_tt39_unenriched',
    'probe_cpg_tt39_unenriched_lasso': 'probe_meth_tt39_unenriched',
    'probe_cpg_tt39_enriched': 'probe_meth_tt39_enriched',
    'probe_cpg_tt39_enriched_lasso': 'probe_meth_tt39_enriched',
}

# Clean granularity test (same data, protocol): per-CpG row vs its
# aggregated-to-probes row. Same-data only — the tt39 per-CpG rows (39-probe
# subsets) have no aggregated counterpart in the run matrix, so they are
# paired against their probe-averaged tt39 anchor (above) and not here.
AGGREGATED_ANCHOR = {
    'probe_cpg': 'probe_cpg_agg_unenriched',
    'probe_cpg_enriched': 'probe_cpg_agg_enriched',
}

# Summary columns the report reads from phase1_summary.csv
_REQUIRED_SUMMARY_COLS = ['Scope', 'Modality', 'macro_F1_mean',
                          'macro_F1_std', 'median_C', 'total_features']


def tag_for(mod_name, scope):
    """Output-file tag for a (modality, scope) pair (mirrors the pipeline:
    'TOO ' / 'TOO EDTA ' label prefixes become 'TOO_' / 'TOO_EDTA_' tags)."""
    prefix = {'full': '', 'too': 'TOO_', 'too-edta': 'TOO_EDTA_'}[scope]
    return f"{prefix}{mod_name}"


def anchor_for(mod_name):
    """Same-capture probe-averaged anchor for a per-CpG row, else None."""
    return SAME_CAPTURE_ANCHOR.get(mod_name)


def paired_fold_deltas(mod_metrics, anchor_metrics):
    """Per-fold macro-F1 deltas of mod vs anchor, joined on fold.

    Returns a DataFrame with one row per common fold plus summary columns
    broadcast to every row: mean_delta and n_folds_positive
    (fold-consistency).
    """
    m = mod_metrics[['fold', 'macro_f1']].rename(
        columns={'macro_f1': 'mod_macro_f1'})
    a = anchor_metrics[['fold', 'macro_f1']].rename(
        columns={'macro_f1': 'anchor_macro_f1'})
    d = m.merge(a, on='fold', how='inner')
    d['delta'] = d['mod_macro_f1'] - d['anchor_macro_f1']
    d['mean_delta'] = d['delta'].mean()
    d['n_folds_positive'] = int((d['delta'] > 0).sum())
    return d


def missingness_stats(X):
    """(fraction of missing cells, fraction of samples with any missing)."""
    n_cells = X.size
    n_missing = int(np.isnan(X).sum())
    n_samples_any = int((np.isnan(X).any(axis=1)).sum())
    return (n_missing / n_cells, n_samples_any / X.shape[0])


def _load_missingness(mod_name, scope, meta):
    """Missingness for one (modality, scope): NaN fraction over the scope
    universe, computed from the raw feature file via the shared loader
    (impute=False — the pipeline's per-fold imputation is not applied)."""
    cfg = data_loading.get_modality_cfg(mod_name)
    scope_meta, _ = pipeline.get_scope(scope, meta)
    X, _, _ = data_loading.load_modality(
        cfg, scope_meta['TWIST_ID'].values, impute=False)
    return missingness_stats(X)


def _per_source_accuracy(metrics_df):
    """Mean per-source accuracy across folds, as {src_acc_<Source>: value}
    (same column naming as the pipeline's cv_metrics CSVs)."""
    return {c: float(metrics_df[c].mean())
            for c in metrics_df.columns if c.startswith('src_acc_')}


def build_report(out_dir=OUT, meta=None, missingness_fn=None):
    """Assemble the hand-back report from run outputs.

    Parameters
    ----------
    out_dir : Path
        Directory holding phase1_summary.csv, {tag}_cv_metrics.csv and
        {tag}_hyperparameters.json.
    meta : pd.DataFrame or None
        Cleaned metadata (load_metadata()); required for missingness when
        missingness_fn is None.
    missingness_fn : callable or None
        (mod_name, scope) -> (cell_frac, sample_frac). Defaults to the
        loader-based computation over the scope universe.

    Returns
    -------
    (summary_df, deltas_df) : (pd.DataFrame, pd.DataFrame)
    """
    if missingness_fn is None:
        if meta is None:
            meta = data_loading.load_metadata()
        missingness_fn = lambda m, s: _load_missingness(m, s, meta)

    summary = pd.read_csv(out_dir / 'phase1_summary.csv')
    missing_cols = [c for c in _REQUIRED_SUMMARY_COLS if c not in summary.columns]
    if missing_cols:
        raise ValueError('phase1_summary.csv missing columns: '
                         + ', '.join(missing_cols))

    rows = []
    delta_rows = []
    for scope in SCOPES:
        scope_label = SCOPE_TO_SUMMARY[scope]
        for mod_name in pipeline.PHASE1_MODALITIES:
            tag = tag_for(mod_name, scope)
            hp_path = out_dir / f'{tag}_hyperparameters.json'
            metrics_path = out_dir / f'{tag}_cv_metrics.csv'
            if not hp_path.exists():
                log(f"[skip] {tag}: no hyperparameters file")
                continue
            if not metrics_path.exists():
                log(f"[skip] {tag}: no cv_metrics file")
                continue
            hp = json.loads(hp_path.read_text())
            metrics_df = pd.read_csv(metrics_path)
            cell_frac, sample_frac = missingness_fn(mod_name, scope)

            hit = summary[(summary['Scope'] == scope_label) &
                          (summary['Modality'] == hp.get('label', ''))]
            if len(hit) == 0:
                # fall back to summary rows whose label matches prefix + config label
                cfg_label = data_loading.MODALITY_CONFIGS[mod_name]['label']
                hit = summary[(summary['Scope'] == scope_label) &
                              (summary['Modality'] == f"{'' if scope == 'full' else scope_label + ' '}{cfg_label}")]
            if len(hit) == 0:
                log(f"[skip] {tag}: no summary row")
                continue
            hit = hit.iloc[0]

            row = {
                'scope': scope_label,
                'row': pipeline.PHASE1_MODALITIES.index(mod_name) + 1,
                'modality': mod_name,
                'label': hit['Modality'],
                'macro_F1_mean': hit['macro_F1_mean'],
                'macro_F1_std': hit['macro_F1_std'],
                'macro_F1_CI95_lower': hit.get('macro_F1_CI95_lower'),
                'macro_F1_CI95_upper': hit.get('macro_F1_CI95_upper'),
                'balanced_accuracy': metrics_df['balanced_accuracy'].mean(),
                'n_samples': hp.get('n_samples'),
                'n_features': hit.get('total_features'),
                'n_original_features': hit.get('n_original_features'),
                'median_C': hit['median_C'],
                'missing_cell_frac': float(cell_frac),
                'missing_sample_frac': float(sample_frac),
            }
            row.update(_per_source_accuracy(metrics_df))
            rows.append(row)

            if scope == 'full':
                continue  # paired deltas are TOO-scope only (protocol)
            for kind, anchor_map in [('probe-averaged', SAME_CAPTURE_ANCHOR),
                                     ('probe-aggregation', AGGREGATED_ANCHOR)]:
                anchor = anchor_map.get(mod_name)
                if anchor is None:
                    continue
                a_metrics = out_dir / f'{tag_for(anchor, scope)}_cv_metrics.csv'
                if not a_metrics.exists():
                    log(f"[skip] {tag}: anchor {anchor} metrics missing")
                    continue
                d = paired_fold_deltas(metrics_df, pd.read_csv(a_metrics))
                d.insert(0, 'scope', scope_label)
                d.insert(1, 'row', pipeline.PHASE1_MODALITIES.index(mod_name) + 1)
                d.insert(2, 'modality', mod_name)
                d.insert(3, 'anchor', anchor)
                d.insert(4, 'anchor_type', kind)
                delta_rows.append(d)

    summary_df = pd.DataFrame(rows)
    deltas_df = (pd.concat(delta_rows, ignore_index=True)
                 if delta_rows else pd.DataFrame())
    return summary_df, deltas_df


def main():
    log("Version screen report (ticket #12)")
    log("=" * 70)
    meta = data_loading.load_metadata()
    summary, deltas = build_report(out_dir=OUT, meta=meta)

    if summary.empty:
        log("No results found. Run the screen first "
            "(run_phase1_pipeline.py --scopes full too too-edta).")
        return 1

    summary_path = OUT / 'version_screen_report.csv'
    summary.to_csv(summary_path, index=False)
    log(f"Saved: {summary_path}")
    print()
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    if not deltas.empty:
        deltas_path = OUT / 'version_screen_paired_deltas.csv'
        deltas.to_csv(deltas_path, index=False)
        log(f"Saved: {deltas_path}")
        # One row per pairing: the broadcast summary columns are identical
        # across a pairing's fold rows
        per_row = (deltas[['scope', 'row', 'modality', 'anchor', 'anchor_type',
                           'mean_delta', 'n_folds_positive']]
                   .drop_duplicates()
                   .reset_index(drop=True))
        print()
        print(per_row.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    else:
        log("No paired deltas (need TOO-scope run outputs).")
    return 0


if __name__ == '__main__':
    sys.exit(main())
