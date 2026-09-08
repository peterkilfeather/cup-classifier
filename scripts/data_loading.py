"""
Shared data loading for cup-classifier analyses.
Reused by diagnostic protocol and Phase 1 modelling pipeline.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
INPUT = BASE / 'input'
META = INPUT / 'metadata' / 'metadata_cleaned.csv'

PROBE_METH_DIR = INPUT / 'methylation' / 'probe_meth'
PROBE_CPG_DIR = INPUT / 'methylation' / 'probe_cpg'
PROBE_ANNOTATIONS = INPUT / 'methylation' / 'probe_annotations_450k.csv'

# Assumed probe window half-width (bp): the 122 bp probe region centered on
# the 450K mapinfo target (protocol addendum, issue #18). Pending the
# client's pipeline scripts (requested 2026-09-07); the tt39 in-window guard
# re-checks membership under any future change of this default.
DEFAULT_HALF_WIDTH = 61

# tt39 probe panel (docs/references/tt39-probes.txt): 39 Illumina 450K probes
# from the tumortype39 set. Provenance unrecorded in this repo (selection
# basis unknown); treated as leakage-adjacent — see CONTEXT.md and
# docs/version-screen-protocol.md (version screen rows 8-13).
_tt39_text = (BASE / 'docs' / 'references' / 'tt39-probes.txt').read_text()
TT39_PROBES = re.findall(r'\bcg\d{8,}\b', _tt39_text)
if len(TT39_PROBES) != 39:
    raise ValueError(
        f"Expected 39 tt39 probes in docs/references/tt39-probes.txt, "
        f"found {len(TT39_PROBES)}")

MODALITY_CONFIGS = {
    'probe_meth': {
        'file': PROBE_METH_DIR
                / 'all_samples.probe_meth_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'label': 'Methylation (probe-avg)',
    },
    'fem4': {
        'file': INPUT / 'fragmentomic' / 'ALL_fem4_features.tsv',
        'sep': '\t', 'sample_col': 'sample',
        'drop_cols': ['tissue', 'N_motifs'],
        'high_dim': False, 'label': 'FEM4 (256)',
    },
    'fragment_length': {
        'file': INPUT / 'fragmentomic' / 'fragment_length_features_qc.csv',
        'sep': ',', 'sample_col': 'Sample_ID', 'drop_cols': ['Tissue'],
        'high_dim': False, 'label': 'Fragment length (369)',
    },
    'cnvkit': {
        'file': INPUT / 'cnvkit' / 'cnvkit_features_thr0.10.tsv',
        'sep': '\t', 'sample_col': 'sample',
        'drop_cols': ['cns_path', 'tissue', 'cnr_path'],
        'high_dim': False, 'label': 'CNVkit thr0.10',
    },
    'probe_cpg': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'label': 'Per-CpG methylation (32K->PCs)',
    },
    'end_density': {
        'file': INPUT / 'fragmentomic' / 'end_density'
                / 'all_samples_ends_100kb_CPM_matrix.tsv',
        'sep': '\t', 'sample_col': None, 'drop_cols': [],
        'high_dim': True, 'label': 'End density (31K bins->PCs)',
    },
    # ── Version screen rows (issue #16, protocol run matrix) ──
    'probe_meth_unenriched': {
        'file': PROBE_METH_DIR / 'all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'long_format': True,
        'label': 'Methylation (probe-avg, unenriched)',
    },
    'probe_meth_unfiltered_qc': {
        'file': PROBE_METH_DIR / 'all_samples.probe_meth_unfiltered.long_v2.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'long_format': True,
        'label': 'Methylation (probe-avg, unfiltered QC)',
    },
    'probe_cpg_enriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'label': 'Per-CpG methylation (enriched)',
    },
    'probe_cpg_agg_unenriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'aggregate': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (aggregated to probes, unenriched)',
    },
    'probe_cpg_agg_enriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'aggregate': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (aggregated to probes, enriched)',
    },
    'probe_meth_tt39_enriched': {
        'file': PROBE_METH_DIR
                / 'all_samples.probe_meth_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'probe_subset': TT39_PROBES,
        'label': 'Methylation (probe-avg, tt39, enriched)',
    },
    'probe_meth_tt39_unenriched': {
        'file': PROBE_METH_DIR / 'all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'long_format': True, 'probe_subset': TT39_PROBES,
        'label': 'Methylation (probe-avg, tt39, unenriched)',
    },
    'probe_cpg_tt39_unenriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv',
        'probe_subset': TT39_PROBES,
        'label': 'Per-CpG methylation (tt39, unenriched)',
    },
    'probe_cpg_tt39_unenriched_lasso': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'dr': 'lasso',
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv',
        'probe_subset': TT39_PROBES,
        'label': 'Per-CpG methylation (tt39, unenriched, raw LASSO)',
    },
    'probe_cpg_tt39_enriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv',
        'probe_subset': TT39_PROBES,
        'label': 'Per-CpG methylation (tt39, enriched)',
    },
    'probe_cpg_tt39_enriched_lasso': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'dr': 'lasso',
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv',
        'probe_subset': TT39_PROBES,
        'label': 'Per-CpG methylation (tt39, enriched, raw LASSO)',
    },
    # ── Version screen addendum rows (issue #18, protocol rows 14-17) ──
    # In-probe-only rows: feature columns restricted to sites with
    # |pos - mapinfo| <= half_width (default 61), derived at runtime from
    # the manifest x probe_annotations_450k.csv join — nothing hardcoded.
    'probe_cpg_inprobe_unenriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'dr': 'lasso', 'in_probe_only': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (in-probe, unenriched, raw LASSO)',
    },
    'probe_cpg_inprobe_enriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': True, 'dr': 'lasso', 'in_probe_only': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (in-probe, enriched, raw LASSO)',
    },
    'probe_cpg_agg_inprobe_unenriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'in_probe_only': True, 'aggregate': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (aggregated to probes, in-probe, unenriched)',
    },
    'probe_cpg_agg_inprobe_enriched': {
        'file': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
        'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
        'high_dim': False, 'in_probe_only': True, 'aggregate': True,
        'manifest': PROBE_CPG_DIR
                / 'all_samples.probe_cpg_enriched_filtered.mNonCpGlt4.manifest.tsv',
        'label': 'Per-CpG methylation (aggregated to probes, in-probe, enriched)',
    },
}


def load_metadata():
    """Load cleaned metadata as DataFrame with standard normalizations."""
    meta = pd.read_csv(META)
    meta['Tissue'] = meta['Tissue'].str.lower()
    meta['Sex'] = meta['Sex'].str.strip()
    return meta


def get_modality_cfg(mod_name):
    """Return config dict for a modality name."""
    return MODALITY_CONFIGS[mod_name].copy()


def _long_to_wide(df, sample_col):
    """Pivot LONG (sample, probe_id, CpG_frac) rows to a wide sample x probe
    frame. Raises ValueError on duplicate sample x probe rows."""
    dup_mask = df.duplicated(subset=[sample_col, 'probe_id'])
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        ex = df.loc[dup_mask, [sample_col, 'probe_id']].iloc[0]
        raise ValueError(
            f"Duplicate sample x probe rows in LONG feature file: {n_dup} "
            f"duplicate(s); first: {ex[sample_col]} x {ex['probe_id']}.")
    wide = df.pivot(index=sample_col, columns='probe_id', values='CpG_frac')
    return wide.reset_index()


def _aggregate_sites_to_probes(df, sample_col, manifest_path):
    """Mean of observed site betas per probe (all measured sites per probe).

    Maps each per-CpG feature column to its probe via the manifest
    (feature_id -> probe_id). Unweighted mean over observed (non-NaN) sites;
    a probe with no observed sites stays NaN (imputed downstream).
    """
    man = pd.read_csv(manifest_path, sep='\t')
    feat_to_probe = dict(zip(man['feature_id'], man['probe_id']))
    missing = [c for c in df.columns if c != sample_col and c not in feat_to_probe]
    if missing:
        raise ValueError(f"Site columns missing from manifest: {missing[:5]}...")
    probe_cols = {}
    for col in df.columns:
        if col == sample_col:
            continue
        probe_cols.setdefault(feat_to_probe[col], []).append(col)
    out = {sample_col: df[sample_col].values}
    for probe_id, cols in probe_cols.items():
        out[probe_id] = df[cols].mean(axis=1, skipna=True).values
    return pd.DataFrame(out)


def _in_probe_feature_ids(man, half_width):
    """feature_ids whose site lies within half_width of its probe's 450K
    mapinfo (the in-probe site set).

    Joins the manifest chr/pos with probe_annotations_450k.csv chr/mapinfo,
    normalizing chr naming (manifest uses 'chr1', annotations use '1').
    Membership is derived at runtime — nothing probe- or count-specific is
    hardcoded. Raises ValueError when a manifest probe has no annotation row.
    """
    ann = pd.read_csv(PROBE_ANNOTATIONS)
    ann['chr'] = ann['chr'].astype(str)
    ann['mapinfo'] = ann['mapinfo'].astype(float)
    man = man.copy()
    man['chr'] = man['chr'].str.replace(r'^chr', '', regex=True)
    merged = man.merge(ann[['probe_id', 'chr', 'mapinfo']],
                       on=['probe_id', 'chr'], how='left')
    missing = merged.loc[merged['mapinfo'].isna(), 'probe_id'].unique()
    if len(missing):
        raise ValueError(
            f"Probes in manifest missing from {PROBE_ANNOTATIONS.name}: "
            f"{missing.tolist()}")
    in_probe = (merged['pos'] - merged['mapinfo']).abs() <= half_width
    return merged.loc[in_probe, 'feature_id'].tolist()


def load_modality(cfg, meta_samples, impute=True):
    """Load a modality's feature matrix, aligned to a set of sample IDs.

    Parameters
    ----------
    cfg : dict
        Modality config from get_modality_cfg().
    meta_samples : array-like
        Sample ID values to align to (typically metadata TWIST_ID).
    impute : bool
        If True, mean-impute NaN values for well-conditioned (non-high-dim)
        modalities before returning. Set to False to defer imputation to the
        caller (e.g. for per-fold imputation inside CV loops).

    Returns
    -------
    X : ndarray (n_samples, n_features)
        Feature matrix, filtered to samples present in meta_samples.
    sample_ids : ndarray
        Sample IDs for each row of X.
    feat_names : ndarray or None
        Feature names for each column.
    """
    fpath = cfg['file']

    if cfg['sample_col'] is None:
        # End-density format: columns are sample IDs, rows are bins
        df = pd.read_csv(fpath, sep=cfg['sep'])
        sample_cols = [c for c in df.columns if c in meta_samples]
        non_sample = ['chr', 'start', 'end']
        feat_names = df[non_sample].astype(str).agg('_'.join, axis=1).values
        X = df[sample_cols].values.T
        sample_ids = np.array(sample_cols)
    else:
        # Per-CpG rows restricted to subset probes and/or in-probe sites:
        # read only their site columns (feature_id per manifest), skipping
        # the full 32K/54K read
        usecols = None
        if cfg.get('manifest') is not None:
            man = pd.read_csv(cfg['manifest'], sep='\t')
            keep_feats = None
            if cfg.get('in_probe_only'):
                keep_feats = _in_probe_feature_ids(
                    man, cfg.get('half_width', DEFAULT_HALF_WIDTH))
            if cfg.get('probe_subset') is not None:
                missing = [p for p in cfg['probe_subset'] if p not in set(man['probe_id'])]
                if missing:
                    raise ValueError(f"Subset probes missing from manifest: {missing}")
                subset_feats = man.loc[man['probe_id'].isin(cfg['probe_subset']),
                                       'feature_id'].tolist()
                if keep_feats is None:
                    keep_feats = subset_feats
                else:
                    keep_feats = [f for f in subset_feats if f in set(keep_feats)]
            if keep_feats is not None:
                usecols = [cfg['sample_col']] + keep_feats

        df = pd.read_csv(fpath, sep=cfg['sep'], usecols=usecols)
        for col in cfg.get('drop_cols', []):
            if col in df.columns:
                df.drop(columns=[col], inplace=True)

        if cfg.get('long_format'):
            df = _long_to_wide(df, cfg['sample_col'])

        if cfg.get('probe_subset') is not None and cfg.get('manifest') is None:
            # Probe-level (probe_meth) rows: restrict feature columns to subset
            missing = [p for p in cfg['probe_subset'] if p not in df.columns]
            if missing:
                raise ValueError(f"Subset probes missing from feature file: {missing}")
            df = df[[cfg['sample_col']] + list(cfg['probe_subset'])]

        if cfg.get('aggregate'):
            # Aggregation rows (mean of observed site betas per probe);
            # the site set is whatever columns were kept above (all
            # measured sites, or the in-probe set for in_probe_only rows)
            df = _aggregate_sites_to_probes(df, cfg['sample_col'], cfg['manifest'])

        if cfg['sample_col'] in df.columns:
            df.set_index(cfg['sample_col'], inplace=True)
        # Deduplicate: keep first occurrence per sample ID
        dup_mask = df.index.duplicated(keep='first')
        if dup_mask.any():
            df = df[~dup_mask]
        X = df.values.astype(np.float64)
        sample_ids = df.index.values
        feat_names = df.columns.values

    # Inner join with requested sample set
    mask = np.isin(sample_ids, meta_samples)
    X, sample_ids = X[mask], sample_ids[mask]

    # Drop all-NaN features
    nan_feats = np.all(np.isnan(X), axis=0)
    if nan_feats.any():
        X = X[:, ~nan_feats]
        if feat_names is not None:
            feat_names = feat_names[~nan_feats]

    # Mean-impute remaining NaNs (for well-conditioned modalities)
    if impute and not cfg.get('high_dim', False):
        col_mean = np.nanmean(X, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(col_mean, inds[1])

    return X, sample_ids, feat_names
