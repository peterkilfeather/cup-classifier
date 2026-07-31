"""
Shared data loading for cup-classifier analyses.
Reused by diagnostic protocol and Phase 1 modelling pipeline.
"""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
INPUT = BASE / 'input'
META = INPUT / 'metadata' / 'metadata_cleaned.csv'

MODALITY_CONFIGS = {
    'probe_meth': {
        'file': INPUT / 'methylation' / 'probe_meth'
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
        'file': INPUT / 'methylation' / 'probe_cpg'
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


def load_modality(cfg, meta_samples, impute=True):
    """Load a modality's feature matrix, aligned to a set of sample IDs.

    Parameters
    ----------
    cfg : dict
        Modality config from get_modality_cfg().
    meta_samples : array-like
        Sample ID values to align to (typically metadata TWIST_ID).

    Returns
    -------
    X : ndarray (n_samples, n_features)
        Feature matrix, filtered to samples present in meta_samples.
    sample_ids : ndarray
        Sample IDs for each row of X.
    feat_names : ndarray or None
        Feature names for each column.
    impute : bool
        If True, mean-impute NaN values for well-conditioned (non-high-dim)
        modalities before returning. Set to False to defer imputation to the
        caller (e.g. for per-fold imputation inside CV loops).
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
        df = pd.read_csv(fpath, sep=cfg['sep'])
        for col in cfg.get('drop_cols', []):
            if col in df.columns:
                df.drop(columns=[col], inplace=True)
        if cfg['sample_col'] in df.columns:
            df.set_index(cfg['sample_col'], inplace=True)
        # Deduplicate: keep first occurrence per sample ID
        dup_mask = df.index.duplicated(keep='first')
        if dup_mask.any():
            n_dup = dup_mask.sum()
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
