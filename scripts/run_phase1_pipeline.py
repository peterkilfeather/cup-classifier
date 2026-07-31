#!/usr/bin/env python3
"""
Phase 1: Baseline tissue classifier on current 164 samples.

Pipeline:
  For each modality (FEM4 primary, probe_meth secondary, fragment_length tertiary):
    1. Inner join feature matrix with cleaned metadata
    2. Source-stratified CV: StratifiedKFold (by Tissue) + source-coverage check
    3. Per fold: StandardScaler → L1-logreg with inner GridSearchCV for C
    4. Report: macro-F1 (primary), balanced accuracy (secondary), per-source accuracy
    5. EDTA-only sensitivity on 96-sample subset
    6. Refit on all 164 samples → serialized model for Phase 2a

Design per protocol (docs/multi-source-validation-protocol.md §2).
"""

import sys, os, warnings, json, time, argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score
from sklearn.exceptions import ConvergenceWarning
import joblib

warnings.filterwarnings('ignore', category=ConvergenceWarning)
warnings.filterwarnings('ignore', category=UserWarning)
np.random.seed(42)

# ── Paths ─────────────────────────────────────────────

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'scripts'))
from data_loading import load_metadata, get_modality_cfg, load_modality

OUT = BASE / 'output' / 'phase1'
OUT.mkdir(parents=True, exist_ok=True)
FIGS = OUT / 'figures'
FIGS.mkdir(exist_ok=True)
MODELS = OUT / 'models'
MODELS.mkdir(exist_ok=True)

RANDOM_STATE = 42
N_SPLITS = 5
N_PCS = 20  # fixed PCs for high-dim modalities
C_GRID = np.logspace(-3, 1, 6)  # [0.001, 0.01, 0.1, 1, 10]
MAX_ITER = 5000
MAX_RESHUF = 100

# Modalities to run (in priority order)
PHASE1_MODALITIES = ['fem4', 'probe_meth', 'fragment_length', 'probe_cpg', 'end_density']

PALETTE_TISSUE = {
    'colon': '#e41a1c', 'liver': '#377eb8', 'pancreas': '#4daf4a',
    'prostate': '#984ea3', 'stomach': '#ff7f00', 'healthyblood': '#a65628',
}
PALETTE_SOURCE = {
    'Fox Chase': '#1f78b4', 'Audubon': '#33a02c',
    'Sowalsky': '#e31a1c', 'NIH_Clinical_Center': '#ff7f00',
}


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


# ── Source-covering CV ────────────────────────────────

def source_covering_split(meta, n_splits=N_SPLITS, random_state=RANDOM_STATE,
                          max_retries=MAX_RESHUF):
    """StratifiedKFold stratified by Tissue, rejecting splits where any
    training fold lacks all sources present in the dataset.
    Number of required sources is determined from the data (e.g. 4 for full
    dataset, 3 for EDTA subset). Retries with incremented seed.

    Returns
    -------
    splits : list of (train_idx, test_idx) arrays
    n_retries : int
    """
    y = meta['Tissue'].values
    sources = meta['Source'].values
    required_sources = set(meta['Source'].unique())

    for attempt in range(max_retries):
        seed = random_state + attempt
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(skf.split(np.zeros(len(meta)), y))
        # Verify each training fold
        ok = True
        for train_idx, _ in splits:
            fold_sources = set(sources[train_idx])
            if not required_sources.issubset(fold_sources):
                ok = False
                break
        if ok:
            return splits, attempt

    log(f"WARNING: Could not satisfy source coverage after {max_retries} retries. "
        "Using best-effort split.")
    return splits, max_retries


# ── PCA for high-dim modalities ────────────────────────

PCA_N_COMPONENTS = N_PCS


def pca_fit_transform(X_tr, X_te=None):
    """Impute NaN, scale, PCA-fit on training, transform training (and test).

    NaN imputation uses training column means (no leakage). PCA is fit on
    training data only. Number of components = min(N_PCS, n_samples, n_features).

    Returns (dual arity based on X_te)
    -----------------------------------
    If X_te is None:
        X_tr_pc : ndarray (n_tr_samples, k)  — PCA-reduced training data
        pca : PCA object                      — fitted PCA
        scaler : StandardScaler               — fitted scaler
        k : int                               — number of components

    If X_te is not None:
        X_tr_pc : ndarray (n_tr_samples, k)  — PCA-reduced training data
        X_te_pc : ndarray (n_te_samples, k)  — PCA-reduced test data
        pca : PCA object                      — fitted PCA
        scaler : StandardScaler               — fitted scaler
        k : int                               — number of components
    """
    # Impute training NaN using training column means
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        col_mean = np.nanmean(X_tr, axis=0)
    col_mean = np.nan_to_num(col_mean, nan=0.0)
    X_tr_imp = X_tr.copy()
    inds_tr = np.where(np.isnan(X_tr_imp))
    X_tr_imp[inds_tr] = np.take(col_mean, inds_tr[1])

    # Scale
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X_tr_imp)

    # Fixed K PCs (L1 will zero uninformative components)
    k = min(PCA_N_COMPONENTS, X_tr_s.shape[0], X_tr_s.shape[1])
    pca = PCA(n_components=k)
    X_tr_pc = pca.fit_transform(X_tr_s)
    cumvar = pca.explained_variance_ratio_.sum() * 100
    log(f"  PCA K={k}, cum var={cumvar:.1f}%")

    if X_te is not None:
        X_te_imp = X_te.copy()
        inds_te = np.where(np.isnan(X_te_imp))
        X_te_imp[inds_te] = np.take(col_mean, inds_te[1])
        X_te_s = scaler.transform(X_te_imp)
        X_te_pc = pca.transform(X_te_s)
        return X_tr_pc, X_te_pc, pca, scaler, k

    return X_tr_pc, pca, scaler, k


def _impute_with_mean(X, col_mean):
    """Return a copy of X with NaN values replaced by col_mean."""
    X_out = X.copy()
    inds = np.where(np.isnan(X_out))
    X_out[inds] = np.take(col_mean, inds[1])
    return X_out


# ── Inner CV: C tuning ────────────────────────────────

def tune_C(X_tr, y_tr, random_state=RANDOM_STATE):
    """Tune L1 penalty C via inner GridSearchCV, return fitted model.

    Falls back to C=0.1 without grid search when min class count < 2.
    Grid: C_GRID (np.logspace(-3, 1, 6)).
    Inner CV: StratifiedKFold(n_splits=min(3, min_class_count)).

    Parameters
    ----------
    X_tr : ndarray (n_train, n_features)
        Training features, already scaled.
    y_tr : ndarray (n_train,)
        Training labels, already encoded.
    random_state : int
        Seed for inner CV shuffling.

    Returns
    -------
    LogisticRegression
        Fitted best estimator (tuned C + fit on full X_tr).
        Access selected C via `.C`.
    """
    n_inner = min(3, np.min(np.bincount(y_tr)))
    if n_inner < 2:
        clf = LogisticRegression(
            l1_ratio=1, solver='saga',
            C=0.1, class_weight='balanced', max_iter=MAX_ITER,
            random_state=random_state)
        clf.fit(X_tr, y_tr)
        return clf

    inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=random_state)
    clf = LogisticRegression(
        l1_ratio=1, solver='saga',
        class_weight='balanced', max_iter=MAX_ITER,
        random_state=random_state)
    gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
    gs.fit(X_tr, y_tr)
    return gs.best_estimator_


# ── Per-modality pipeline ─────────────────────────────

def run_modality_pipeline(mod_name, meta, label_prefix=''):
    """Run Phase 1 pipeline for one modality.

    Parameters
    ----------
    mod_name : str
        Key into MODALITY_CONFIGS.
    meta : pd.DataFrame
        Cleaned metadata (filtered for scope, e.g. EDTA-only).
    label_prefix : str
        Prefix for output labels (e.g. 'EDTA ').

    Returns
    -------
    dict with results, or None on failure.
    """
    t0 = time.time()
    cfg = get_modality_cfg(mod_name)
    modality_label = cfg['label']
    full_label = f"{label_prefix}{modality_label}"
    log(f"{'='*60}")
    log(f"Starting {full_label} ({len(meta)} samples)")

    # Load features
    try:
        X, sample_ids, feat_names = load_modality(cfg, meta['TWIST_ID'].values, impute=False)
    except Exception as e:
        log(f"FAILED load: {e}")
        return None

    if X.shape[0] < 5 or X.shape[1] < 2:
        log(f"Too few samples/features: {X.shape}")
        return None

    # Align metadata to loaded samples
    meta_aligned = meta.set_index('TWIST_ID').loc[sample_ids].reset_index()
    log(f"X shape: {X.shape}")

    # ── Source-stratified CV ──────────────────────
    splits, n_retries = source_covering_split(meta_aligned)
    log(f"CV: {N_SPLITS}-fold (re-shuffled {n_retries}x)")

    le = LabelEncoder()
    y = le.fit_transform(meta_aligned['Tissue'].values)
    source_labels = meta_aligned['Source'].values
    is_high_dim = cfg.get('high_dim', False)

    fold_metrics = []
    fold_models = []
    fold_Cs = []

    for fold, (train_idx, test_idx) in enumerate(splits):
        fold_start = time.time()
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        src_tr, src_te = source_labels[train_idx], source_labels[test_idx]

        pca = None
        pca_scaler = None
        fold_k = None
        if is_high_dim:
            # PCA reduction: impute → scale → PCA (fit on train, transform both)
            X_tr_pc, X_te_pc, pca, pca_scaler, fold_k = pca_fit_transform(X_tr, X_te)
            fold_X_tr, fold_X_te = X_tr_pc, X_te_pc
        else:
            # Per-fold NaN imputation (fit on train, transform test)
            col_mean = np.nanmean(X_tr, axis=0)
            col_mean = np.nan_to_num(col_mean, nan=0.0)
            X_tr = _impute_with_mean(X_tr, col_mean)
            X_te = _impute_with_mean(X_te, col_mean)
            # Standard scaling
            scaler = StandardScaler()
            fold_X_tr = scaler.fit_transform(X_tr)
            fold_X_te = scaler.transform(X_te)

        # Tune C on training fold
        model = tune_C(fold_X_tr, y_tr)
        fold_Cs.append(model.C)

        # Predict
        y_pred = model.predict(fold_X_te)

        # Metrics
        macro_f1 = f1_score(y_te, y_pred, average='macro')
        bal_acc = balanced_accuracy_score(y_te, y_pred)
        acc = accuracy_score(y_te, y_pred)

        # Per-source accuracy
        source_acc = {}
        for src in np.unique(source_labels):
            src_mask = src_te == src
            if src_mask.sum() > 0:
                source_acc[src] = accuracy_score(y_te[src_mask], y_pred[src_mask])

        n_nonzero = np.any(model.coef_ != 0, axis=0).sum()
        log(f"  Fold {fold+1}/{N_SPLITS}: C={model.C:.4f}, "
            f"macro-F1={macro_f1:.4f}, bal-acc={bal_acc:.4f}, "
            f"non-zero coeffs={n_nonzero}, [{time.time()-fold_start:.0f}s]")

        fold_metrics.append({
            'fold': fold + 1,
            'C': model.C,
            'macro_f1': macro_f1,
            'balanced_accuracy': bal_acc,
            'accuracy': acc,
            'n_nonzero_coef': n_nonzero,
            'train_samples': len(train_idx),
            'test_samples': len(test_idx),
            **{f'src_acc_{s}': source_acc.get(s, None) for s in np.unique(source_labels)},
        })

        fold_models.append({
            'scaler': pca_scaler if is_high_dim else scaler,
            'pca': pca if is_high_dim else None,
            'model': model,
            'C': model.C,
            'n_pcs': fold_k,
            'train_idx': train_idx,
            'test_idx': test_idx,
            'selected_features': feat_names,
            'feature_indices': np.where(np.any(model.coef_ != 0, axis=0))[0],
        })

    # ── Aggregate metrics ─────────────────────────
    metrics_df = pd.DataFrame(fold_metrics)
    macro_f1_mean = metrics_df['macro_f1'].mean()
    macro_f1_std = metrics_df['macro_f1'].std()

    log(f"\n  {full_label} CV results:")
    log(f"  macro-F1: {macro_f1_mean:.4f} ± {macro_f1_std:.4f} "
        f"(95% CI: [{macro_f1_mean-1.96*macro_f1_std:.4f}, "
        f"{macro_f1_mean+1.96*macro_f1_std:.4f}])")
    log(f"  balanced accuracy: {metrics_df['balanced_accuracy'].mean():.4f}")
    log(f"  chance (uniform): {1/len(le.classes_):.4f}")

    # Per-source accuracy across folds
    source_cols = [c for c in metrics_df.columns if c.startswith('src_acc_')]
    for col in source_cols:
        vals = metrics_df[col].dropna()
        if len(vals) > 0:
            log(f"  {col.replace('src_acc_','')}: {vals.mean():.4f} ± {vals.std():.4f}")

    # ── Refit on all data for Phase 2a ─────────────
    log(f"  Refitting on all {len(meta_aligned)} samples for Phase 2a...")
    median_C = np.median(fold_Cs)

    if is_high_dim:
        # Impute, scale, PCA on all data
        col_mean_all = np.nanmean(X, axis=0)
        col_mean_all = np.nan_to_num(col_mean_all, nan=0.0)
        X_imp = X.copy()
        inds_all = np.where(np.isnan(X_imp))
        X_imp[inds_all] = np.take(col_mean_all, inds_all[1])

        pca_scaler_full = StandardScaler()
        X_s = pca_scaler_full.fit_transform(X_imp)
        k_full = min(N_PCS, X_s.shape[0], X_s.shape[1])
        pca_full = PCA(n_components=k_full)
        X_full_pc = pca_full.fit_transform(X_s)
        cumvar_full = pca_full.explained_variance_ratio_.sum() * 100
        log(f"  Full PCA: K={k_full}, cum var={cumvar_full:.1f}%")

        X_full_for_clf = X_full_pc
        pc_names = np.array([f'PC{i+1}' for i in range(k_full)])
        full_feat_names = pc_names
    else:
        # Impute NaN on all data (safe for full refit; no test set involved)
        col_mean_all = np.nanmean(X, axis=0)
        col_mean_all = np.nan_to_num(col_mean_all, nan=0.0)
        X_imp = _impute_with_mean(X, col_mean_all)
        scaler_full = StandardScaler()
        X_full_for_clf = scaler_full.fit_transform(X_imp)
        full_feat_names = feat_names

    model_full = LogisticRegression(
        l1_ratio=1, solver='saga',
        C=median_C, class_weight='balanced', max_iter=MAX_ITER,
        random_state=RANDOM_STATE)
    model_full.fit(X_full_for_clf, y)

    n_nonzero_full = int(np.any(model_full.coef_ != 0, axis=0).sum())
    log(f"  Full model: C={median_C:.4f}, non-zero coeffs={n_nonzero_full}/{X_full_for_clf.shape[1]}")

    # Selected features (PC names for high-dim, original names for well-conditioned)
    selected_mask = np.any(model_full.coef_ != 0, axis=0)
    selected_features = full_feat_names[selected_mask] if full_feat_names is not None else None

    # ── Save artifacts ────────────────────────────
    safe_name = mod_name.replace(' ', '_')
    tag = f"{label_prefix}{safe_name}".strip('_').replace(' ', '_')

    model_dict = {
        'model': model_full,
        'C': median_C,
        'classes': le.classes_.tolist(),
        'modality': mod_name,
        'modality_label': full_label,
        'timestamp': datetime.now().isoformat(),
        'n_samples': len(meta_aligned),
        'cv_macro_f1_mean': macro_f1_mean,
        'cv_macro_f1_std': macro_f1_std,
        'is_high_dim': is_high_dim,
    }

    if is_high_dim:
        # PCA path: scaler + PCA + model
        model_dict.update({
            'pca_scaler': pca_scaler_full,
            'pca': pca_full,
            'n_pcs': k_full,
            'n_original_features': X.shape[1],
            'pca_var_ratio': pca_full.explained_variance_ratio_.tolist(),
        })
        # For high-dim, "feature_names" refers to original raw features
        model_dict['feature_names'] = feat_names
        model_dict['pc_names'] = pc_names.tolist()
        model_dict['selected_features'] = selected_features.tolist() if selected_features is not None else None
    else:
        model_dict.update({
            'scaler': scaler_full,
            'feature_names': feat_names,
            'selected_features': selected_features,
            'selected_mask': selected_mask,
        })

    model_path = MODELS / f'{tag}_full_model.joblib'
    joblib.dump(model_dict, model_path)
    log(f"  Saved: {model_path}")

    # Save CV metrics
    csv_path = OUT / f'{tag}_cv_metrics.csv'
    metrics_df.to_csv(csv_path, index=False)
    log(f"  Saved: {csv_path}")

    # Set k_full (only meaningful for high_dim)
    k_full = locals().get('k_full', None)

    hp = {
        'modality': mod_name,
        'label': full_label,
        'timestamp': datetime.now().isoformat(),
        'n_splits': N_SPLITS,
        'max_retries': n_retries,
        'C_grid': C_GRID.tolist(),
        'median_C': median_C,
        'per_fold_C': fold_Cs,
        'macro_f1_mean': macro_f1_mean,
        'macro_f1_std': macro_f1_std,
        'n_samples': len(meta_aligned),
        'n_features': X.shape[1],
        'n_pcs': k_full,
        'n_nonzero_coef': int(n_nonzero_full),
        'random_state': RANDOM_STATE,
        'classifier': 'LogisticRegression(L1, saga, multinomial)',
        'cv_strategy': 'StratifiedKFold(5, stratified_by_Tissue, source_coverage_check)',
        'high_dim': is_high_dim,
    }
    hp_path = OUT / f'{tag}_hyperparameters.json'
    with open(hp_path, 'w') as f:
        json.dump(hp, f, indent=2)
    log(f"  Saved: {hp_path}")

    # ── Plots ─────────────────────────────────────
    _plot_cv_results(metrics_df, full_label, tag, chance=1/len(le.classes_))
    _plot_coefficient_heatmap(model_full, full_feat_names, le.classes_, full_label, tag)

    elapsed = time.time() - t0
    log(f"Done {full_label} in {elapsed:.0f}s")
    log(f"{'='*60}\n")

    return {
        'mod_name': mod_name,
        'label': full_label,
        'metrics_df': metrics_df,
        'macro_f1_mean': macro_f1_mean,
        'macro_f1_std': macro_f1_std,
        'fold_Cs': fold_Cs,
        'median_C': median_C,
        'n_nonzero_full': n_nonzero_full,
        'n_total_features': X.shape[1],
        'n_pcs': k_full,
        'model_full': model_full,
        'selected_features': selected_features,
        'is_high_dim': is_high_dim,
    }


# ── Plots ─────────────────────────────────────────────

def _plot_cv_results(metrics_df, label, tag, chance):
    """Per-fold macro-F1 bar chart (left) + per-source accuracy grouped bars (right).

    Parameters
    ----------
    metrics_df : pd.DataFrame
        With columns: fold, macro_f1, balanced_accuracy, accuracy,
        src_acc_<SourceName> for each source present.
    label : str
        Modality label for plot title (e.g. 'FEM4 (256)').
    tag : str
        File-safe tag for output filename.
    chance : float
        Chance-level baseline (e.g. 1/6 for 6 classes).

    Saves
    -----
    {FIGS}/{tag}_cv_results.png
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Left: macro-F1 per fold
    ax = axes[0]
    folds = metrics_df['fold'].values
    f1_vals = metrics_df['macro_f1'].values
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(folds)))
    ax.bar(folds, f1_vals, color=colors, edgecolor='k', linewidth=0.5, width=0.6)
    mean_f1 = f1_vals.mean()
    ax.axhline(y=mean_f1, color='red', linestyle='--', alpha=0.7,
               label=f'Mean = {mean_f1:.3f}')
    ax.axhline(y=chance, color='gray', linestyle=':', alpha=0.5,
               label=f'Chance ({chance:.3f})')
    ax.set_xlabel('Fold')
    ax.set_ylabel('Macro-F1')
    ax.set_title(f'{label} — Per-Fold Macro-F1')
    ax.set_xticks(folds)
    ax.legend(fontsize=8)

    # Right: per-source accuracy
    ax = axes[1]
    src_cols = [c for c in metrics_df.columns if c.startswith('src_acc_')]
    if src_cols:
        src_names = [c.replace('src_acc_', '') for c in src_cols]
        x = np.arange(len(src_names))
        width = 0.25
        for i, fold in enumerate(folds):
            vals = [metrics_df.loc[metrics_df['fold'] == fold, c].values[0]
                    if pd.notna(metrics_df.loc[metrics_df['fold'] == fold, c].values[0])
                    else 0 for c in src_cols]
            offset = (i - len(folds)/2 + 0.5) * width
            ax.bar(x + offset, vals, width, label=f'Fold {fold}',
                   color=plt.cm.viridis(i / len(folds)), edgecolor='k', linewidth=0.3)
        ax.set_xlabel('Source')
        ax.set_ylabel('Accuracy')
        ax.set_title(f'{label} — Per-Source Accuracy')
        ax.set_xticks(x)
        ax.set_xticklabels(src_names, fontsize=8)
        ax.legend(fontsize=7, loc='lower left')
        ax.set_ylim(0, 1.05)

    fig.tight_layout()
    fig.savefig(FIGS / f'{tag}_cv_results.png', dpi=150)
    plt.close(fig)


def _plot_coefficient_heatmap(model, feat_names, classes, label, tag):
    """Heatmap of non-zero L1 coefficients per class (top 30 features).

    Parameters
    ----------
    model : LogisticRegression
        Fitted model with .coef_ attribute (n_classes, n_features).
    feat_names : ndarray or None
        Feature names for x-axis labels. Falls back to indices if None.
    classes : ndarray
        Class labels for y-axis.
    label : str
        Modality label for plot title.
    tag : str
        File-safe tag for output filename.

    Saves
    -----
    {FIGS}/{tag}_coefficients.png
    Skips if all coefficients are zero.
    """
    coef = model.coef_  # (n_classes, n_features)
    non_zero = np.any(coef != 0, axis=0)
    if non_zero.sum() == 0:
        log(f"  All coefficients zero — skipping coefficient heatmap")
        return
    # Limit to top features per class for readability
    top_k = min(30, non_zero.sum())
    coef_subset = coef[:, non_zero]
    if feat_names is not None:
        feat_subset = feat_names[non_zero]
    else:
        feat_subset = np.arange(coef_subset.shape[1]).astype(str)

    # Select top features by max absolute coefficient across classes
    max_abs = np.max(np.abs(coef_subset), axis=0)
    top_idx = np.argsort(max_abs)[-top_k:]

    fig, ax = plt.subplots(figsize=(max(8, top_k * 0.4), max(4, len(classes) * 0.5)))
    vmax = max(np.abs(coef_subset[:, top_idx]).max(), 1e-10)
    sns.heatmap(coef_subset[:, top_idx], ax=ax,
                xticklabels=feat_subset[top_idx], yticklabels=classes,
                cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
                cbar_kws={'label': 'Coefficient'})
    ax.set_xlabel('Feature')
    ax.set_title(f'{label} — L1 Coefficients (top {top_k} features)')
    plt.setp(ax.get_xticklabels(), fontsize=6, rotation=45, ha='right')
    fig.tight_layout()
    fig.savefig(FIGS / f'{tag}_coefficients.png', dpi=150)
    plt.close(fig)


# ── Combined modality pipeline ─────────────────────────

COMBINE_MODALITIES = ['fem4', 'probe_meth', 'fragment_length']

def run_combined_pipeline(mod_names, meta, label_prefix=''):
    """Run Phase 1 pipeline with concatenated features from multiple modalities.

    Each modality is loaded raw, then scaled INDEPENDENTLY per fold
    (no across-fold leakage). Concatenated after scaling.
    """
    t0 = time.time()
    full_label = f"{label_prefix}Combined ({'+'.join(mod_names)})"
    log(f"{'='*60}")
    log(f"Starting {full_label} ({len(meta)} samples)")

    # Load raw features for each modality, aligned to metadata order
    X_raw = {}
    all_feat_names = {}
    n_features_list = []
    first_ids = None

    for mod_name in mod_names:
        cfg = get_modality_cfg(mod_name)
        try:
            X_mod, ids_mod, feats_mod = load_modality(cfg, meta['TWIST_ID'].values, impute=False)
        except Exception as e:
            log(f"  Failed to load {mod_name}: {e}")
            return None

        if first_ids is None:
            first_ids = ids_mod

        # Align rows to first modality's order (FEM4 as reference)
        if not np.array_equal(ids_mod, first_ids):
            order_map = {s: i for i, s in enumerate(first_ids)}
            idx = np.array([order_map[s] for s in ids_mod])
            X_mod = X_mod[idx]
            ids_mod = first_ids

        prefix = mod_name.replace('_meth', '').replace('_length', '')
        prefixed = np.array([f'{prefix}_{f}' for f in (feats_mod if feats_mod is not None
                                                         else [str(i) for i in range(X_mod.shape[1])])])
        X_raw[mod_name] = X_mod
        all_feat_names[mod_name] = prefixed
        n_features_list.append(X_mod.shape[1])

    combined_feat_names = np.concatenate(list(all_feat_names.values()))
    total_n = sum(n_features_list)
    log(f"Raw features loaded: {dict(zip(mod_names, n_features_list))}, total={total_n}")

    # Align metadata
    meta_aligned = meta.set_index('TWIST_ID').loc[first_ids].reset_index()

    # ── Source-stratified CV (per-fold scaling) ───
    splits, n_retries = source_covering_split(meta_aligned)
    log(f"CV: {N_SPLITS}-fold (re-shuffled {n_retries}x)")

    le = LabelEncoder()
    y = le.fit_transform(meta_aligned['Tissue'].values)
    source_labels = meta_aligned['Source'].values

    fold_metrics = []
    fold_Cs = []

    for fold, (train_idx, test_idx) in enumerate(splits):
        fold_start = time.time()

        # Scale each modality per-fold, then concatenate
        fold_Xs_tr = []
        fold_Xs_te = []
        for mod_name in mod_names:
            Xm = X_raw[mod_name]
            Xm_tr = Xm[train_idx]
            Xm_te = Xm[test_idx]
            # Per-fold NaN imputation (fit on train, transform test)
            col_mean = np.nanmean(Xm_tr, axis=0)
            col_mean = np.nan_to_num(col_mean, nan=0.0)
            Xm_tr = _impute_with_mean(Xm_tr, col_mean)
            Xm_te = _impute_with_mean(Xm_te, col_mean)
            scaler = StandardScaler()
            Xm_tr_s = scaler.fit_transform(Xm_tr)
            Xm_te_s = scaler.transform(Xm_te)
            fold_Xs_tr.append(Xm_tr_s)
            fold_Xs_te.append(Xm_te_s)

        X_tr = np.hstack(fold_Xs_tr)
        X_te = np.hstack(fold_Xs_te)
        y_tr, y_te = y[train_idx], y[test_idx]

        model = tune_C(X_tr, y_tr)
        fold_Cs.append(model.C)
        y_pred = model.predict(X_te)

        macro_f1 = f1_score(y_te, y_pred, average='macro')
        bal_acc = balanced_accuracy_score(y_te, y_pred)
        acc = accuracy_score(y_te, y_pred)

        source_acc = {}
        for src in np.unique(source_labels):
            src_mask = (source_labels[test_idx] == src)
            if src_mask.sum() > 0:
                source_acc[src] = accuracy_score(y_te[src_mask], y_pred[src_mask])

        n_nonzero = np.any(model.coef_ != 0, axis=0).sum()
        log(f"  Fold {fold+1}/{N_SPLITS}: C={model.C:.4f}, macro-F1={macro_f1:.4f}, "
            f"non-zero={n_nonzero}/{total_n}, [{time.time()-fold_start:.0f}s]")

        fold_metrics.append({
            'fold': fold + 1,
            'C': model.C,
            'macro_f1': macro_f1,
            'balanced_accuracy': bal_acc,
            'accuracy': acc,
            'n_nonzero_coef': n_nonzero,
            'train_samples': len(train_idx),
            'test_samples': len(test_idx),
            **{f'src_acc_{s}': source_acc.get(s, None) for s in np.unique(source_labels)},
        })

    # ── Aggregate metrics ─────────────────────────
    metrics_df = pd.DataFrame(fold_metrics)
    macro_f1_mean = metrics_df['macro_f1'].mean()
    macro_f1_std = metrics_df['macro_f1'].std()

    log(f"\n  {full_label} CV results:")
    log(f"  macro-F1: {macro_f1_mean:.4f} ± {macro_f1_std:.4f}")
    log(f"  chance (uniform): {1/len(le.classes_):.4f}")

    source_cols = [c for c in metrics_df.columns if c.startswith('src_acc_')]
    for col in source_cols:
        vals = metrics_df[col].dropna()
        if len(vals) > 0:
            log(f"  {col.replace('src_acc_','')}: {vals.mean():.4f}")

    # ── Refit on all data for Phase 2a ─────────────
    log(f"  Refitting on all {len(meta_aligned)} samples for Phase 2a...")
    median_C = np.median(fold_Cs)

    # Impute NaN, then scale each modality on all data, concatenate
    full_scalers = {}
    X_full_pieces = []
    for mod_name in mod_names:
        Xm = X_raw[mod_name]
        col_mean = np.nanmean(Xm, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        Xm = _impute_with_mean(Xm, col_mean)
        scaler = StandardScaler()
        X_full_pieces.append(scaler.fit_transform(Xm))
        full_scalers[mod_name] = scaler
    X_full = np.hstack(X_full_pieces)

    model_full = LogisticRegression(
        l1_ratio=1, solver='saga',
        C=median_C, class_weight='balanced', max_iter=MAX_ITER,
        random_state=RANDOM_STATE)
    model_full.fit(X_full, y)

    n_nonzero_full = int(np.any(model_full.coef_ != 0, axis=0).sum())
    log(f"  Full model: C={median_C:.4f}, non-zero={n_nonzero_full}/{total_n}")

    selected_mask = np.any(model_full.coef_ != 0, axis=0)
    selected_features = combined_feat_names[selected_mask]

    # ── Save artifacts ────────────────────────────
    safe_mods = '_'.join(mod_names)
    tag = f"{label_prefix}combined_{safe_mods}".strip('_')

    model_path = MODELS / f'{tag}_full_model.joblib'
    joblib.dump({
        'model': model_full,
        'scaler_per_modality': full_scalers,
        'modalities': mod_names,
        'C': median_C,
        'classes': le.classes_.tolist(),
        'feature_names': combined_feat_names,
        'selected_features': selected_features,
        'selected_mask': selected_mask,
        'modality': 'combined',
        'modality_label': full_label,
        'timestamp': datetime.now().isoformat(),
        'n_samples': len(meta_aligned),
        'cv_macro_f1_mean': macro_f1_mean,
        'cv_macro_f1_std': macro_f1_std,
        'n_features_per_modality': dict(zip(mod_names, n_features_list)),
    }, model_path)
    log(f"  Saved: {model_path}")

    csv_path = OUT / f'{tag}_cv_metrics.csv'
    metrics_df.to_csv(csv_path, index=False)
    log(f"  Saved: {csv_path}")

    hp = {
        'modality': 'combined',
        'label': full_label,
        'timestamp': datetime.now().isoformat(),
        'modalities': mod_names,
        'n_splits': N_SPLITS,
        'max_retries': n_retries,
        'C_grid': C_GRID.tolist(),
        'median_C': median_C,
        'per_fold_C': fold_Cs,
        'macro_f1_mean': macro_f1_mean,
        'macro_f1_std': macro_f1_std,
        'n_samples': len(meta_aligned),
        'n_features': total_n,
        'n_features_per_modality': dict(zip(mod_names, n_features_list)),
        'n_nonzero_coef': int(n_nonzero_full),
        'random_state': RANDOM_STATE,
        'classifier': 'LogisticRegression(L1, saga, multinomial)',
        'cv_strategy': 'StratifiedKFold(5, stratified_by_Tissue, source_coverage_check)',
    }
    hp_path = OUT / f'{tag}_hyperparameters.json'
    with open(hp_path, 'w') as f:
        json.dump(hp, f, indent=2)
    log(f"  Saved: {hp_path}")

    # ── Plots ─────────────────────────────────────
    _plot_cv_results(metrics_df, full_label, tag, chance=1/len(le.classes_))
    _plot_coefficient_heatmap(model_full, combined_feat_names, le.classes_, full_label, tag)

    elapsed = time.time() - t0
    log(f"Done {full_label} in {elapsed:.0f}s")
    log(f"{'='*60}\n")

    return {
        'mod_name': 'combined',
        'label': full_label,
        'metrics_df': metrics_df,
        'macro_f1_mean': macro_f1_mean,
        'macro_f1_std': macro_f1_std,
        'fold_Cs': fold_Cs,
        'median_C': median_C,
        'n_nonzero_full': n_nonzero_full,
        'n_total_features': total_n,
        'n_pcs': None,
        'model_full': model_full,
        'selected_features': selected_features,
        'is_high_dim': False,
    }


# ── Summary table ─────────────────────────────────────

def build_summary(all_results):
    """Aggregate per-modality results into a summary DataFrame.

    Infers scope ('Full' / 'EDTA') from label prefix. Skips None results
    (failed modalities). For high-dim modalities, reports PCA components
    as total_features and n_original_features as additional column.

    Parameters
    ----------
    all_results : list of dict or None
        Each result dict from run_modality_pipeline() or run_combined_pipeline().

    Returns
    -------
    pd.DataFrame
        Columns: Scope, Modality, macro_F1_mean, macro_F1_std,
        macro_F1_CI95_lower, macro_F1_CI95_upper, median_C,
        selected_features, total_features.
        Additional column n_original_features for high-dim modalities.
    """
    rows = []
    for r in all_results:
        if r is None:
            continue
        # Infer scope from label_prefix if provided, else 'Full'
        label = r['label']
        scope = 'EDTA' if label.startswith('EDTA ') else 'Full'
        total = r['n_pcs'] if r.get('is_high_dim') else r['n_total_features']
        row = {
            'Scope': scope,
            'Modality': label,
            'macro_F1_mean': round(r['macro_f1_mean'], 4),
            'macro_F1_std': round(r['macro_f1_std'], 4),
            'macro_F1_CI95_lower': round(r['macro_f1_mean'] - 1.96*r['macro_f1_std'], 4),
            'macro_F1_CI95_upper': round(r['macro_f1_mean'] + 1.96*r['macro_f1_std'], 4),
            'median_C': round(r['median_C'], 4),
            'selected_features': r['n_nonzero_full'],
            'total_features': total,
        }
        if r.get('is_high_dim'):
            row['n_original_features'] = r['n_total_features']
        rows.append(row)
    return pd.DataFrame(rows)


# ── Main ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Phase 1: baseline tissue classifier')
    parser.add_argument('--modalities', nargs='+', default=PHASE1_MODALITIES,
                        choices=list(PHASE1_MODALITIES),
                        help='Modalities to run (default: all Phase 1 modalities)')
    parser.add_argument('--skip-edta', action='store_true',
                        help='Skip EDTA-only sensitivity analysis')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip generating plots')
    parser.add_argument('--combine', action='store_true',
                        help='Also run combined model (fem4 + probe_meth + fragment_length)')
    args = parser.parse_args()

    log("Phase 1: Baseline tissue classifier")
    log(f"Modalities: {args.modalities}")
    log(f"{'='*60}\n")

    # Load metadata
    meta = load_metadata()
    log(f"Loaded {len(meta)} samples, {meta['Tissue'].nunique()} tissues, "
        f"{meta['Source'].nunique()} sources")

    # Chance baseline
    tissue_counts = meta['Tissue'].value_counts()
    chance_uniform = 1 / len(tissue_counts)
    log(f"Chance baseline (uniform): {chance_uniform:.4f}")
    log(f"Chance baseline (stratified): {(tissue_counts / len(meta)).mean():.4f}\n")

    # ── Full dataset (164 samples) ────────────────
    all_results = []
    for mod_name in args.modalities:
        r = run_modality_pipeline(mod_name, meta)
        if r:
            all_results.append(r)

    # ── EDTA-only sensitivity (96 samples) ────────
    meta_edta = None
    if not args.skip_edta:
        log("\n" + "="*60)
        log("EDTA-only sensitivity analysis")
        log("="*60)
        meta_edta = meta[meta['BCT'] == 'EDTA'].copy()
        log(f"EDTA subset: {len(meta_edta)} samples, "
            f"{meta_edta['Tissue'].nunique()} tissues")
        # Note: pancreas has 0 EDTA samples — will produce 5-class problem
        log(f"Tissues in EDTA: {sorted(meta_edta['Tissue'].unique())}")

        for mod_name in args.modalities:
            r = run_modality_pipeline(mod_name, meta_edta, label_prefix='EDTA ')
            if r:
                all_results.append(r)

    # ── Combined modality ─────────────────────────
    if args.combine:
        log("\n" + "="*60)
        log("Combined modality (fem4 + probe_meth + fragment_length)")
        log("="*60)
        r = run_combined_pipeline(COMBINE_MODALITIES, meta)
        if r:
            all_results.append(r)

        if not args.skip_edta and meta_edta is not None:
            r = run_combined_pipeline(COMBINE_MODALITIES, meta_edta, label_prefix='EDTA ')
            if r:
                all_results.append(r)

    # ── Summary ────────────────────────────────────
    log("\n" + "="*60)
    log("SUMMARY")
    log("="*60)
    summary = build_summary(all_results)
    if not summary.empty:
        summary_path = OUT / 'phase1_summary.csv'
        summary.to_csv(summary_path, index=False)
        log(f"\nSaved: {summary_path}")
        print()
        print(summary.to_string(index=False))
        print()

    # Summary figure
    if not args.skip_plots and not summary.empty:
        fig, ax = plt.subplots(figsize=(10, 5))
        modalities = summary['Modality'].values
        means = summary['macro_F1_mean'].values
        cis = summary['macro_F1_std'].values * 1.96
        colors = plt.cm.tab10(np.linspace(0, 0.8, len(modalities)))
        bars = ax.barh(range(len(modalities)), means, xerr=cis,
                       color=colors, edgecolor='k', linewidth=0.5,
                       capsize=3)
        ax.axvline(x=chance_uniform, color='gray', linestyle=':',
                   alpha=0.7, label=f'Chance ({chance_uniform:.3f})')
        ax.set_yticks(range(len(modalities)))
        ax.set_yticklabels(modalities, fontsize=9)
        ax.set_xlabel('Macro-F1')
        ax.set_title('Phase 1: CV Macro-F1 by Modality')
        for bar, v in zip(bars, means):
            ax.text(v + 0.01, bar.get_y() + bar.get_height()/2,
                    f'{v:.3f}', va='center', fontsize=8)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(FIGS / 'phase1_summary.png', dpi=150)
        plt.close(fig)
        log(f"Saved: {FIGS / 'phase1_summary.png'}")

    log(f"\nAll outputs in: {OUT}/")
    log(f"Models in: {MODELS}/")
    log(f"Figures in: {FIGS}/")
    log(f"Re-run: python3 scripts/run_phase1_pipeline.py")
    log("Done.")


if __name__ == '__main__':
    main()
