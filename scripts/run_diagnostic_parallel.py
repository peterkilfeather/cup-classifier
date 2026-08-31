#!/usr/bin/env python3
"""
Diagnostic Protocol — Parallel execution across modalities.
Each modality-scope combination runs independently via ProcessPoolExecutor.
"""

import sys, os, warnings, subprocess, time, multiprocessing as mp
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV, train_test_split
from sklearn.metrics import f1_score
from sklearn.utils import shuffle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path('/xscratch/farney/cup-classifier')
INPUT = BASE / 'input'
META = INPUT / 'metadata' / 'metadata_cleaned.csv'
OUT = BASE / 'output' / 'diagnostic-protocol'
OUT.mkdir(parents=True, exist_ok=True)
TMP = OUT / '.tmp'
TMP.mkdir(exist_ok=True)

COVARIATES = ['Tissue', 'Source', 'BCT', 'Sex']
PALETTES = {
    'Tissue': {'healthyblood': '#4daf4a', 'colon': '#e41a1c', 'liver': '#377eb8',
               'pancreas': '#ff7f00', 'prostate': '#984ea3', 'stomach': '#f781bf'},
    'Source': {'Fox Chase': '#e41a1c', 'Audubon': '#377eb8',
               'Sowalsky': '#4daf4a', 'NIH_Clinical_Center': '#ff7f00'},
    'BCT': {'EDTA': '#377eb8', 'Citrate': '#e41a1c', 'Streck': '#4daf4a', 'ACD': '#ff7f00'},
    'Sex': {'Male': '#377eb8', 'Female': '#e41a1c'},
}

def log(s):
    print(f"[{mp.current_process().name} {datetime.now().strftime('%H:%M:%S')}] {s}", flush=True)

# ── Data Loading ───────────────────────────────────────

def load_metadata():
    meta = pd.read_csv(META)
    meta['Tissue'] = meta['Tissue'].str.lower()
    meta['Sex'] = meta['Sex'].str.strip()
    return meta

def get_modality_cfg(mod_name):
    """Return config dict for a modality name."""
    configs = {
        'probe_meth': {
            'file': INPUT / 'methylation' / 'probe_meth' / 'all_samples.probe_meth_enriched_filtered.mNonCpGlt4_frac.wide.tsv',
            'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
            'high_dim': False, 'label': 'Methylation (probe-avg)',
        },
        'fem4': {
            'file': INPUT / 'fragmentomic' / 'ALL_fem4_features.tsv',
            'sep': '\t', 'sample_col': 'sample', 'drop_cols': ['tissue', 'N_motifs'],
            'high_dim': False, 'label': 'FEM4 (256)',
        },
        'fragment_length': {
            'file': INPUT / 'fragmentomic' / 'fragment_length_features_qc.csv',
            'sep': ',', 'sample_col': 'Sample_ID', 'drop_cols': ['Tissue'],
            'high_dim': False, 'label': 'Fragment length (369)',
        },
        'cnvkit': {
            'file': INPUT / 'cnvkit' / 'cnvkit_features_thr0.10.tsv',
            'sep': '\t', 'sample_col': 'sample', 'drop_cols': ['cns_path', 'tissue', 'cnr_path'],
            'high_dim': False, 'label': 'CNVkit thr0.10',
        },
        'probe_cpg': {
            'file': INPUT / 'methylation' / 'probe_cpg' / 'all_samples.probe_cpg_unenriched_filtered.mNonCpGlt4_frac.wide.tsv',
            'sep': '\t', 'sample_col': 'sample', 'drop_cols': [],
            'high_dim': True, 'label': 'Per-CpG methylation (32K->PCs)',
        },
        'end_density': {
            'file': INPUT / 'fragmentomic' / 'end_density' / 'all_samples_ends_100kb_CPM_matrix.tsv',
            'sep': '\t', 'sample_col': None, 'drop_cols': [],
            'high_dim': True, 'label': 'End density (31K bins->PCs)',
        },
    }
    return configs[mod_name]

def load_modality(cfg, meta_samples):
    fpath = cfg['file']
    if cfg['sample_col'] is None:
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
        # Deduplicate: keep first occurrence per sample ID (matches data_loading.py)
        dup_mask = df.index.duplicated(keep='first')
        if dup_mask.any():
            log(f"  Dropping {dup_mask.sum()} duplicate sample rows")
            df = df[~dup_mask]
        X = df.values.astype(np.float64)
        sample_ids = df.index.values
        feat_names = df.columns.values

    mask = np.isin(sample_ids, meta_samples)
    X, sample_ids = X[mask], sample_ids[mask]

    nan_feats = np.all(np.isnan(X), axis=0)
    if nan_feats.any():
        log(f"  Dropping {nan_feats.sum()} all-NaN features")
        X = X[:, ~nan_feats]
        if feat_names is not None:
            feat_names = feat_names[~nan_feats]

    if not cfg.get('high_dim', False):
        col_mean = np.nanmean(X, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(col_mean, inds[1])

    return X, sample_ids, feat_names

def elbow_k(X, max_k=30, min_k=5):
    """Determine number of PCs via elbow of scree plot, with minimum floor.
    
    The second-derivative elbow finds the sharpest drop in singular values.
    For high-dim biological data, this can be too aggressive (K=1-2).
    `min_k` ensures enough PCs to capture distributed biological signal.
    Capped at min(max_k, 20, n-1).
    """
    n = min(X.shape)
    max_k = min(max_k, n - 1)
    if max_k < 2:
        return max(min_k, n - 1)
    Xc = X.copy()
    col_mean = np.nanmean(Xc, axis=0)
    col_mean = np.nan_to_num(col_mean, nan=0.0)
    inds = np.where(np.isnan(Xc))
    Xc[inds] = np.take(col_mean, inds[1])
    from sklearn.utils.extmath import randomized_svd
    U, s, Vt = randomized_svd(Xc, n_components=max_k, random_state=42)
    d2 = np.abs(np.diff(s, 2))
    if len(d2) == 0:
        return min(max(max_k, min_k), 20, n - 1)
    elbow = np.argmax(d2) + 1
    return min(max(elbow + 1, min_k), 20, n - 1)

# ── 1. PCA ─────────────────────────────────────────────

def run_pca(X, meta, modality_label, scope_label):
    log(f"PCA...")
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    n_comp = min(Xs.shape[0], Xs.shape[1], 10)
    pca = PCA(n_components=n_comp)
    scores = pca.fit_transform(Xs)
    var_exp = pca.explained_variance_ratio_

    safe = modality_label.replace(' ', '_').replace('->', 'to')
    sc = scope_label.replace(' ', '_')

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, n_comp+1), var_exp*100, 'o-', color='#333')
    ax.set_xlabel('PC'); ax.set_ylabel('% Variance')
    ax.set_title(f'{modality_label} - {scope_label}')
    fig.tight_layout()
    fig.savefig(OUT / f'scree_{safe}_{sc}.png', dpi=150)
    plt.close(fig)

    for cov in ['Tissue', 'Source', 'BCT', 'Sex']:
        fig, ax = plt.subplots(figsize=(7, 5.5))
        pal = PALETTES.get(cov, {})
        for val in sorted(meta[cov].unique()):
            mask = meta[cov] == val
            c = pal.get(val, '#999')
            ax.scatter(scores[mask, 0], scores[mask, 1], c=c, label=val,
                      s=40, edgecolors='k', linewidth=0.3, alpha=0.8)
        ax.set_xlabel(f'PC1 ({var_exp[0]*100:.1f}%)')
        ax.set_ylabel(f'PC2 ({var_exp[1]*100:.1f}%)')
        ax.set_title(f'{modality_label} - colored by {cov}')
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(OUT / f'pca_{safe}_{cov}_{sc}.png', dpi=150)
        plt.close(fig)

    return var_exp

# ── 2. PERMANOVA ───────────────────────────────────────

def run_permanova(X, meta, scope_label, modality_label):
    log(f"PERMANOVA...")
    n = X.shape[0]
    if n < 10:
        return None

    Xs = StandardScaler().fit_transform(X)

    D = squareform(pdist(Xs, 'euclidean'))

    safe = f"{modality_label.replace(' ', '_').replace('->', 'to')}_{scope_label.replace(' ', '_')}"
    td = TMP / safe
    td.mkdir(exist_ok=True)

    np.savetxt(td / 'dist.txt', D)
    meta[COVARIATES].to_csv(td / 'meta.csv', index=False)

    rscript = f'''
    library(vegan)
    dist <- as.matrix(read.table("{td}/dist.txt"))
    meta <- read.csv("{td}/meta.csv")
    for (c in colnames(meta)) meta[[c]] <- as.factor(meta[[c]])
    # Drop single-level factors (e.g., BCT in EDTA subset)
    keep <- sapply(meta, function(x) nlevels(x) > 1)
    keep_names <- names(keep)[keep]
    if (length(keep_names) < 1) {{ cat("No multi-level factors\n"); q() }}
    formula_str <- paste("dist ~", paste(keep_names, collapse=" + "))
    fit <- adonis2(as.formula(formula_str), data=meta, by="margin", permutations=999)
    write.table(as.data.frame(fit), file="{td}/marginal.tsv", sep="\t", quote=FALSE)
    fit_seq <- adonis2(as.formula(formula_str), data=meta, by="terms", permutations=999)
    write.table(as.data.frame(fit_seq), file="{td}/sequential.tsv", sep="\t", quote=FALSE)
    # Also run with BCT-first ordering if BCT is present
    if ("BCT" %in% keep_names) {{
        others <- setdiff(keep_names, "BCT")
        bct_first_str <- paste("dist ~ BCT +", paste(others, collapse=" + "))
        fit_bct <- adonis2(as.formula(bct_first_str), data=meta, by="terms", permutations=999)
        write.table(as.data.frame(fit_bct), file="{td}/bct_first.tsv", sep="\t", quote=FALSE)
    }}
    '''

    with open(td / 'run.R', 'w') as f:
        f.write(rscript)

    try:
        result = subprocess.run(
            ['Rscript', str(td / 'run.R')],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode != 0:
            log(f"  R failed: {result.stderr[:200]}")
            return None
    except subprocess.TimeoutExpired:
        log(f"  R timed out")
        return None

    results = {}
    for key, fname in [('marginal', 'marginal.tsv'),
                        ('sequential', 'sequential.tsv'),
                        ('bct_first', 'bct_first.tsv')]:
        fpath = td / fname
        if not fpath.exists():
            continue
        try:
            df = pd.read_csv(fpath, sep='\t', index_col=0)
            results[key] = {}
            for t in COVARIATES:
                if t in df.index:
                    results[key][t] = {
                        'R2': float(df.loc[t, 'R2']),
                        'F': float(df.loc[t, 'F']),
                        'p': float(df.loc[t, 'Pr(>F)']),
                    }
        except Exception as e:
            log(f"  Parse error {fname}: {e}")

    return results

# ── 3. Classifier ──────────────────────────────────────

def run_single_classifier(X, y):
    n_classes = len(np.unique(y))
    if n_classes <= 1:
        return 0.0

    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    Xs = StandardScaler().fit_transform(X)

    clf = LogisticRegression(l1_ratio=1, solver='saga', max_iter=5000,
                             random_state=42, class_weight='balanced')
    param_grid = {'C': np.logspace(-3, 1, 6)}

    n_splits = min(5, np.min(np.bincount(y_enc)))
    if n_splits < 2:
        return _fallback_classifier(Xs, y_enc)

    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = []

    for train_idx, test_idx in outer_cv.split(Xs, y_enc):
        X_tr, X_te = Xs[train_idx], Xs[test_idx]
        y_tr, y_te = y_enc[train_idx], y_enc[test_idx]

        inner_n = min(3, np.min(np.bincount(y_tr)))
        if inner_n < 2:
            gs = clf
            gs.fit(X_tr, y_tr)
        else:
            inner_cv = StratifiedKFold(n_splits=inner_n, shuffle=True, random_state=42)
            gs = GridSearchCV(clf, param_grid, cv=inner_cv, scoring='f1_macro')
            gs.fit(X_tr, y_tr)
            gs = gs.best_estimator_

        y_pred = gs.predict(X_te)
        scores.append(f1_score(y_te, y_pred, average='macro'))

    return np.mean(scores) if scores else 0.0

def _fallback_classifier(X, y):
    X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
    clf = LogisticRegression(l1_ratio=1, solver='saga', max_iter=5000,
                             C=0.1, random_state=42, class_weight='balanced')
    clf.fit(X_tr, y_tr)
    y_pred = clf.predict(X_te)
    return f1_score(y_te, y_pred, average='macro')

def run_classifier_analysis(X, meta, modality_label, scope_label):
    log(f"Classifier...")

    targets = [('Tissue', 'Tissue'), ('BCT', 'BCT'), ('Source', 'Source'),
               ('Shuffled_Tissue', 'Tissue')]

    f1_scores = {}
    for name, col in targets:
        y = meta[col].values
        if name == 'Shuffled_Tissue':
            y = shuffle(y, random_state=42)
        f1_scores[name] = run_single_classifier(X, y)

    chance = f1_scores['Shuffled_Tissue']

    fig, ax = plt.subplots(figsize=(6, 4))
    targets_plot = ['Tissue', 'BCT', 'Source', 'Shuffled_Tissue']
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#cccccc']
    vals = [f1_scores.get(t, 0) for t in targets_plot]
    bars = ax.bar(targets_plot, vals, color=colors, edgecolor='k', linewidth=0.5)
    ax.axhline(y=vals[0]*0.8, color='red', linestyle=':', alpha=0.5, label=f'0.8xF1_T={vals[0]*0.8:.3f}')
    ax.set_ylabel('Macro F1')
    ax.set_title(f'{modality_label} - {scope_label}')
    ax.legend(fontsize=7)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.01,
                f'{v:.3f}', ha='center', va='bottom', fontsize=8)
    ax.set_ylim(0, max(vals + [chance, 0.1]) * 1.25 + 0.05)
    fig.tight_layout()
    safe = modality_label.replace(' ', '_').replace('->', 'to')
    fig.savefig(OUT / f'classifier_{safe}_{scope_label.replace(" ","_")}.png', dpi=150)
    plt.close(fig)

    for t in targets_plot:
        log(f"  {t}: F1={f1_scores.get(t,0):.3f} (chance={chance:.3f})")

    return f1_scores, chance

# ── 4. Variance Partitioning ───────────────────────────

def run_variance_partitioning(X, meta, feat_names, modality_label, scope_label):
    log(f"Variance partitioning ({X.shape[1]} features)...")
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm

    results = {c: [] for c in COVARIATES}
    n = X.shape[1]
    step = max(1, n // 10)

    # Build a data frame with features + covariates for formula API
    # Formula API needs the feature as a named column
    cov_df_base = pd.DataFrame({c: meta[c].values for c in COVARIATES})

    # Drop single-level factors (e.g., BCT in EDTA subset) — they break OLS
    active_covariates = [c for c in COVARIATES if meta[c].nunique() > 1]
    inactive = [c for c in COVARIATES if c not in active_covariates]

    t0 = time.time()
    for i in range(n):
        if i % step == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            log(f"  {i}/{n} ({rate:.0f} feats/s)")

        y = X[:, i]
        mask = ~np.isnan(y)
        if mask.sum() < 10:
            for c in COVARIATES:
                results[c].append(0.0)
            continue

        try:
            # Build data frame for this feature
            y_name = 'feat'
            data_i = cov_df_base[mask].copy()
            data_i[y_name] = y[mask]
            formula = f'{y_name} ~ {" + ".join(f"C({c})" for c in active_covariates)}'
            model = smf.ols(formula, data=data_i).fit()
            anova = anova_lm(model, typ=2)
            ss_total = anova['sum_sq'].sum()
            # Active covariates: read from ANOVA table
            for c in active_covariates:
                anova_idx = f'C({c})'
                if anova_idx in anova.index:
                    ss = anova.loc[anova_idx, 'sum_sq']
                    results[c].append(ss / ss_total if ss_total > 0 else 0)
                else:
                    results[c].append(0.0)
            # Inactive (single-level) covariates: zero by definition
            for c in inactive:
                results[c].append(0.0)
        except Exception as e:
            # If model fails, record zeros
            for c in COVARIATES:
                results[c].append(0.0)

    fig, ax = plt.subplots(figsize=(8, 5))
    data_violin = []
    colors_violin = []
    for cov in COVARIATES:
        v = np.array(results[cov])
        v = np.clip(v, 0, np.percentile(v, 98)*1.5)
        data_violin.append(v)
        pal = list(PALETTES.get(cov, {'x': '#999'}).values())
        colors_violin.append(pal[0] if pal else '#999')

    parts = ax.violinplot(data_violin, showmeans=True, showmedians=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors_violin[i]); pc.set_alpha(0.6)
    parts['cmeans'].set_color('#333')

    ax.set_xticks(range(1, len(COVARIATES)+1))
    ax.set_xticklabels(COVARIATES, fontsize=9)
    ax.set_ylabel('Variance fraction (SS / total SS)')
    ax.set_title(f'{modality_label} - {scope_label}')
    ax.axhline(y=0.10, color='red', linestyle='--', alpha=0.5)
    fig.tight_layout()
    safe = modality_label.replace(' ', '_').replace('->', 'to')
    fig.savefig(OUT / f'varpart_{safe}_{scope_label.replace(" ","_")}.png', dpi=150)
    plt.close(fig)

    summary = {}
    for cov in COVARIATES:
        v = np.array(results[cov])
        summary[cov] = {
            'mean': float(np.mean(v)),
            'median': float(np.median(v)),
            'pct_gt_0.10': float(np.mean(v > 0.10) * 100),
            'pct_gt_0.25': float(np.mean(v > 0.25) * 100),
        }
    return results, summary

# ── Single Modality Worker ─────────────────────────────

def process_modality(args):
    """Process one modality-scope combination. Returns result dict."""
    mod_name, scope_label = args

    meta = pd.read_csv(META)
    meta['Tissue'] = meta['Tissue'].str.lower()
    meta['Sex'] = meta['Sex'].str.strip()

    # Filter to EDTA-only for EDTA scope
    if 'EDTA' in scope_label:
        meta = meta[meta['BCT'] == 'EDTA'].copy()
        if len(meta) < 10:
            return {'mod_name': mod_name, 'scope': scope_label,
                    'error': f'Too few EDTA samples: {len(meta)}', 'result': None}

    cfg = get_modality_cfg(mod_name)

    log(f"Starting {cfg['label']} ({scope_label}) with {len(meta)} samples")

    try:
        X, sample_ids, feat_names = load_modality(cfg, meta['TWIST_ID'].values)
    except Exception as e:
        log(f"FAILED load: {e}")
        return {'mod_name': mod_name, 'scope': scope_label, 'error': str(e), 'result': None}

    if X.shape[0] < 5 or X.shape[1] < 2:
        return {'mod_name': mod_name, 'scope': scope_label,
                'error': f'Too few: {X.shape}', 'result': None}

    meta_aligned = meta.set_index('TWIST_ID').loc[sample_ids].reset_index()

    # High-dim reduction
    if cfg.get('high_dim', False):
        log(f"Reducing {X.shape[1]} features via PCA...")
        X_imp = X.copy()
        col_mean = np.nanmean(X_imp, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        inds = np.where(np.isnan(X_imp))
        X_imp[inds] = np.take(col_mean, inds[1])
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_imp)
        k = elbow_k(X_scaled)
        log(f"Elbow K = {k}")
        pca = PCA(n_components=k)
        X_pc = pca.fit_transform(X_scaled)
        log(f"PC1: {pca.explained_variance_ratio_[0]*100:.1f}%, Cum: {pca.explained_variance_ratio_.sum()*100:.1f}%")
        X_clf = X_pc
    else:
        X_clf = X

    result = {'mod_name': mod_name, 'scope': scope_label, 'error': None,
              'label': cfg['label'], 'classifier': None, 'varpart': None,
              'permanova': None, 'var_exp': None}

    # 1. PCA
    try:
        result['var_exp'] = run_pca(X if not cfg.get('high_dim') else X_clf,
                                    meta_aligned, cfg['label'], scope_label)
    except Exception as e:
        log(f"PCA failed: {e}")

    # 2. PERMANOVA
    try:
        perm_res = run_permanova(X_clf, meta_aligned, scope_label, cfg['label'])
        result['permanova'] = perm_res
        if perm_res:
            m = perm_res.get('marginal', {})
            log(f"Marginal R2: T={m.get('Tissue',{}).get('R2',0):.4f} "
                f"B={m.get('BCT',{}).get('R2',0):.4f} S={m.get('Source',{}).get('R2',0):.4f}")
    except Exception as e:
        log(f"PERMANOVA failed: {e}")

    # 3. Classifier
    try:
        f1_scores, chance = run_classifier_analysis(X_clf, meta_aligned, cfg['label'], scope_label)
        result['classifier'] = f1_scores
        result['chance'] = chance
    except Exception as e:
        log(f"Classifier failed: {e}")

    # 4. Variance partitioning (on raw features, subsampled for high-dim)
    try:
        vp_X = X
        vp_feats = feat_names
        if cfg.get('high_dim', False) and X.shape[1] > 5000:
            np.random.seed(42)
            idx = np.random.choice(X.shape[1], min(5000, X.shape[1]), replace=False)
            vp_X = X[:, idx]
            vp_feats = feat_names[idx] if feat_names is not None else None
            log(f"Subsampled to {vp_X.shape[1]} features for variance partitioning")
        vp_res, vp_summary = run_variance_partitioning(
            vp_X, meta_aligned, vp_feats, cfg['label'], scope_label)
        result['varpart'] = vp_summary
    except Exception as e:
        log(f"Variance partitioning failed: {e}")

    log(f"Done")
    return result

# ── Adjudication ───────────────────────────────────────

def build_adjudication(all_results, scope_label):
    rows = []
    for r in all_results:
        if r is None or r.get('error') or r.get('result') is None and r.get('classifier') is None:
            continue
        label = r.get('label', '?')
        perm = (r.get('permanova') or {}).get('marginal', {})
        r2_tissue = perm.get('Tissue', {}).get('R2', 0) or 0
        r2_bct = perm.get('BCT', {}).get('R2', 0) or 0
        r2_source = perm.get('Source', {}).get('R2', 0) or 0

        clf = r.get('classifier') or {}
        f1_tissue = clf.get('Tissue', 0) or 0
        f1_bct = clf.get('BCT', 0) or 0
        f1_source = clf.get('Source', 0) or 0
        chance = r.get('chance', 0) or 0
        thresh = f1_tissue * 0.8

        vp = r.get('varpart') or {}
        vp_t = vp.get('Tissue', {}).get('mean', 0) if isinstance(vp.get('Tissue'), dict) else 0
        vp_b = vp.get('BCT', {}).get('mean', 0) if isinstance(vp.get('BCT'), dict) else 0

        def verdict(r2, f1, target_f1, chance):
            if r2 >= 0.10 and f1 >= target_f1: return 'DOMINANT'
            if r2 >= 0.10: return 'Material'
            if f1 >= target_f1: return 'Borderline'
            if f1 < 2 * chance: return 'Weak'
            return 'Detectable'

        rows.append({
            'Modality': label,
            'R2_Tissue': round(r2_tissue, 4), 'R2_BCT': round(r2_bct, 4),
            'R2_Source': round(r2_source, 4),
            'F1_Tissue': round(f1_tissue, 4), 'F1_BCT': round(f1_bct, 4),
            'F1_Source': round(f1_source, 4),
            'Chance': round(chance, 4), '0.8xF1_Tissue': round(thresh, 4),
            'BCT_Verdict': verdict(r2_bct, f1_bct, thresh, chance),
            'Source_Verdict': verdict(r2_source, f1_source, thresh, chance),
            'VP_Tissue': round(vp_t, 4), 'VP_BCT': round(vp_b, 4),
        })

    if not rows:
        return None
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f'adjudication_{scope_label.replace(" ","_")}.csv', index=False)

    print(f"\n  ADJUDICATION - {scope_label}")
    for _, row in df.iterrows():
        print(f"  {row['Modality']:35s}  BCT:{row['BCT_Verdict']:15s}  Source:{row['Source_Verdict']:15s}  "
              f"R2_T={row['R2_Tissue']:.3f} R2_B={row['R2_BCT']:.3f}  F1_T={row['F1_Tissue']:.3f} F1_B={row['F1_BCT']:.3f}")
    return df

# ── Main ───────────────────────────────────────────────

if __name__ == '__main__':
    print("Loading metadata...")
    meta = load_metadata()
    print(f"  {len(meta)} samples, {meta['Tissue'].nunique()} tissues")

    MODALITY_NAMES = ['probe_meth', 'fem4', 'fragment_length', 'cnvkit', 'probe_cpg', 'end_density']

    # Build task list: (mod_name, scope_label)
    tasks = [(m, 'Full (164)') for m in MODALITY_NAMES]
    tasks += [(m, 'EDTA (96)') for m in MODALITY_NAMES]

    print(f"Launching {len(tasks)} tasks with up to 4 parallel workers...")

    n_workers = min(4, mp.cpu_count())
    all_results = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_modality, t): t for t in tasks}
        for future in as_completed(futures):
            t = futures[future]
            try:
                result = future.result()
                all_results.append(result)
                status = 'OK' if not result.get('error') else f"ERR: {result['error']}"
                print(f"  [{result.get('mod_name','?')} {t[1]}] {status}")
            except Exception as e:
                print(f"  [{t[0]} {t[1]}] CRASH: {e}")

    # Separate results by scope
    for scope in ['Full (164)', 'EDTA (96)']:
        scope_results = [r for r in all_results if r and r.get('scope') == scope]
        build_adjudication(scope_results, scope)

    # Generate summary
    print(f"\n{'='*60}")
    print(f"  All output in: {OUT}/")
    print(f"  Re-run: python3 scripts/run_diagnostic_parallel.py")
    print(f"{'='*60}")
