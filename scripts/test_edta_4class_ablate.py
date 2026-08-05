"""Isolate which modality destabilizes the combined model on EDTA 4-class."""
import sys, warnings, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import f1_score

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'scripts'))
from data_loading import load_metadata, get_modality_cfg, load_modality

OUT = BASE / 'output' / 'phase1' / 'edta_4class_test'
OUT.mkdir(parents=True, exist_ok=True)

N_SPLITS = 5
C_GRID = np.logspace(-3, 1, 6)
MAX_ITER = 5000
MAX_RESHUF = 100
RANDOM_STATE = 42

meta = load_metadata()
meta_edta = meta[meta['BCT'] == 'EDTA'].copy()
meta_4c = meta_edta[meta_edta['Tissue'] != 'healthyblood'].copy()

def run_modality(mod_name, meta):
    cfg = get_modality_cfg(mod_name)
    X, sample_ids, feat_names = load_modality(cfg, meta['TWIST_ID'].values, impute=False)
    meta_aligned = meta.set_index('TWIST_ID').loc[sample_ids].reset_index()
    
    y = meta_aligned['Tissue'].values
    sources = meta_aligned['Source'].values
    required_sources = set(sources)
    for attempt in range(MAX_RESHUF):
        seed = RANDOM_STATE + attempt
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        splits = list(skf.split(np.zeros(len(meta_aligned)), y))
        if all(required_sources.issubset(set(sources[tr_idx])) for tr_idx, _ in splits):
            break

    le = LabelEncoder()
    y_enc = le.fit_transform(meta_aligned['Tissue'].values)
    fold_f1s = []
    for tr_idx, te_idx in splits:
        X_tr, X_te = X[tr_idx], X[te_idx]
        y_tr, y_te = y_enc[tr_idx], y_enc[te_idx]
        col_mean = np.nanmean(X_tr, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        X_tr_imp = X_tr.copy()
        X_tr_imp[np.where(np.isnan(X_tr_imp))] = np.take(col_mean, np.where(np.isnan(X_tr_imp))[1])
        X_te_imp = X_te.copy()
        X_te_imp[np.where(np.isnan(X_te_imp))] = np.take(col_mean, np.where(np.isnan(X_te_imp))[1])
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_tr_imp)
        X_te_s = scaler.transform(X_te_imp)
        n_inner = min(3, np.min(np.bincount(y_tr)))
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=RANDOM_STATE)
        clf = LogisticRegression(l1_ratio=1, solver='saga', class_weight='balanced', max_iter=MAX_ITER, random_state=RANDOM_STATE)
        gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
        gs.fit(X_tr_s, y_tr)
        model = gs.best_estimator_
        y_pred = model.predict(X_te_s)
        fold_f1s.append(f1_score(y_te, y_pred, average='macro'))
    return np.mean(fold_f1s), np.std(fold_f1s)

def run_combined(mod_names, meta):
    X_raw, first_ids = {}, None
    for mod_name in mod_names:
        cfg = get_modality_cfg(mod_name)
        X_mod, ids_mod, _ = load_modality(cfg, meta['TWIST_ID'].values, impute=False)
        if first_ids is None:
            first_ids = ids_mod
        if not np.array_equal(ids_mod, first_ids):
            order_map = {s: i for i, s in enumerate(first_ids)}
            X_mod = X_mod[np.array([order_map[s] for s in ids_mod])]
        X_raw[mod_name] = X_mod
    meta_aligned = meta.set_index('TWIST_ID').loc[first_ids].reset_index()
    
    y = meta_aligned['Tissue'].values
    sources = meta_aligned['Source'].values
    required_sources = set(sources)
    for attempt in range(MAX_RESHUF):
        seed = RANDOM_STATE + attempt
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        splits = list(skf.split(np.zeros(len(meta_aligned)), y))
        if all(required_sources.issubset(set(sources[tr_idx])) for tr_idx, _ in splits):
            break
    
    le = LabelEncoder()
    y_enc = le.fit_transform(meta_aligned['Tissue'].values)
    fold_f1s = []
    for tr_idx, te_idx in splits:
        fold_Xs_tr, fold_Xs_te = [], []
        for mod_name in mod_names:
            Xm = X_raw[mod_name]
            Xm_tr, Xm_te = Xm[tr_idx], Xm[te_idx]
            col_mean = np.nanmean(Xm_tr, axis=0)
            col_mean = np.nan_to_num(col_mean, nan=0.0)
            Xm_tr_imp = Xm_tr.copy()
            Xm_tr_imp[np.where(np.isnan(Xm_tr_imp))] = np.take(col_mean, np.where(np.isnan(Xm_tr_imp))[1])
            Xm_te_imp = Xm_te.copy()
            Xm_te_imp[np.where(np.isnan(Xm_te_imp))] = np.take(col_mean, np.where(np.isnan(Xm_te_imp))[1])
            scaler = StandardScaler()
            fold_Xs_tr.append(scaler.fit_transform(Xm_tr_imp))
            fold_Xs_te.append(scaler.transform(Xm_te_imp))
        X_tr = np.hstack(fold_Xs_tr)
        X_te = np.hstack(fold_Xs_te)
        y_tr, y_te = y_enc[tr_idx], y_enc[te_idx]
        n_inner = min(3, np.min(np.bincount(y_tr)))
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=RANDOM_STATE)
        clf = LogisticRegression(l1_ratio=1, solver='saga', class_weight='balanced', max_iter=MAX_ITER, random_state=RANDOM_STATE)
        gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
        gs.fit(X_tr, y_tr)
        y_pred = gs.predict(X_te)
        fold_f1s.append(f1_score(y_te, y_pred, average='macro'))
    return np.mean(fold_f1s), np.std(fold_f1s)

# Single modalities
for mod in ['fem4', 'probe_meth', 'fragment_length']:
    m, s = run_modality(mod, meta_4c)
    print(f"{mod:20s}: F1 = {m:.4f} ± {s:.4f}")

# Pairs
print()
for pair in [['fem4','probe_meth'], ['fem4','fragment_length'], ['probe_meth','fragment_length']]:
    m, s = run_combined(pair, meta_4c)
    print(f"{' + '.join(pair):30s}: F1 = {m:.4f} ± {s:.4f}")
