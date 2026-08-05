"""Diagnose why probe_meth collapses: check per-fold confusion patterns."""
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import f1_score, confusion_matrix

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'scripts'))
from data_loading import load_metadata, get_modality_cfg, load_modality

N_SPLITS = 5
C_GRID = np.logspace(-3, 1, 6)
MAX_ITER = 5000
MAX_RESHUF = 100
RANDOM_STATE = 42

meta = load_metadata()
meta_edta = meta[meta['BCT'] == 'EDTA'].copy()
meta_4c = meta_edta[meta_edta['Tissue'] != 'healthyblood'].copy()
cfg = get_modality_cfg('probe_meth')
X, ids, _ = load_modality(cfg, meta_4c['TWIST_ID'].values, impute=False)
meta_aligned = meta_4c.set_index('TWIST_ID').loc[ids].reset_index()

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

for fold, (tr_idx, te_idx) in enumerate(splits):
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
    
    y_pred = gs.predict(X_te_s)
    
    cm = confusion_matrix(y_te, y_pred)
    
    print(f"\n=== Fold {fold+1} (F1={f1_score(y_te, y_pred, average='macro'):.4f}) ===")
    print(f"  Test set:")
    for _, row in meta_aligned.iloc[te_idx].iterrows():
        true_label = row['Tissue']
        pred_label = le.inverse_transform([y_pred[len(meta_aligned.iloc[te_idx].index) - len(meta_aligned.iloc[te_idx].index) + list(te_idx).index(_)]])[0]
    # Better: map predictions back
    te_df = meta_aligned.iloc[te_idx].copy()
    te_df['pred'] = le.inverse_transform(y_pred)
    te_df['correct'] = te_df['Tissue'] == te_df['pred']
    for _, row in te_df.iterrows():
        mark = '✓' if row['correct'] else '✗'
        print(f"    {mark} true={row['Tissue']:12s} src={row['Source']:25s} pred={row['pred']:12s}")
    
    # Confusion matrix as table
    print(f"\n  Confusion matrix (rows=true, cols=pred):")
    print(f"    {'':>12s} {' '.join(f'{c:12s}' for c in le.classes_)}")
    for i, true_cls in enumerate(le.classes_):
        row_str = ' '.join(f'{cm[i,j]:>12d}' for j in range(len(le.classes_)))
        print(f"    {true_cls:>12s} {row_str}")
