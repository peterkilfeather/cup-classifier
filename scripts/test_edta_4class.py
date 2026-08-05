"""Quick test: EDTA 4-class (cancer-only) using Phase 1 classifier."""
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

# Load metadata, filter to EDTA + cancer-only
meta = load_metadata()
meta_edta = meta[meta['BCT'] == 'EDTA'].copy()
meta_4c = meta_edta[meta_edta['Tissue'] != 'healthyblood'].copy()
print(f"EDTA 4-class: {len(meta_4c)} samples, "
      f"{meta_4c['Tissue'].nunique()} tissues: {sorted(meta_4c['Tissue'].unique())}")
print(meta_4c['Tissue'].value_counts().to_string())

# -- FEM4 only --
for mod_name in ['fem4']:  # could extend to ['fem4','probe_meth','fragment_length']
    cfg = get_modality_cfg(mod_name)
    X, sample_ids, feat_names = load_modality(cfg, meta_4c['TWIST_ID'].values, impute=False)
    meta_aligned = meta_4c.set_index('TWIST_ID').loc[sample_ids].reset_index()

    # Source-covering CV (same as Phase 1)
    y = meta_aligned['Tissue'].values
    sources = meta_aligned['Source'].values
    required_sources = set(sources)
    for attempt in range(MAX_RESHUF):
        seed = RANDOM_STATE + attempt
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        splits = list(skf.split(np.zeros(len(meta_aligned)), y))
        ok = True
        for tr_idx, _ in splits:
            if not required_sources.issubset(set(sources[tr_idx])):
                ok = False
                break
        if ok:
            break

    le = LabelEncoder()
    y_enc = le.fit_transform(meta_aligned['Tissue'].values)

    fold_f1s = []
    for fold, (tr_idx, te_idx) in enumerate(splits):
        X_tr, X_te = X[tr_idx], X[te_idx]
        y_tr, y_te = y_enc[tr_idx], y_enc[te_idx]

        # Per-fold imputation
        col_mean = np.nanmean(X_tr, axis=0)
        col_mean = np.nan_to_num(col_mean, nan=0.0)
        X_tr_imp = X_tr.copy()
        X_tr_imp[np.where(np.isnan(X_tr_imp))] = np.take(col_mean, np.where(np.isnan(X_tr_imp))[1])
        X_te_imp = X_te.copy()
        X_te_imp[np.where(np.isnan(X_te_imp))] = np.take(col_mean, np.where(np.isnan(X_te_imp))[1])

        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_tr_imp)
        X_te_s = scaler.transform(X_te_imp)

        # Inner CV for C
        n_inner = min(3, np.min(np.bincount(y_tr)))
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=RANDOM_STATE)
        clf = LogisticRegression(l1_ratio=1, solver='saga', class_weight='balanced',
                                 max_iter=MAX_ITER, random_state=RANDOM_STATE)
        gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
        gs.fit(X_tr_s, y_tr)
        model = gs.best_estimator_

        y_pred = model.predict(X_te_s)
        f1 = f1_score(y_te, y_pred, average='macro')
        fold_f1s.append(f1)
        print(f"  Fold {fold+1}: C={model.C:.4f}, F1={f1:.4f}")

    f1_mean = np.mean(fold_f1s)
    f1_std = np.std(fold_f1s)
    print(f"\nEDTA 4-class {mod_name}: F1 = {f1_mean:.4f} ± {f1_std:.4f}")
    print(f"  (vs other repo's FEM4-only on 4-class EDTA: 0.7596 ± 0.0368)")
    print(f"  (vs current repo's 5-class EDTA FEM4: 0.7153 ± 0.1365)")

    # Save results
    result = {
        'modality': mod_name,
        'n_samples': len(meta_aligned),
        'n_features': X.shape[1],
        'n_classes': len(le.classes_),
        'classes': le.classes_.tolist(),
        'macro_f1_mean': f1_mean,
        'macro_f1_std': f1_std,
        'per_fold_f1': fold_f1s,
        'per_fold_C': [model.C for model in [gs.best_estimator_]],  # last fold's model
        'classifier': 'LogisticRegression(L1, saga, multinomial)',
        'cv': f'{N_SPLITS}-fold StratifiedKFold(Tissue, source_coverage_check)',
    }
    out_path = OUT / f'{mod_name}_4class_edta.json'
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"Saved: {out_path}")
