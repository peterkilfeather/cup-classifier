"""Quick test: EDTA 4-class (cancer-only) — combined modalities."""
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

# Modalities to combine (matching other repo's "3-family": FEM4 + Fragment + CNVkit)
# But our CNVkit path may differ; use what we have: fem4 + probe_meth + fragment_length
COMBINE_MODALITIES = ['fem4', 'probe_meth', 'fragment_length']

def run_combined_4class(mod_names, meta, label=''):
    """Combined modality pipeline for 4-class EDTA, matching Phase 1 approach."""
    # Load raw features for each modality
    X_raw = {}
    all_feat_names = {}
    n_features_list = []
    first_ids = None

    for mod_name in mod_names:
        cfg = get_modality_cfg(mod_name)
        X_mod, ids_mod, feats_mod = load_modality(cfg, meta['TWIST_ID'].values, impute=False)
        if first_ids is None:
            first_ids = ids_mod
        # Align rows to first modality's order
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
    print(f"\nCombined {label}: {dict(zip(mod_names, n_features_list))}, total={total_n} features")

    meta_aligned = meta.set_index('TWIST_ID').loc[first_ids].reset_index()

    # Source-covering CV
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
            print(f"  Source-covering CV found after {attempt} retries")
            break

    le = LabelEncoder()
    y_enc = le.fit_transform(meta_aligned['Tissue'].values)

    fold_f1s = []
    fold_Cs = []
    for fold, (tr_idx, te_idx) in enumerate(splits):
        # Scale each modality per-fold, then concatenate
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
        clf = LogisticRegression(l1_ratio=1, solver='saga', class_weight='balanced',
                                 max_iter=MAX_ITER, random_state=RANDOM_STATE)
        gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
        gs.fit(X_tr, y_tr)
        model = gs.best_estimator_

        y_pred = model.predict(X_te)
        f1 = f1_score(y_te, y_pred, average='macro')
        fold_f1s.append(f1)
        fold_Cs.append(model.C)
        n_nonzero = np.any(model.coef_ != 0, axis=0).sum()
        print(f"  Fold {fold+1}: C={model.C:.4f}, F1={f1:.4f}, nonzero={n_nonzero}/{total_n}")

    f1_mean = np.mean(fold_f1s)
    f1_std = np.std(fold_f1s)
    print(f"\nEDTA 4-class {label}: F1 = {f1_mean:.4f} ± {f1_std:.4f}")
    print(f"  (vs other repo's 3-family 4-class EDTA: 0.7395 ± 0.0345)")

    result = {
        'label': label,
        'modalities': mod_names,
        'n_samples': len(meta_aligned),
        'n_features': total_n,
        'n_features_per_modality': dict(zip(mod_names, n_features_list)),
        'n_classes': len(le.classes_),
        'classes': le.classes_.tolist(),
        'macro_f1_mean': f1_mean,
        'macro_f1_std': f1_std,
        'per_fold_f1': fold_f1s,
        'per_fold_C': fold_Cs,
        'classifier': 'LogisticRegression(L1, saga, multinomial)',
        'cv': f'{N_SPLITS}-fold StratifiedKFold(Tissue, source_coverage_check)',
    }
    tag = '_'.join(mod_names)
    out_path = OUT / f'combined_{tag}_4class_edta.json'
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"Saved: {out_path}")
    return result

run_combined_4class(COMBINE_MODALITIES, meta_4c, label='fem4+probe_meth+fragment_length')
