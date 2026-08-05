"""Controlled ablation: pin CV splits once, test all combos on same folds.
FIX: align metadata to FEM4 sample ordering before split creation."""
import sys, warnings, json
from pathlib import Path
import numpy as np
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

# ── Preload all modalities ──
ALL_MODS = ['fem4', 'probe_meth', 'fragment_length']
raw_data = {}
for mod in ALL_MODS:
    cfg = get_modality_cfg(mod)
    X, ids, feat_names = load_modality(cfg, meta_4c['TWIST_ID'].values, impute=False)
    raw_data[mod] = {'X': X, 'ids': ids, 'feat_names': feat_names}

# Align all modalities to fem4's sample ordering (reference)
ref_ids = raw_data['fem4']['ids']
for mod in ALL_MODS:
    ids = raw_data[mod]['ids']
    if not np.array_equal(ids, ref_ids):
        order = {s: i for i, s in enumerate(ref_ids)}
        idx = np.array([order[s] for s in ids])
        raw_data[mod]['X'] = raw_data[mod]['X'][idx]
        raw_data[mod]['ids'] = ref_ids

# ── Align metadata to FEM4 order (CRITICAL: matches original test_edta_4class.py) ──
meta_aligned = meta_4c.set_index('TWIST_ID').loc[ref_ids].reset_index()
print(f"EDTA 4-class: {len(meta_aligned)} samples, {meta_aligned['Tissue'].nunique()} tissues")

# ── Pin CV splits once, using FEM4-aligned metadata ──
y = meta_aligned['Tissue'].values
sources = meta_aligned['Source'].values
required_sources = set(sources)
for attempt in range(MAX_RESHUF):
    seed = RANDOM_STATE + attempt
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
    splits = list(skf.split(np.zeros(len(meta_aligned)), y))
    if all(required_sources.issubset(set(sources[tr_idx])) for tr_idx, _ in splits):
        print(f"  Split found at seed={seed} (attempt {attempt})")
        break

print("\n  Fold compositions (test set):")
for fold, (_, te_idx) in enumerate(splits):
    comp = meta_aligned.iloc[te_idx].groupby(['Tissue','Source']).size()
    print(f"    Fold {fold+1}: {dict(comp)}")

le = LabelEncoder()
y_enc = le.fit_transform(meta_aligned['Tissue'].values)

def evaluate(combo_name, mod_names):
    """Run CV on pinned splits for given modality subset."""
    fold_f1s = []
    for tr_idx, te_idx in splits:
        fold_Xs_tr, fold_Xs_te = [], []
        for mod in mod_names:
            Xm = raw_data[mod]['X']
            Xm_tr, Xm_te = Xm[tr_idx], Xm[te_idx]
            col_mean = np.nanmean(Xm_tr, axis=0)
            col_mean = np.nan_to_num(col_mean, nan=0.0)
            for X_src in [Xm_tr, Xm_te]:
                inds = np.where(np.isnan(X_src))
                X_src[inds] = np.take(col_mean, inds[1])
            scaler = StandardScaler()
            fold_Xs_tr.append(scaler.fit_transform(Xm_tr))
            fold_Xs_te.append(scaler.transform(Xm_te))
        X_tr = np.hstack(fold_Xs_tr)
        X_te = np.hstack(fold_Xs_te)
        y_tr = y_enc[tr_idx]
        y_te = y_enc[te_idx]

        n_inner = min(3, np.min(np.bincount(y_tr)))
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=RANDOM_STATE)
        clf = LogisticRegression(l1_ratio=1, solver='saga', class_weight='balanced',
                                 max_iter=MAX_ITER, random_state=RANDOM_STATE)
        gs = GridSearchCV(clf, {'C': C_GRID}, cv=inner_cv, scoring='f1_macro', n_jobs=1)
        gs.fit(X_tr, y_tr)
        model = gs.best_estimator_
        y_pred = model.predict(X_te)
        fold_f1s.append(f1_score(y_te, y_pred, average='macro'))
    return np.mean(fold_f1s), np.std(fold_f1s), fold_f1s

# ── Run all combos on pinned splits ──
combos = [
    ('FEM4', ['fem4']),
    ('Fragment', ['fragment_length']),
    ('Probe_meth', ['probe_meth']),
    ('FEM4+Fragment', ['fem4', 'fragment_length']),
    ('FEM4+Probe_meth', ['fem4', 'probe_meth']),
    ('Fragment+Probe_meth', ['fragment_length', 'probe_meth']),
    ('All 3', ['fem4', 'probe_meth', 'fragment_length']),
]

print(f"\n{'Combo':25s} {'F1 mean':>10s} {'F1 std':>8s}  {'Per-fold F1s':>30s}")
print('-' * 75)
results = {}
for name, mods in combos:
    m, s, fold_f1s = evaluate(name, mods)
    results[name] = {'mean': m, 'std': s, 'per_fold': fold_f1s}
    f1s_str = ' '.join(f'{f:.3f}' for f in fold_f1s)
    print(f"{name:25s} {m:10.4f} {s:8.4f}  {f1s_str}")

out_path = OUT / 'controlled_ablation_fixed.json'
with open(out_path, 'w') as f:
    json.dump({'splits_seed': seed, 'fold_compositions': [
        dict(meta_aligned.iloc[te_idx].groupby(['Tissue','Source']).size())
        for _, te_idx in splits
    ], 'results': results}, f, indent=2)
print(f"\nSaved: {out_path}")
