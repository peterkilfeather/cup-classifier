"""Check fold composition for the combined 4-class EDTA runs."""
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings('ignore')
np.random.seed(42)

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / 'scripts'))
from data_loading import load_metadata, get_modality_cfg, load_modality

N_SPLITS = 5
MAX_RESHUF = 100
RANDOM_STATE = 42

meta = load_metadata()
meta_edta = meta[meta['BCT'] == 'EDTA'].copy()
meta_4c = meta_edta[meta_edta['Tissue'] != 'healthyblood'].copy()

COMBINE_MODALITIES = ['fem4', 'probe_meth', 'fragment_length']

def get_fold_composition(meta, seed_offset=0):
    """Return per-fold test-set Tissue × Source breakdown for one CV split."""
    y = meta['Tissue'].values
    sources = meta['Source'].values
    required_sources = set(sources)
    for attempt in range(MAX_RESHUF):
        seed = RANDOM_STATE + seed_offset + attempt
        skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        splits = list(skf.split(np.zeros(len(meta)), y))
        if all(required_sources.issubset(set(sources[tr_idx])) for tr_idx, _ in splits):
            break
    fold_data = []
    for fold, (_, te_idx) in enumerate(splits):
        te_meta = meta.iloc[te_idx]
        comp = te_meta.groupby(['Tissue', 'Source']).size().to_dict()
        fold_data.append((fold, comp, len(te_idx)))
    return fold_data

# For the all-3 combined run (seed=42, retries until source coverage satisfied)
meta_aligned = meta_4c  # all samples pass through
print("=== Combined 3-way (FEM4+probe_meth+fragment_length) fold composition ===")
folds = get_fold_composition(meta_4c, seed_offset=0)
for fold, comp, n in folds:
    print(f"  Fold {fold+1} (n_test={n}):")
    for (t, s), cnt in sorted(comp.items()):
        print(f"    {t:12s} / {s:25s}: {cnt}")

# The ablation runs each have their own seed (depends on retry count)
# Let me also check what seed the FEM4+Fragment run landed on
# Actually, let's just check multiple seed offsets to see variation
print("\n=== Fold composition variability across seeds ===")
for off in range(3):
    folds = get_fold_composition(meta_4c, seed_offset=off)
    print(f"\nSeed offset {off}:")
    for fold, comp, n in folds:
        min_class = min(comp.values())
        tissues = list(set(t for (t,s) in comp.keys()))
        print(f"  Fold {fold+1}: {n} test, {len(tissues)} tissues, min_class={min_class}, "
              f"classes: {sorted(comp.keys())}")
