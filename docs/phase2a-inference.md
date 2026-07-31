# Phase 2a: External Validation Inference

**Goal:** Apply frozen Phase 1 models to ~41 incoming samples and evaluate
whether the tissue classifier generalizes to new data.

**Success criterion:** macro-F1 on incoming data ≥ 0.7 × Phase 1 CV macro-F1
AND > 0.17 (6-class chance).

## Frozen model structure

Models are at `output/phase1/models/{scope}_{modality}_full_model.joblib`.
Each is a dict with these keys depending on modality type:

### Low-dim models (fem4, probe_meth, fragment_length)

```python
{
    'model': LogisticRegression,         # fitted L1-logreg
    'scaler': StandardScaler,            # fitted on all Phase 1 data
    'C': float,                          # median C across CV folds
    'classes': list,                     # class labels in order
    'feature_names': ndarray,            # original feature names
    'selected_features': ndarray,        # names of non-zero coef features
    'selected_mask': ndarray,            # bool mask of non-zero coefs
    'is_high_dim': False,
    'modality': str,
    'modality_label': str,
    'n_samples': int,
    'cv_macro_f1_mean': float,
    'cv_macro_f1_std': float,
    'timestamp': str,
}
```

### High-dim models (probe_cpg, end_density)

```python
{
    'model': LogisticRegression,         # fitted L1-logreg on PCs
    'pca_scaler': StandardScaler,        # fitted on all raw data
    'pca': PCA,                          # fitted on all scaled data
    'n_pcs': int,                        # number of PCs (≤20)
    'feature_names': ndarray,            # original raw feature names
    'pc_names': list,                    # PC1..PCN labels
    'selected_features': list,           # PC names with non-zero coefs
    'is_high_dim': True,
    # ... same metadata keys as low-dim
}
```

## Inference pipeline

### Low-dim modality

```python
import joblib
import numpy as np

model_dict = joblib.load('output/phase1/models/fem4_full_model.joblib')
model = model_dict['model']
scaler = model_dict['scaler']
classes = model_dict['classes']

# X_new: (n_incoming, n_features) numpy array, same features as Phase 1
X_scaled = scaler.transform(X_new)
y_pred = model.predict(X_scaled)
y_prob = model.predict_proba(X_scaled)
```

### High-dim modality

```python
model_dict = joblib.load('output/phase1/models/end_density_full_model.joblib')
model = model_dict['model']
pca_scaler = model_dict['pca_scaler']
pca = model_dict['pca']
classes = model_dict['classes']

# X_new: (n_incoming, n_original_features) — must match Phase 1 feature set
X_imp = np.nan_to_num(X_new, nan=0.0)   # crude NaN handling
X_scaled = pca_scaler.transform(X_imp)
X_pc = pca.transform(X_scaled)
y_pred = model.predict(X_pc)
y_prob = model.predict_proba(X_pc)
```

**Note:** High-dim NaN imputation above is crude (replace with 0). If incoming
data has extensive missingness, impute using Phase 1 training column means
(stored implicitly in the PCA scaler's `mean_` attribute isn't directly
available; reconstruct from the reference model's training data or use
modality-appropriate imputation).

## Evaluation

```python
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score

y_true = ...  # tissue labels for incoming samples
macro_f1 = f1_score(y_true, y_pred, average='macro')
bal_acc = balanced_accuracy_score(y_true, y_pred)
acc = accuracy_score(y_true, y_pred)

# Per-tissue breakdown (required per ticket #6)
for tissue in set(y_true):
    mask = y_true == tissue
    f1 = f1_score(y_true[mask], y_pred[mask], average='macro')
    # or accuracy for single-class per tissue
    acc_tissue = accuracy_score(y_true[mask], y_pred[mask])
```

Report per-source if source labels are known for incoming data.

## Per-modality inference script

A convenience script `scripts/run_phase2a_inference.py` does not yet exist.
The recommended structure:

1. Accept `--model` (path to joblib), `--input` (feature CSV), `--metadata` (labels)
2. Detect `is_high_dim` from model dict
3. Run the appropriate inference path
4. Print macro-F1, balanced accuracy, per-tissue accuracy
5. Compare against the Phase 1 CV performance from `model_dict['cv_macro_f1_mean']`
