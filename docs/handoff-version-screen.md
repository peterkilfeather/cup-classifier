# Handoff: Methylation version screen (issue #8, step 3)

## Repo

`github.com/peterkilfeather/cup-classifier` at `/xscratch/farney/cup-classifier`, branch `main`, HEAD `75defc3`.
Issue: https://github.com/peterkilfeather/cup-classifier/issues/8 (steps 1–2 committed; this handoff covers **step 3 only**).

## Already done (steps 1–2, committed — read, do not redo)

- `docs/methylation-inventory.md` — per-version stats (all 5 versions join 164/164 to `metadata_cleaned.csv`; probe_meth missingness ≤0.7%, probe_cpg unenriched 6.0%, enriched 37.1%; long-file layout; overlap table)
- `docs/related-work.md` — 4 cited literature sections (granularity, capture, missingness/imputation, aggregation)
- `docs/figures/methylation_01..06_*.png` — client-facing figures, regenerable via `scripts/make_methylation_figures.py`
- `CONTEXT.md` — vocabulary incl. Per-Sample/Per-Feature Missingness, Coverage; use its language throughout

## Task — version screen (issue #8 step 3)

Client requirement: run each methylation version ALONE through the **same Phase 1 protocol** (5-fold source-covering CV stratified by Tissue, StandardScaler → L1 logistic regression with inner GridSearchCV over C = logspace(-3, 1, 6), per-fold mean imputation, macro-F1 primary; Full 164-sample scope + EDTA 96-sample scope). Record per-version macro-F1 + n. Report results; **do not** interpret/decide (that is step 4's grilling session with Peter).

### Add 3 versions to `scripts/data_loading.py` `MODALITY_CONFIGS`

| Config key | File (input/methylation/) | Format | Label | high_dim |
|---|---|---|---|---|
| `probe_meth_unenriched` | `probe_meth/all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv` | LONG | `Methylation (probe-avg, unenriched)` | False |
| `probe_meth_unfiltered_qc` | `probe_meth/all_samples.probe_meth_unfiltered.long_v2.tsv` | LONG | `Methylation (probe-avg, unfiltered QC)` | False |
| `probe_cpg_enriched` | `probe_cpg/all_samples.probe_cpg_enriched_filtered.mNonCpGlt4_frac.wide.tsv` | wide | `Per-CpG methylation (enriched)` | True |

Follow the existing config pattern (`sep`, `sample_col`, `drop_cols`, `high_dim`, `label`). Existing keys `probe_meth` and `probe_cpg` (unenriched) stay unchanged.

### Loader change (data_loading.py)

- Add `long_format: True` to the two long configs.
- In `load_modality`, when `cfg.get('long_format')`: **raise ValueError on duplicate sample×probe rows** (verified: 0 duplicates in both files — 29,297 rows each), then pivot to wide (index `sample`, columns `probe_id`, values `CpG_frac`) **before** the standard path. Missing (sample, probe) pairs become NaN, matching the inventory's missingness definition.
- Wide path unchanged (note: existing silent `keep='first'` sample dedup stays as-is).

### CLI change (run_phase1_pipeline.py)

- Add the 3 names to `PHASE1_MODALITIES` (argparse `--modalities` choices derive from it). Default run then covers 8 modalities; for step 3 run explicitly:

```
python scripts/run_phase1_pipeline.py --modalities probe_meth_unenriched probe_meth_unfiltered_qc probe_cpg_enriched
```

Long-running job (probe_cpg enriched is 54,300 features, high_dim → PCA-20PC path; per-fold imputation) — run in background/with `hub`. Both Full and EDTA scopes run automatically. Do NOT pass `--combine` (that is step 8).

### Outputs / acceptance

- `output/phase1/{tag}_cv_metrics.csv` per version (tag derives from label), `{tag}_full_model.joblib`, `{tag}_cv_results.png`, `{tag}_coefficients.png`; `phase1_summary.csv` is overwritten (currently 2 stale rows — step 8 regenerates the final 12-row version; don't chase it now).
- Table of per-version macro-F1 (Full + EDTA) with n alongside every score (n = 164/164 for all versions after metadata join — inventory table is the control).
- Protocol controls (decided in the issue — apply without re-deciding): compare on the common sample intersection (trivially all 164 here); do NOT crown a winner from CV alone (provisional until Phase 2a); report missingness differences, don't impute them away.
- Sanity: probe_cpg enriched has 1,156 all-NaN features (of 54,300) — the loader's standard all-NaN drop removes them; per-sample missingness 37% goes through per-fold mean imputation (same as Phase 1 protocol). Expected to underperform; the grilling interprets.
- Unfiltered vs filtered probe_meth: identical sample×probe keys, read-level difference only (CpG_frac shifted −0.017 mean) — near-identical scores expected; treat as the QC control showing the filter's effect, per the issue.

## Don't

- Modify `input/`, `metadata`, or `CONTEXT.md`.
- Run steps 4–8 (grilling, deeper exploration, combination check, folded summary) — stop after recording the screen results.
- Redo steps 1–2.

## Hand back

Report: per-version macro-F1 + n table (Full/EDTA), any loader/pipeline surprises, and the metrics file paths. Commit scripts + outputs (repo convention: pipeline outputs are committed; the big probe_cpg files are already gitignored).
