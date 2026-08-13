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

Client requirement: run each methylation version ALONE through the **same Phase 1 protocol** (5-fold source-covering CV stratified by Tissue, StandardScaler → L1 logistic regression with inner GridSearchCV over C = logspace(-3, 1, 6), per-fold mean imputation, macro-F1 primary; Full 164-sample scope + EDTA 96-sample scope). Record per-version macro-F1 + n.

### Step 3A — Grill the protocol FIRST (this is why you were spawned)

Peter flagged two protocol questions that must be grilled (grill-with-docs session, stateful — record resolutions in CONTEXT.md/ADRs) **before any pipeline run**:

1. **Why no tt39 panels?** The tt39 panel (39 probes, `docs/references/tt39-probes.txt`, "previously identified as effective for tissue-of-origin discrimination") is **fully contained in the 148-probe panel (verified: 39/39)** — a tt39-restricted screen variant is free data-wise (subset 148 columns / the per-CpG sites of those 39 probes) for every version and granularity. Open questions to grill: should the screen include a tt39-restricted variant (and for which granularities)? What is tt39's provenance — was it selected on this cohort (→ label-informed, leakage-adjacent like the archived 100kb bins in CONTEXT.md) or independently (→ legitimate fixed-panel prior)? Decide and record.
2. **Who said 20 PCs?** `N_PCS = 20` is a Phase 1 constant ("fixed PCs for high-dim modalities") with no documented rationale. The per-CpG versions (32K unenriched, 54K enriched) go through PCA-20PC → L1. Open questions: is 20 the right PC count for these versions; should the screen vary n_PCs (e.g. 5/10/20/40 or variance-explained), or compare PCA vs per-probe CpG aggregation vs sparse LASSO inside CV? The issue's step-4 candidate list (PCA-20PC, PCA variants, per-probe CpG aggregation, tt39 restriction, sparse LASSO inside CV; RFE-CV ruled out) is the raw material — pulling the DR-route choice into this grilling is in scope.

Natural satellites that hang off these and should also be grilled (evidence in `docs/related-work.md` §c and `docs/methylation-inventory.md`):

- **Imputation at 37% missingness** (probe_cpg enriched): per-fold mean imputation on 37% per-site NaN is the Phase 1 protocol but is suspect at that level (literature: coverage-driven missingness is MNAR). Matched-protocol vs alternative handling — the issue says missingness handling is matched across versions; revisit given the inventory.
- **Success criteria** for "per-CpG adds signal" — agreeing these up front sharpens the screen even though the formal call is step 4's.
- **Unfiltered probe_meth as QC control**: identical sample×probe keys to filtered, read-level difference only — near-identical scores expected; what would a *large* gap mean, and is unfiltered worth running or just informative as a control?

### 3A→3B context management (grilling may outrun the session)

The grilling and the execution may not fit one window. Procedure:

1. **Write-as-you-go during 3A.** After *each* resolution, record it immediately (CONTEXT.md term, ADR, or this handoff) — never accumulate decisions only in the conversation. Each grilling round ends at a checkpoint; the repo is always the complete copy. This is what makes the session disposable without losing the why.
2. **At the 3A→3B boundary, walk the phase-boundary tree (ask-matt):**
   - Enough smart zone left (~150k) and this handoff now states a fully decided, executable protocol (no open questions)? → **Continue** (default — the grilling's reasoning is the primary source for execution).
   - Tight but not spent? → **`/compact`** with the instruction "3A done, decisions recorded, execute 3B" — lossy, but the decisions live in files, so the loss is tolerable.
   - Near the limit? → stop at the nearest round boundary, finish recording, and hand to a **fresh session** via this file (handoff = the portable carrier; the next session needs zero inference). Exit criterion for 3A either way: this handoff is executable from the file alone.
3. **During 3B**, the pipeline runs are long jobs — run them in the background (hub) or via subagent, not inline, so the window isn't consumed by waiting. Compact only at run boundaries.

Never push on degraded context mid-round or mid-run.

### Step 3B — Execution (protocol per grilling; mechanics below stand unless grilling overrides)

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

- Modify `input/`, `metadata`, or `CONTEXT.md` except for terms/decisions resolved in the 3A grilling (record those; domain-modeling discipline).
- Run steps 4–8 (post-screen grilling, deeper exploration, combination check, folded summary) — after the screen, report results and stop; interpretation of results is step 4's grilling with Peter.
- Redo steps 1–2.

## Hand back

Report: grilling decisions made (with rationale), per-version macro-F1 + n table (Full/EDTA, incl. any added variants), any loader/pipeline surprises, and the metrics file paths. Commit scripts + outputs (repo convention: pipeline outputs are committed; the big probe_cpg files are already gitignored).
