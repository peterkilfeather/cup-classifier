# Handoff: Methylation version screen execution (issue #8, step 3)

## Repo

`github.com/peterkilfeather/cup-classifier` at `/xscratch/farney/cup-classifier`, branch `main`.
Issue: https://github.com/peterkilfeather/cup-classifier/issues/8.

## Grilling is DONE — decisions are recorded, do not re-grill

The protocol was grilled with the user 2026-08-13 and confirmed. The executable
spec is **`docs/version-screen-protocol.md`** — read it first; it overrides this
file on any conflict. Key decisions (full table in the protocol doc):

- **5-class TOO is the default evaluation.** 3 scopes: TOO-Full (124, 5 classes),
  TOO-EDTA (56, 4 classes — no pancreas), 6-class Full (164, comparability anchor).
- **13 rows × 3 scopes = 39 runs, all fresh** (no reuse of existing outputs).
- Per-CpG DR: PCA-20PC (protocol) + aggregation-to-probes rows; tt39 per-CpG rows
  additionally run raw LASSO (the pilot).
- Matched per-fold mean imputation everywhere; missingness reported, not imputed away.
- Success criteria: per-fold paired deltas + mean delta + fold-consistency vs the
  **same-capture** probe-avg row; step-5 warrant = per-CpG ≥ probe-avg in TOO-Full
  and not reversed in TOO-EDTA. No winner crowned; formal call stays post-screen grilling.
- tt39 treated as leakage-adjacent (flag lightly, don't overweight).

## Change list (mechanics)

1. `scripts/data_loading.py` — add configs to `MODALITY_CONFIGS` (existing pattern):
   - `probe_meth_unenriched` (LONG file `probe_meth/all_samples.probe_meth_filtered.mNonCpGlt4.long_v2.tsv`, label 'Methylation (probe-avg, unenriched)')
   - `probe_meth_unfiltered_qc` (LONG file `all_samples.probe_meth_unfiltered.long_v2.tsv`, label 'Methylation (probe-avg, unfiltered QC)')
   - `probe_cpg_enriched` (wide 54K, `high_dim: True`, label 'Per-CpG methylation (enriched)')
   - tt39 subsets: `probe_meth_tt39_enriched`, `probe_meth_tt39_unenriched` (39 probe columns from `docs/references/tt39-probes.txt`), `probe_cpg_tt39_unenriched`, `probe_cpg_tt39_enriched` (manifest sites of the 39 probes)
   - aggregation rows: `probe_cpg_agg_unenriched`, `probe_cpg_agg_enriched` (mean of observed site betas per probe — groupby on manifest `probe_id`; 148 features, NOT high_dim)
   - `long_format: True` on the two LONG configs: in `load_modality`, when set → **raise ValueError on duplicate sample×probe rows**, then pivot to wide (index `sample`, columns `probe_id`, values `CpG_frac`) before the standard path.
   - raw-LASSO variants: a DR flag (e.g. `dr: 'lasso'`) that skips the PCA step in the high-dim path (mean-impute + scale + L1 directly on the ~5–8K tt39 sites).
2. `scripts/run_phase1_pipeline.py`:
   - Add TOO scopes to `main()`: TOO = `meta[meta['Tissue'] != 'healthyblood']` (label prefix `TOO `), TOO-EDTA = `meta[(meta['BCT'] == 'EDTA') & (meta['Tissue'] != 'healthyblood')]` (prefix `TOO EDTA `). Keep 6-class Full (no prefix).
   - `build_summary`: extend scope inference — check `TOO ` prefix BEFORE `EDTA ` (note: 'TOO EDTA ' also contains 'EDTA').
   - Add all 13 row names to `PHASE1_MODALITIES` (CLI choices derive from it).
3. Run: `python scripts/run_phase1_pipeline.py --modalities <all 13>` per scope family (or add a `--scopes` flag). Long jobs — **background them (hub) or use a subagent; never burn the window waiting**. Heavy rows: per-CpG PCA (54K) × 3 scopes, tt39 per-CpG × 2 DR × 3 scopes. Est. total 4–8h.

## Report (hand back)

Per-scope table: row, macro-F1 ± std, n samples, n features, median C, missingness.
Per-fold paired deltas vs same-capture probe-avg row (TOO scopes). Loader/pipeline
surprises. Metrics file paths. Do NOT crown a winner.

## Don't

- Modify `input/`, `metadata/`. Don't re-grill decisions recorded in `docs/version-screen-protocol.md`.
- Run steps 4–8 (deeper exploration, combination, log, folded summary). Step 8's
  "12 rows" is obsolete — the issue text needs updating when step 8 is scheduled.
- Delete `docs/version-screen-protocol.md` or this handoff without the user's instruction.
- Commit per-run; commit scripts + outputs at the end (pipeline outputs are committed by convention; large probe_cpg files are gitignored).
