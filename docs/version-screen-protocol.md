# Methylation Version Screen Protocol (issue #8, step 3)

Screen: run each methylation feature version ALONE through the Phase 1 protocol,
per the client requirement (260731 follow-up). This document records the
protocol as grilled with the user (2026-08-13). Inputs: `docs/methylation-inventory.md`
(step 1), `docs/related-work.md` (step 2).

## Design decisions (from grilling session)

| Question | Decision |
|----------|----------|
| tt39 provenance | Probe set provided for consideration; no provenance record exists in repo. Treated as leakage-adjacent (flagged lightly in docs only — do not let it overweight decision-making). |
| tt39 variants | Included: both granularities × both captures (4 variants: probe-avg tt39 × {enriched, unenriched}, per-CpG tt39 × {enriched, unenriched}). |
| per-CpG DR route | Protocol control PCA-20PC for both per-CpG versions, PLUS per-CpG → probe aggregation rows (both captures). Raw LASSO and PC-count variants deferred to step 5. Rationale: PCA keeps variance not class signal (PC1 = 34%/52% of variance; PCA-20PC per-CpG is the weakest Phase 1 modality, 0.399); aggregation is the same-data granularity test, cheap, and inventory shows it is missingness-lossless. |
| Aggregation site set | All measured sites per probe (flanking-inclusive; 89/148 probes have site span > 122bp, median 529bp, max 8,975bp). Unweighted mean of observed site betas. NOT equivalent to probe_meth files (read-weighted, 122bp-window, verified: `CpG_frac == CpG_meth/CpG_total` exactly). |
| Flanking-site design question | Open client question (why sites extend beyond the 122bp window; deliberate targets vs capture spillover). Record in exploration log; do not block the screen. Optional evidence check: coverage-by-distance on enriched long file. |
| Imputation | Matched per-fold mean imputation for ALL rows (Phase 1 protocol, no leakage). No site filtering, no coverage thresholding (unenriched has no per-CpG long file → 10x filtering cannot be matched across captures). Missingness reported alongside every score. |
| Unfiltered probe_meth QC | Run as a screen version (Full + EDTA). Inventory corrected the issue's "expected worse": filtered/unfiltered differ at read level only (CpG_frac shift −0.017 mean); expect near-identical classification. The row demonstrates the QC's classification impact empirically. |
| Success criteria for "per-CpG adds signal" | No hard threshold. Report per-fold paired deltas, mean delta, and fold-consistency (folds are identical across rows: seed 42, metadata-only) for both scopes. Step-5 warrant: per-CpG ≥ probe-avg in Full AND not reversed in EDTA. Formal call stays the post-screen grilling (issue protocol: no winner crowned from CV alone; provisional until Phase 2a). |
| Reuse existing outputs | Reuse `probe_meth` (enriched) and `probe_cpg` (unenriched) Phase 1 outputs — produced by current code (pipeline untouched since 4e17307), identical protocol/folds. Step 8 `--combine` rerun is deterministic and covers reproducibility. |
| Variant scope matrix | Full (164) + EDTA (96) for ALL rows, including variants (tt39 ×4, aggregation ×2). |
| tt39 per-CpG DR | BOTH PCA-20PC (protocol) and raw LASSO (direct L1 on ~5–8K tt39 sites — cheap at this size; doubles as the raw-LASSO pilot for step 5). |

## Run matrix

Rows marked (reuse) are existing `output/phase1/` outputs; all others are new runs.

| # | Row | Features | Path | New code |
|---|-----|----------|------|----------|
| 1 | probe_meth enriched filtered | 148 | wide, existing config | — (reuse) |
| 2 | probe_meth unenriched filtered | 148 | LONG, pivot | `probe_meth_unenriched` config, `long_format: True` |
| 3 | probe_meth unfiltered QC | 148 | LONG, pivot | `probe_meth_unfiltered_qc` config, `long_format: True` |
| 4 | probe_cpg unenriched PCA-20PC | 32K→20 PCs | wide, existing config | — (reuse) |
| 5 | probe_cpg enriched PCA-20PC | 54K→20 PCs | wide, new config | `probe_cpg_enriched` config |
| 6 | per-CpG aggregated → probes, unenriched | 148 | derived: mean of observed sites per probe | loader transform |
| 7 | per-CpG aggregated → probes, enriched | 148 | derived: mean of observed sites per probe | loader transform |
| 8 | probe-avg tt39 enriched | 39 | subset of row 1 columns | config |
| 9 | probe-avg tt39 unenriched | 39 | subset of row 2 columns | config |
| 10 | per-CpG tt39 unenriched PCA-20PC | ~5–8K→20 PCs | subset of row 4 columns | config |
| 11 | per-CpG tt39 unenriched raw LASSO | ~5–8K | same subset, no PCA | config + DR flag |
| 12 | per-CpG tt39 enriched PCA-20PC | ~5–8K→20 PCs | subset of row 5 columns | config |
| 13 | per-CpG tt39 enriched raw LASSO | ~5–8K | same subset, no PCA | config + DR flag |

Every row runs Full (164) + EDTA (96) scopes. tt39 per-CpG feature counts are
approximate (39 probes' manifest sites).

## Protocol (unchanged from Phase 1 unless noted)

- 5-fold Source-covering CV (StratifiedKFold by Tissue, seed 42, source-coverage
  verification, up to 100 reshuffles)
- L1-logreg (saga, multinomial, class_weight='balanced', max_iter 5000),
  C tuned via inner GridSearchCV over `np.logspace(-3, 1, 6)`
- Per-fold mean imputation (training statistics only) for non-high-dim;
  mean-impute → scale → PCA-20PC (fit on train) for high-dim rows
- Raw-LASSO rows (11, 13): skip the PCA step, run L1 directly on the subset;
  same imputation and scaling as the non-high-dim path
- Report per row: macro-F1 ± std (mean + CI), balanced accuracy, n samples,
  n features, median C, per-source accuracy, per-version missingness.
  Per-fold paired deltas vs probe-avg row where applicable.
- Headline comparisons use the common sample intersection (all rows join
  164/164 — intersection is full); missingness differences reported, not
  imputed away.

## Client questions to record (exploration log)

- Why do per-CpG site sets extend beyond the 122bp probe window (89/148 probes),
  and why do site sets differ between capture conditions? Deliberate targets
  or capture spillover?

## Outputs

`output/phase1/{tag}_cv_metrics.csv`, `{tag}_hyperparameters.json`,
`{tag}_full_model.joblib`, figures; per-run `phase1_summary.csv` overwrite
stays truncated until step 8 (`--combine`).

## Out of scope (deferred)

- Raw LASSO on full per-CpG versions, PC-count variants, variance-explained
  selection → step 5 deeper exploration, informed by screen results.
- Combination with FEM4 → step 6. Exploration log → step 7. Folded summary → step 8.
