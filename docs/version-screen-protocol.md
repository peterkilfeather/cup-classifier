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
| Unfiltered probe_meth QC | Run as a screen version (all scopes). Inventory corrected the issue's "expected worse": filtered/unfiltered differ at read level only (CpG_frac shift −0.017 mean); expect near-identical classification. The row demonstrates the QC's classification impact empirically. |
| **Evaluation target** | **5-class TOO is the default evaluation** (CUP clinical task; 6-class macro-F1 is inflated by the trivially separable healthyblood class; Phase 2a external data is cancer patients → 5-class frozen models are what Phase 2a needs). |
| Scope structure | 3 scopes: **TOO-Full** (124 samples, 5 classes: colon 22, liver 24, pancreas 19, prostate 38, stomach 21; 3 sources — NIH drops out), **TOO-EDTA** (56 samples, 4 classes — pancreas has 0 EDTA samples; BCT control), **6-class Full** (164, 6 classes; Phase 1 comparability anchor only). |
| Reuse existing outputs | **No reuse — all rows rerun fresh** (user decision; coherent single run history). Existing 6-class outputs stay as historical Phase 1 artifacts. |
| Success criteria for "per-CpG adds signal" | No hard threshold. Report per-fold paired deltas, mean delta, and fold-consistency (folds are identical across rows: seed 42, metadata-only) for the TOO scopes. Step-5 warrant: per-CpG ≥ probe-avg in TOO-Full AND not reversed in TOO-EDTA. **Comparisons pair same-capture rows**: per-CpG unenriched vs probe_meth unenriched; per-CpG enriched vs probe_meth enriched; per-CpG vs per-CpG-aggregated is the clean granularity test (same data). Formal call stays the post-screen grilling (issue protocol: no winner crowned from CV alone; provisional until Phase 2a). |
| tt39 per-CpG DR | BOTH PCA-20PC (protocol) and raw LASSO (direct L1 on ~5–8K tt39 sites — cheap at this size; doubles as the raw-LASSO pilot for step 5). |

## Run matrix

13 rows × 3 scopes (TOO-Full, TOO-EDTA, 6-class Full) = 39 runs, all fresh.

| # | Row | Features | Notes |
|---|-----|----------|-------|
| 1 | probe_meth enriched filtered | 148 | wide, existing config |
| 2 | probe_meth unenriched filtered | 148 | LONG, pivot (`probe_meth_unenriched` config, `long_format: True`) |
| 3 | probe_meth unfiltered QC | 148 | LONG, pivot (`probe_meth_unfiltered_qc` config, `long_format: True`) |
| 4 | probe_cpg unenriched PCA-20PC | 32K→20 PCs | wide, existing config |
| 5 | probe_cpg enriched PCA-20PC | 54K→20 PCs | wide, new config `probe_cpg_enriched` |
| 6 | per-CpG aggregated → probes, unenriched | 148 | mean of observed sites per probe (all sites) |
| 7 | per-CpG aggregated → probes, enriched | 148 | mean of observed sites per probe (all sites) |
| 8 | probe-avg tt39 enriched | 39 | subset of row 1 columns |
| 9 | probe-avg tt39 unenriched | 39 | subset of row 2 columns |
| 10 | per-CpG tt39 unenriched PCA-20PC | ~5–8K→20 PCs | subset of row 4 columns |
| 11 | per-CpG tt39 unenriched raw LASSO | ~5–8K | same subset, no PCA |
| 12 | per-CpG tt39 enriched PCA-20PC | ~5–8K→20 PCs | subset of row 5 columns |
| 13 | per-CpG tt39 enriched raw LASSO | ~5–8K | same subset, no PCA |

tt39 per-CpG feature counts are approximate (39 probes' manifest sites).

## Scope definitions (pipeline change)

- Add a `TOO` scope: `meta[meta['Tissue'] != 'healthyblood']` (124), label prefix `TOO `.
- `TOO-EDTA`: `meta[(meta['BCT'] == 'EDTA') & (meta['Tissue'] != 'healthyblood')]` (56), label prefix `TOO EDTA `.
- 6-class Full: existing Full scope (164), no prefix.
- `build_summary` scope inference must handle the new prefixes.
- Chance baselines: TOO-Full 1/5, TOO-EDTA 1/4, 6-class Full 1/6.
- Source-covering CV applies as-is (sources present in each scope derived from data).

## Protocol (unchanged from Phase 1 unless noted)

- 5-fold Source-covering CV (StratifiedKFold by Tissue, seed 42, source-coverage
  verification, up to 100 reshuffles)
- L1-logreg (saga, multinomial, class_weight='balanced', max_iter 5000),
  C tuned via inner GridSearchCV over `np.logspace(-3, 1, 6)`
- Per-fold mean imputation (training statistics only) for non-high-dim;
  mean-impute → scale → PCA-20PC (fit on train) for high-dim rows
- Raw-LASSO rows (11, 13): skip the PCA step, run L1 directly on the subset;
  same imputation and scaling as the non-high-dim path
- Report per row: scope, macro-F1 ± std (mean + CI), balanced accuracy,
  n samples, n features, median C, per-source accuracy, per-version missingness.
  Per-fold paired deltas vs the same-capture probe-avg row.
- Headline comparisons use the common sample intersection (all rows join the
  full universe — intersection is full at each scope); missingness differences
  reported, not imputed away.

## Client questions to record (exploration log)

- Why do per-CpG site sets extend beyond the 122bp probe window (89/148 probes),
  and why do site sets differ between capture conditions? Deliberate targets
  or capture spillover?

## Outputs

`output/phase1/{tag}_cv_metrics.csv`, `{tag}_hyperparameters.json`,
`{tag}_full_model.joblib`, figures; per-run `phase1_summary.csv` overwrite
stays truncated until step 8. Frozen full models are now 5-class TOO models
for TOO scopes (Phase 2a relevant), 6-class for the comparability scope.

## Step 8 note

The issue's "expect 12 rows" (10 modality×scope + 2 combined) is obsolete:
the screen produces 51 rows (39 + 12 addendum rows, issue #18) + combined
rows. Step 8 must be rewritten with the new scope structure and row counts.

## Out of scope (deferred)

- Raw LASSO on full per-CpG versions, PC-count variants, variance-explained
  selection → step 5 deeper exploration, informed by screen results.
- Combination with FEM4 → step 6 (needs 5-class FEM4 rerun for comparability).
- Exploration log → step 7. Folded summary → step 8 (rewritten).

## Addendum: in-probe-only rows (issue #18) — mechanics landed, execution pending

Extends the screen with rows restricted to the **In-Probe Site Set** (sites
within the 122 bp probe window) and a concrete decomposition of the
difference between client probe-averaged values (`probe_meth`: read-weighted
`CpG_meth / CpG_total` over the window) and the per-CpG aggregation
(unweighted mean of observed site betas over the flanking-inclusive set).
Client-committed (2026-09-01 hand-back): "I will be running a parallel
version that excludes those CpGs we talked about that are located beyond
the probe bounds."

### Window definition

- Center: probe mapinfo from `probe_annotations_450k.csv` (the 450K target).
- Half-width: config key `half_width`, default **61 bp** (122 bp region
  centered on the target).
- Membership derived at runtime: manifest chr/pos joined with annotation
  chr/mapinfo (`|pos − mapinfo| ≤ half_width`), chr naming normalized
  (manifest `chr1`, annotations `1`). Nothing probe- or count-specific is
  hardcoded; per-probe site counts are computed at load.
- Sanity anchor: the assumed ±61 bp window is checked against the client's
  pipeline scripts + one-sample raw Bismark output (requested 2026-09-07;
  reply accepted 2026-09-07) before the execution runs.
- Under the assumed window (verified 2026-09-08): 3,546/32,084 sites
  in-probe (unenriched), 5,865/54,300 (enriched); 147/148 probes have ≥1
  in-probe site — cg14861089 has none (derived; report the actual counts
  per run). The tt39 probe set sits fully in-window (max site distance
  61 bp, both captures; guard test re-checks under any window change).

### Decisions (grilled 2026-09-08)

- **Per-CpG in-probe rows use raw LASSO, NOT PCA-20PC.** In-window filtering
  cuts features ~9x (32K→3.5K / 54K→5.9K), which brings raw LASSO back on
  the table (the tt39 LASSO pilot already outperformed PCA in TOO-EDTA).
  PCA runs are deferred until PC-count appropriateness is examined.
- **No winner crowned** from this addendum; the step-4 grilling session
  (winner call, issue #8) stays open and now gates on these rows (the
  per-CpG deltas from #12 are 89% flanking sites — the call could flip
  in-window).
- 12 fresh runs (4 rows × 3 scopes); no reuse. Execution waits on the client
  window definition — mechanics land first; the ticket stays open.

### Run matrix (addendum rows 14–17)

| # | Row | Features | Notes |
|---|-----|----------|-------|
| 14 | per-CpG in-probe unenriched raw LASSO | ~3.5K | wide subset, `in_probe_only`, `dr: lasso` |
| 15 | per-CpG in-probe enriched raw LASSO | ~5.9K | wide subset, `in_probe_only`, `dr: lasso` |
| 16 | per-CpG aggregated in-probe → probes, unenriched | 147 | mean of in-probe sites; cg14861089 drops |
| 17 | per-CpG aggregated in-probe → probes, enriched | 147 | mean of in-probe sites; cg14861089 drops |

Feature counts are pre-drop memberships (all-NaN feature drops apply at
load); actual per-run counts are reported from the runs. Probes with zero
in-probe sites drop out as a consequence of the filter — never hardcoded.

### Reporting (paired deltas per scope)

- agg in-probe vs probe_meth — window-matched; isolates **weighting**
  (read-weighted vs unweighted)
- agg in-probe vs agg flanking-inclusive — isolates **site set** (window vs
  flanking-inclusive)
- per-CpG in-probe LASSO vs agg in-probe — same-data granularity, in-window

### Decomposition deliverable (execution)

Per-probe delta table of client probe_meth vs our aggregation, split by
mechanism:

- **weighting** (read-weighted vs unweighted): enriched full 2×2 —
  read-weighted in-window (= client), read-weighted all-sites (from the
  enriched per-CpG long file), unweighted in-window (= agg in-probe),
  unweighted all-sites (= agg). Unenriched: 3 of 4 cells (no per-CpG long
  file).
- **site set** (window vs flanking-inclusive)
- **missingness** (unweighted mean over few observed sites; cg23960736 is
  0% flanking yet the largest delta +0.216)

Headline measured (2026-09-08): median per-probe r 0.990 (enriched) / 0.970
(unenriched); MAD ≈ 0.024; 19/148 probes with |mean Δ| > 0.05; worst
cg23960736 +0.216, cg14743594 +0.213 (91% flanking), cg14394692 +0.184 (95%
flanking). Reconcile with the client scripts when received; verify the
one-sample raw Bismark output reproduces their `CpG_frac`.
